#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compartment: typed replacement for SimCell's per-layer dict representation.

Each SimCell layer (slab, solvent, vacuum) is a Compartment subclass pairing
its construction-time parameters with a build() method that produces (or, for
vacuum, bypasses) the layer's MDAnalysis Universe. SimCell._stack_layers
dispatches on the concrete subclass instead of a `layer["type"]` string.
"""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, List, Optional

import MDAnalysis as mda

from mdinterface.build.regions import FilledRegion
from mdinterface.build.solvent import make_solvent_box

logger = logging.getLogger(__name__)


class Compartment(ABC):
    """One layer of a SimCell."""

    @abstractmethod
    def build(self, xsize: float, ysize: float, all_sp_univs: list,
              layered: bool = False, do_match: bool = True) -> Optional[mda.Universe]:
        """Return this layer's assembled Universe, or None."""
        raise NotImplementedError


@dataclass
class SlabCompartment(Compartment):
    """A solid-interface layer (slab or pre-built), tiled during SimCell._fit_slabs."""
    species: Any
    nlayers: int
    label: str = "slab"
    _slab: Any = field(default=None, repr=False, compare=False)
    _native_xy: Optional[tuple] = field(default=None, repr=False, compare=False)

    def build(self, xsize, ysize, all_sp_univs, layered=False, do_match=True):
        if self._slab is None:
            return None
        return self._slab.to_universe(
            layered=layered, match_cell=do_match, xydim=[xsize, ysize]
        )


@dataclass
class VacuumCompartment(Compartment):
    """An empty vacuum gap; never packed via PACKMOL."""
    zdim: float = 0.0

    def build(self, xsize, ysize, all_sp_univs, layered=False, do_match=True):
        return None


@dataclass
class SolventCompartment(Compartment):
    """A solvent (liquid) layer, optionally with dissolved species and nested regions."""
    solvent: List[Any]
    solute: List[Any]
    nsolute: Any
    zdim: float
    density: Optional[float]
    nsolvent: Any
    concentration: Optional[float]
    conmodel: Optional[dict]
    solute_pos: Any
    dilate: float
    packmol_tolerance: float
    ratio: Optional[List[float]]
    regions: List[FilledRegion] = field(default_factory=list)
    seed: Optional[int] = None

    def build(self, xsize, ysize, all_sp_univs, layered=False, do_match=True):
        eff_zdim = self.zdim * self.dilate
        eff_rho = self.density / self.dilate if self.density is not None else None

        if self.dilate != 1.0:
            logger.info(
                "  >> dilation x%.2f: packing %.1f Å at %.3f g/cm³  (target: %.1f Å)",
                self.dilate, eff_zdim,
                eff_rho if eff_rho is not None else float("nan"),
                self.zdim,
            )

        return make_solvent_box(
            species=all_sp_univs,
            solvent=self.solvent or None,
            solute=self.solute or None,
            volume=[xsize, ysize, eff_zdim],
            density=eff_rho,
            nsolute=self.nsolute,
            concentration=self.concentration,
            conmodel=self.conmodel,
            solute_pos=self.solute_pos,
            nsolvent=self.nsolvent,
            tolerance=self.packmol_tolerance,
            ratio=self.ratio,
            regions=self.regions,
            seed=self.seed,
        )
