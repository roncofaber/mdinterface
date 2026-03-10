#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
solvent.py — solvent box assembly.

Contains parameter validation, molecule-count computation, and the
high-level ``make_solvent_box`` entry point.  Low-level PACKMOL
orchestration lives in ``box.py`` (``populate_box``).
"""

from typing import Any, Dict, List, Optional, Tuple, Union

import warnings

import numpy as np
import MDAnalysis as mda
from ase import units

from mdinterface.build.box import populate_box
from mdinterface.build.continuum2sim import discretize_concentration

import logging

logger = logging.getLogger("mdinterface.box")

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _validate_solvent_box_parameters(
    nions: Optional[Union[int, List[int]]],
    concentration: Optional[float],
    conmodel: Optional[Dict[int, Tuple[List[float], List[float]]]],
    ions: Optional[List[Any]],
    solvents: List[Any],
    density: Optional[float],
    nsolvent=None,
    ratio: Optional[List[float]] = None,
) -> None:
    """Validate parameter combinations for make_solvent_box."""

    if nions is not None and concentration is not None:
        raise ValueError("Cannot specify both 'nions' and 'concentration'. Use one or the other.")

    if conmodel is not None and (ions is None or len(ions) == 0):
        raise ValueError("When using 'conmodel', 'ions' must be provided and non-empty.")

    if isinstance(nions, (list, tuple)) and ions is not None:
        if len(nions) != len(ions):
            raise ValueError(f"Length of 'nions' ({len(nions)}) must match number of ion species ({len(ions)}).")

    if ratio is not None:
        if len(ratio) != len(solvents):
            raise ValueError(f"Length of 'ratio' ({len(ratio)}) must match number of solvent species ({len(solvents)}).")
        if isinstance(nsolvent, (list, tuple)):
            raise ValueError("Cannot use 'ratio' with 'nsolvent' as a list; use ratio+density or ratio+nsolvent(int).")
        if density is None and nsolvent is None:
            raise ValueError("'ratio' requires either 'density' or 'nsolvent' (total count).")

    if isinstance(nsolvent, (list, tuple)) and solvents:
        if len(nsolvent) != len(solvents):
            raise ValueError(f"Length of 'nsolvent' ({len(nsolvent)}) must match number of solvent species ({len(solvents)}).")

    if len(solvents) > 1 and ratio is None and not isinstance(nsolvent, (list, tuple)) and density is not None and nsolvent is None:
        raise ValueError("For a solvent mixture, specify 'nsolvent' (list), 'ratio'+'density', or 'ratio'+'nsolvent'.")

    if density is not None and nsolvent is not None and ratio is None:
        warnings.warn(
            "Both 'density' and 'nsolvent' are specified. Using 'nsolvent' and ignoring 'density'.",
            UserWarning, stacklevel=4,
        )

    if (density is not None or nsolvent is not None) and not solvents:
        warnings.warn("Density or nsolvent specified but no solvent provided. Will be ignored.",
                      UserWarning, stacklevel=4)

    if not solvents and (ions is None or len(ions) == 0):
        warnings.warn("No solvent or ions specified. Empty box will be created.",
                      UserWarning, stacklevel=3)


# ---------------------------------------------------------------------------
# Ion placement
# ---------------------------------------------------------------------------

def populate_with_ions(ions, nions, volume, ion_pos=False, conmodel=None):

    def place_ion(ion, volume, ion_coords, ion_radii, zpos=None, max_attempts=100):
        ion_radius = ion.estimate_specie_radius()
        for _ in range(max_attempts):
            new_coord = ion_radius + 1 + np.random.rand(3) * (volume - 2 * (ion_radius + 1))
            if zpos is not None:
                new_coord[2] = zpos
            if not ion_coords or np.all(np.linalg.norm(ion_coords - new_coord, axis=1) >= np.array(ion_radii) + ion_radius + 1):
                return new_coord
        print(f"Warning: Failed to place ion {ion} after {max_attempts} attempts")
        return None

    def place_ions_conmodel(ions, conmodel, volume, max_attempts=100):
        instructions = []
        ion_coords = []
        ion_radii = []
        for cc, ion in enumerate(ions):
            z_coords, conc_profile = conmodel[cc]
            z_positions = discretize_concentration(ion, conc_profile, z_coords, volume)
            for zpos in z_positions:
                new_coord = place_ion(ion, volume, ion_coords, ion_radii, zpos=zpos,
                                      max_attempts=max_attempts)
                if new_coord is not None:
                    ion_coords.append(new_coord)
                    ion_radii.append(ion.estimate_specie_radius())
                    instructions.append((ion.to_universe(), new_coord, "fixed"))
        return instructions

    def place_ions_random(ions, nions, volume, to_center, max_attempts=100):
        instructions = []
        ion_coords = []
        ion_radii = []
        for cc, ion in enumerate(ions):
            nrep = nions if isinstance(nions, int) else nions[cc]
            for _ in range(nrep):
                new_coord = volume / 2 if to_center else place_ion(ion, volume,
                                                                    ion_coords,
                                                                    ion_radii,
                                                                    max_attempts=max_attempts)
                if new_coord is not None:
                    ion_coords.append(new_coord)
                    ion_radii.append(ion.estimate_specie_radius())
                    instructions.append((ion.to_universe(), new_coord, "fixed"))
        return instructions

    max_attempts = 100
    volume = np.array(volume)

    if conmodel is not None:
        assert len(conmodel) == len(ions), "Need one profile per specie"
        return place_ions_conmodel(ions, conmodel, volume, max_attempts=max_attempts)

    if ion_pos == "box":
        return [(ion.to_universe(), nions if isinstance(nions, int) else nions[cc], "box") for cc, ion in enumerate(ions)]

    to_center = ion_pos == "center"
    if ion_pos == "left":
        volume[2] /= 2

    return place_ions_random(ions, nions, volume, to_center, max_attempts=max_attempts)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def make_solvent_box(
    species: List[Any],
    solvent: Optional[Any],
    ions: Optional[List[Any]],
    volume: List[float],
    density: Optional[float],
    nions: Optional[Union[int, List[int]]],
    concentration: Optional[float],
    conmodel: Optional[Dict[int, Tuple[List[float], List[float]]]],
    ion_pos: Optional[str],
    nsolvent=None,
    tolerance: float = 2.0,
    ratio: Optional[List[float]] = None,
) -> Optional[mda.Universe]:
    """
    Build a solvent box with optional ionic species.

    Parameters
    ----------
    species : list
        All registered species (used for topology look-up after PACKMOL).
    solvent : Specie, list of Specie, or None
        Solvent molecule(s). Pass a list for mixed-solvent boxes.
    ions : list or None
        Ionic species to dissolve.
    volume : list of float
        Box dimensions [x, y, z] in Å.
    density : float or None
        Solvent density in g/cm³.  Ignored when *nsolvent* is given.
    nions : int, list, or None
        Number of each ionic species.
    concentration : float or None
        Ionic concentration in Molar (alternative to *nions*).
    conmodel : dict or None
        Spatially varying concentration model.
    ion_pos : str or None
        Ion placement strategy: ``"box"``, ``"center"``, ``"left"``, or
        ``None`` for random fixed placement.
    nsolvent : int, list, or None
        Number of solvent molecules.  List length must match *solvent* list.
    tolerance : float
        PACKMOL minimum inter-molecule distance (Å).  Default 2.0.
    ratio : list of float or None
        Molar mixing ratio for multi-solvent boxes.  Must be used with
        either *density* or *nsolvent* (int).

    Returns
    -------
    MDAnalysis.Universe or None
    """

    # Normalise solvent to a list (backward-compatible).
    if solvent is None:
        solvents = []
    elif isinstance(solvent, (list, tuple)):
        solvents = list(solvent)
    else:
        solvents = [solvent]

    _validate_solvent_box_parameters(nions, concentration, conmodel, ions,
                                     solvents, density, nsolvent, ratio)

    # convert concentration to number of ions
    if concentration is not None:
        nions = int(concentration * np.prod(volume) * units.mol / ((units.m / 10) ** 3))

    instructions = []

    # ions
    if (conmodel is not None) or (nions is not None and ions is not None):
        ion_instr = populate_with_ions(ions, nions, volume, ion_pos=ion_pos,
                                       conmodel=conmodel)
        instructions.extend(ion_instr)

    # solvents
    if solvents:
        solvent_volume = 1e-24 * np.prod(volume)

        if ratio is not None:
            ratio_arr = np.array(ratio, dtype=float)
            if nsolvent is not None:
                counts = [max(1, int(r / ratio_arr.sum() * nsolvent)) for r in ratio_arr]
            else:
                masses = np.array([s.atoms.get_masses().sum() for s in solvents])
                mass_per_unit = float(np.dot(ratio_arr, masses))
                n_units = units.mol * density * solvent_volume / mass_per_unit
                counts = [max(1, int(r * n_units)) for r in ratio_arr]
        elif isinstance(nsolvent, (list, tuple)):
            counts = list(nsolvent)
        elif nsolvent is not None:
            counts = [int(nsolvent)]
        else:
            mass = solvents[0].atoms.get_masses().sum()
            counts = [int(units.mol * density * (1.0 / mass) * solvent_volume)]

        for sp, n in zip(solvents, counts):
            if n > 0:
                instructions.append([sp.to_universe(), n, "box"])

    universe = populate_box(volume, instructions, tolerance=tolerance)

    if universe is None:
        return None

    species_dict = {specie.residues.resnames[0]: specie for specie in species}

    alist = []
    for res in universe.residues:
        resname = res.resname
        if resname in species_dict:
            nmol = species_dict[resname].copy()
            nmol.atoms.positions = res.atoms.positions
            alist.append(nmol.atoms)

    solution = mda.Merge(*alist)
    solution.dimensions = volume + [90, 90, 90]

    return solution
