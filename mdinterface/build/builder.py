#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fluent builder API for assembling multi-layer simulation boxes.
"""

import logging
from typing import List, Optional, Union, Tuple, Any

from mdinterface.utils.auxiliary import find_smallest_missing, label_to_element
from mdinterface.io.lammpswriter import DATAWriter, write_lammps_coefficients
from mdinterface.build.box import make_interface_slab, make_solvent_box, add_component

import ase
import MDAnalysis as mda
import numpy as np
import shutil

#%%

logger = logging.getLogger("mdinterface.builder")


def _configure_logger(level) -> None:
    """
    Attach a StreamHandler to the mdinterface loggers and set their level.

    Parameters
    ----------
    level : bool, int, or str
        ``True`` → INFO, ``False`` → WARNING, a :mod:`logging` integer
        constant (e.g. ``logging.DEBUG``), or a string (``"DEBUG"``,
        ``"INFO"``, ``"WARNING"`` …).
    """
    if isinstance(level, bool):
        level = logging.INFO if level else logging.WARNING
    elif isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    fmt = logging.Formatter("[mdinterface] %(levelname)-8s %(message)s")

    for name in ("mdinterface.builder", "mdinterface.box"):
        lg = logging.getLogger(name)
        if not any(isinstance(h, logging.StreamHandler) for h in lg.handlers):
            handler = logging.StreamHandler()
            handler.setFormatter(fmt)
            lg.addHandler(handler)
        lg.setLevel(level)
        lg.propagate = False  # prevent double-printing via the root logger


class BoxBuilder:
    """
    Fluent builder for assembling layered simulation boxes.

    Layers are stacked in the order they are added. Slabs (solid interfaces),
    solvent regions, and vacuum gaps can be freely combined.

    Example
    -------
    >>> simbox = BoxBuilder(xysize=[15, 15])
    >>> simbox.add_slab(gold, nlayers=3)
    >>> simbox.add_solvent(water, ions=[na, cl], nions=[5, 5], zdim=25, density=1.0)
    >>> simbox.add_slab(gold, nlayers=3)
    >>> simbox.add_vacuum(zdim=5)
    >>> simbox.build()
    >>> simbox.write_lammps("data.lammps")
    """

    def __init__(
        self,
        xysize: Union[List[float], Tuple[float, float]],
        verbose: Union[None, bool, int, str] = None,
    ) -> None:
        if verbose is not None:
            _configure_logger(verbose)
        xsize, ysize = self._validate_xysize(xysize)
        self._xsize = xsize
        self._ysize = ysize
        self._layers: List[dict] = []
        self._all_species: List[Any] = []
        self._universe: Optional[mda.Universe] = None
        logger.debug("BoxBuilder initialised: xysize=[%.2f, %.2f] Å", xsize, ysize)

    # ------------------------------------------------------------------
    # Layer-adding methods (chainable)
    # ------------------------------------------------------------------

    def add_slab(self, species: Any, nlayers: int = 1) -> "BoxBuilder":
        """
        Add a solid-interface layer (slab).

        Parameters
        ----------
        species : Specie
            The unit-cell species to tile into a slab.
        nlayers : int
            Number of unit-cell layers to stack along Z.
        """
        slab_sp = species.copy()
        slab_idx = sum(1 for lay in self._layers if lay["type"] == "slab")
        suffix = f"_s{slab_idx}"
        for atom in slab_sp._stype:
            atom.set_label(atom.label + suffix)
        self._register(slab_sp)
        self._layers.append({"type": "slab", "species": slab_sp, "nlayers": nlayers})
        logger.debug("Layer added — slab: species=%s, nlayers=%d, suffix=%s",
                     getattr(slab_sp, "resname", "?"), nlayers, suffix)
        return self

    def add_solvent(
        self,
        solvent: Optional[Any] = None,
        ions: Optional[List[Any]] = None,
        nions: Optional[Union[int, List[int]]] = None,
        zdim: Optional[float] = None,
        density: Optional[float] = None,
        nsolvent: Optional[int] = None,
        concentration: Optional[float] = None,
        conmodel: Optional[dict] = None,
        ion_pos: Optional[str] = None,
        dilate: float = 1.0,
        packmol_tolerance: float = 2.0,
    ) -> "BoxBuilder":
        """
        Add a solvent (liquid) layer, optionally with dissolved ions.

        Parameters
        ----------
        solvent : Specie or None
            Solvent molecule species. Pass ``None`` for an ion-only region.
        ions : list of Specie, optional
            Ionic species to dissolve.
        nions : int or list of int, optional
            Number of each ionic species (alternative to *concentration*).
        zdim : float
            Thickness of this region in Angstroms.
        density : float, optional
            Solvent density in g/cm³ (alternative to *nsolvent*).
        nsolvent : int, optional
            Explicit number of solvent molecules (alternative to *density*).
        concentration : float, optional
            Ionic concentration in Molar (alternative to *nions*).
        conmodel : dict, optional
            Spatially varying concentration model.
        ion_pos : str, optional
            Ion placement strategy passed to PACKMOL
            (``"random"``, ``"center"`` …).
        dilate : float, optional
            Dilation factor > 1 for concentrated systems. PACKMOL will pack
            into a box ``dilate`` times taller at ``density / dilate``,
            keeping the number of solvent molecules identical. The larger box
            gives PACKMOL more breathing room; subsequent NpT equilibration
            compresses the system to the correct density. Default is ``1.0``
            (no dilation).
        packmol_tolerance : float, optional
            Minimum distance (Å) between atoms of different molecules during
            PACKMOL packing. Reduce from the default of ``2.0`` for very
            concentrated systems if PACKMOL cannot find a valid packing.
        """
        if zdim is None:
            raise ValueError("'zdim' is required for add_solvent()")
        if dilate <= 0:
            raise ValueError(f"'dilate' must be positive, got {dilate}")
        if packmol_tolerance <= 0:
            raise ValueError(f"'packmol_tolerance' must be positive, got {packmol_tolerance}")

        solv_copy = solvent.copy() if solvent is not None else None
        ions_copies = [ion.copy() for ion in (ions or [])]
        self._register(*([solv_copy] if solv_copy else []), *ions_copies)

        self._layers.append({
            "type":               "solvent",
            "solvent":            solv_copy,
            "ions":               ions_copies,
            "nions":              nions,
            "zdim":               zdim,
            "density":            density,
            "nsolvent":           nsolvent,
            "concentration":      concentration,
            "conmodel":           conmodel,
            "ion_pos":            ion_pos,
            "dilate":             dilate,
            "packmol_tolerance":  packmol_tolerance,
        })
        logger.debug(
            "Layer added — solvent: solvent=%s, zdim=%.1f Å, density=%s, "
            "nions=%s, dilate=%.2f, tolerance=%.1f",
            getattr(solv_copy, "resname", "None"), zdim, density, nions,
            dilate, packmol_tolerance,
        )
        return self

    def add_vacuum(self, zdim: float = 0.0) -> "BoxBuilder":
        """
        Add an empty vacuum gap.

        Parameters
        ----------
        zdim : float
            Thickness of the vacuum gap in Angstroms.
        """
        self._layers.append({"type": "vacuum", "zdim": zdim})
        logger.debug("Layer added — vacuum: zdim=%.1f Å", zdim)
        return self

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    def build(
        self,
        padding: float = 1.5,
        center: bool = False,
        layered: bool = False,
        match_cell: bool = False,
        hijack: Optional[ase.Atoms] = None,
    ) -> "BoxBuilder":
        """
        Assemble all layers into a simulation box.

        Parameters
        ----------
        padding : float
            Spacing (Å) inserted between adjacent layers.
        center : bool
            If True, shift the whole system by half-box along Z after assembly.
        layered : bool
            Assign distinct molecule indices to each slab layer for LAMMPS.
        match_cell : bool
            Deform slabs to exactly match XY cell dimensions.
        hijack : ase.Atoms, optional
            After assembly, override both positions and cell dimensions with
            this external ``ase.Atoms`` object. Useful when you have a
            pre-relaxed structure you want to map the topology onto.

        Returns
        -------
        BoxBuilder
            Returns *self* so you can chain output methods.
        """
        n_layers = len(self._layers)
        logger.info("Starting build: %d layers, padding=%.2f Å", n_layers, padding)

        self._update_topology_indexes()

        xsize, ysize = self._xsize, self._ysize

        # ------------------------------------------------------------------
        # First pass: build slab supercells and determine actual cell XY.
        # ------------------------------------------------------------------
        slab_dims = []   # (xi, yi) for each slab, in order
        for layer in self._layers:
            if layer["type"] == "slab":
                tslab = make_interface_slab(
                    layer["species"], xsize, ysize, layers=layer["nlayers"]
                )
                layer["_slab"] = tslab
                if tslab is not None:
                    xi = np.dot([1, 0, 0], tslab.atoms.cell @ [1, 0, 0])
                    yi = np.dot([0, 1, 0], tslab.atoms.cell @ [0, 1, 0])
                    zi = np.dot([0, 0, 1], tslab.atoms.cell @ [0, 0, 1])
                    slab_dims.append((xi, yi))
                    xsize = max(xsize, xi)
                    ysize = max(ysize, yi)
                    logger.debug(
                        "  Slab fitted: %s -> %.3f x %.3f x %.3f A, %d atoms",
                        getattr(layer["species"], "resname", "?"),
                        xi, yi, zi, len(tslab.atoms),
                    )

        if slab_dims:
            logger.info("Cell size after slab fitting: x=%.3f Å, y=%.3f Å", xsize, ysize)
            self._check_xy_mismatch(slab_dims, xsize, ysize)

        # ------------------------------------------------------------------
        # All species as universes for PACKMOL residue-name lookup.
        # ------------------------------------------------------------------
        all_sp_univs = [sp.to_universe() for sp in self._all_species]

        # ------------------------------------------------------------------
        # Second pass: stack layers.
        # ------------------------------------------------------------------
        system = None
        zdim = 0.0

        for ii, layer in enumerate(self._layers):
            ltype = layer["type"]
            logger.info("Stacking layer %d/%d: %s", ii + 1, n_layers, ltype)

            if ltype == "slab":
                slab_u = layer["_slab"].to_universe(
                    layered=layered, match_cell=match_cell, xydim=[xsize, ysize]
                )
                system, zdim = add_component(system, slab_u, zdim, padding=padding)
                logger.info("  -> %d atoms, zdim -> %.3f A",
                            len(slab_u.atoms), zdim)

            elif ltype == "solvent":
                dilate   = layer["dilate"]
                eff_zdim = layer["zdim"] * dilate
                eff_rho  = layer["density"] / dilate if layer["density"] is not None else None

                if dilate != 1.0:
                    logger.info(
                        "  Dilation x%.2f: packing into %.1f A at %.3f g/cm3"
                        " (target zdim=%.1f A)",
                        dilate, eff_zdim,
                        eff_rho if eff_rho is not None else float("nan"),
                        layer["zdim"],
                    )

                ions_list = layer["ions"] if layer["ions"] else None
                solv_box = make_solvent_box(
                    all_sp_univs,
                    layer["solvent"].to_universe() if layer["solvent"] else None,
                    ions_list,
                    [xsize, ysize, eff_zdim],
                    eff_rho,
                    layer["nions"],
                    layer["concentration"],
                    layer["conmodel"],
                    layer["ion_pos"],
                    layer["nsolvent"],
                    layer["packmol_tolerance"],
                )
                if solv_box is not None:
                    logger.info("  -> %d atoms, zdim -> %.3f A",
                                len(solv_box.atoms), zdim + eff_zdim)
                else:
                    logger.warning("  Solvent box is empty - check packmol.log")
                system, zdim = add_component(system, solv_box, zdim, padding=padding)

            elif ltype == "vacuum":
                zdim += layer["zdim"]
                logger.info("  -> zdim -> %.3f A", zdim)

        system.dimensions = [xsize, ysize, zdim] + [90, 90, 90]
        logger.info("Assembly complete: %d atoms, zdim=%.3f A",
                    len(system.atoms), zdim)

        if center:
            system.atoms.translate([0, 0, zdim / 2])
            _ = system.atoms.wrap()
            logger.info("  -> system shifted by half-box along Z")

        if hijack is not None:
            system.dimensions = hijack.get_cell_lengths_and_angles()
            system.atoms.positions = hijack.get_positions()
            logger.info("  -> positions and cell overridden by hijack ase.Atoms")

        self._universe = system
        self._xsize = xsize
        self._ysize = ysize
        return self

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def write_lammps(
        self,
        filename: str = "data.lammps",
        atom_style: str = "full",
        write_coeff: bool = True,
    ) -> "BoxBuilder":
        """
        Write a LAMMPS data file (and optional force-field coefficients).

        Parameters
        ----------
        filename : str
            Output file path.
        atom_style : str
            LAMMPS atom style (``"full"`` or ``"atomic"``).
        write_coeff : bool
            Whether to write force-field coefficient blocks.

        Returns
        -------
        BoxBuilder
            Returns *self* for further chaining.
        """
        if self._universe is None:
            raise RuntimeError("Call build() before write_lammps().")

        logger.info("Writing LAMMPS data file: %s (atom_style=%s, write_coeff=%s)",
                    filename, atom_style, write_coeff)
        system = self._universe.copy()

        if not write_coeff:
            for attr in ["bonds", "angles", "dihedrals", "impropers"]:
                try:
                    system.del_TopologyAttr(attr)
                except Exception:
                    pass

        with DATAWriter(filename) as dt:
            dt.write(system.atoms, atom_style=atom_style)

        if write_coeff:
            temp_file = "tmp_data.lammps"
            sorted_attrs = {
                "atoms":     self.get_sorted_attribute("atoms"),
                "bonds":     self.get_sorted_attribute("bonds"),
                "angles":    self.get_sorted_attribute("angles"),
                "dihedrals": self.get_sorted_attribute("dihedrals"),
                "impropers": self.get_sorted_attribute("impropers"),
            }
            with open(filename, "r") as ffile, open(temp_file, "w") as tfile:
                for fl in ffile:
                    if fl.startswith("Atoms"):
                        write_lammps_coefficients(system, sorted_attrs, fout=tfile)
                    tfile.write(fl)
            shutil.move(temp_file, filename)

        try:
            nbonds = len(system.atoms.bonds)
        except Exception:
            nbonds = 0
        logger.info("Written: %d atoms, %d bonds -> %s", len(system.atoms), nbonds, filename)
        return self

    def to_ase(self) -> ase.Atoms:
        """
        Convert the assembled system to an ``ase.Atoms`` object.

        Returns
        -------
        ase.Atoms
            The simulation box as an ASE Atoms object with cell and PBC set.

        Raises
        ------
        RuntimeError
            If called before :meth:`build`.
        """
        if self._universe is None:
            raise RuntimeError("Call build() before to_ase().")

        system = self._universe
        positions = system.atoms.positions
        masses    = system.atoms.masses
        labels    = system.atoms.types

        symbols = [label_to_element(lab, mas) for lab, mas in zip(labels, masses)]
        atoms = ase.Atoms(symbols=symbols, positions=positions)

        if system.dimensions is not None:
            atoms.set_cell(system.dimensions)
            atoms.set_pbc(True)

        try:
            if system.atoms.charges is not None:
                atoms.set_initial_charges(system.atoms.charges)
        except Exception:
            pass

        return atoms

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def universe(self) -> Optional[mda.Universe]:
        """The assembled MDAnalysis Universe (available after build())."""
        return self._universe

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _validate_xysize(
        xysize: Union[List[float], Tuple[float, float]]
    ) -> Tuple[float, float]:
        if not isinstance(xysize, (list, tuple, np.ndarray)):
            raise TypeError(
                f"xysize must be a list, tuple, or numpy array, "
                f"got {type(xysize).__name__}."
            )
        if len(xysize) != 2:
            raise ValueError(
                f"xysize must have exactly 2 elements [x, y], "
                f"got {len(xysize)}: {xysize}."
            )
        try:
            xsize, ysize = float(xysize[0]), float(xysize[1])
        except (ValueError, TypeError) as e:
            raise ValueError(f"xysize elements must be numeric: {e}") from e
        if xsize <= 0 or ysize <= 0:
            raise ValueError(
                f"xysize values must be positive, got [{xsize}, {ysize}]."
            )
        return xsize, ysize

    def _check_xy_mismatch(
        self, slab_dims: list, xsize: float, ysize: float
    ) -> None:
        """Warn if slabs have inconsistent XY sizes or differ much from requested."""
        threshold = 0.1   # 10 % relative difference triggers a warning

        # Check consistency across slabs.
        if len(slab_dims) > 1:
            xs = [d[0] for d in slab_dims]
            ys = [d[1] for d in slab_dims]
            x_spread = (max(xs) - min(xs)) / xsize
            y_spread = (max(ys) - min(ys)) / ysize
            if x_spread > threshold or y_spread > threshold:
                logger.warning(
                    "Slab XY dimensions are inconsistent: "
                    "x range [%.2f – %.2f Å], y range [%.2f – %.2f Å]. "
                    "Consider using match_cell=True in build().",
                    min(xs), max(xs), min(ys), max(ys),
                )

        # Check against original requested xysize.
        req_x, req_y = self._xsize, self._ysize
        dx = abs(xsize - req_x) / req_x
        dy = abs(ysize - req_y) / req_y
        if dx > threshold or dy > threshold:
            logger.warning(
                "Actual cell size (%.2f × %.2f Å) differs from requested "
                "(%.2f × %.2f Å) by %.1f%% / %.1f%%. "
                "The slab periodicity determines the final XY size.",
                xsize, ysize, req_x, req_y, dx * 100, dy * 100,
            )

    def _register(self, *species: Any) -> None:
        """Add species to the global registry if not already present."""
        for sp in species:
            if sp is not None and not any(s is sp for s in self._all_species):
                self._all_species.append(sp)

    def _update_topology_indexes(self) -> None:
        logger.debug("Updating topology indexes for %d species", len(self._all_species))
        nitems: dict = {
            "_btype": [],
            "_atype": [],
            "_dtype": [],
            "_itype": [],
        }
        for attribute in nitems:
            for specie in self._all_species:
                for attr in specie.__getattribute__(attribute):
                    if attr.id not in nitems[attribute]:
                        nitems[attribute].append(attr.id)
                    else:
                        idx = find_smallest_missing(nitems[attribute], start=1)
                        attr.set_id(idx)
                        nitems[attribute].append(attr.id)

        # Sort atom types alphabetically by extended label.
        atom_types = []
        for specie in self._all_species:
            atom_types.extend([stype.extended_label for stype in specie._stype])
        atom_types.sort()

        for specie in self._all_species:
            for stype in specie._stype:
                idx = np.argwhere(stype.extended_label == np.array(atom_types))[0][0]
                stype.set_id(idx + 1)

        logger.debug(
            "Topology: %d atom type(s), %d bond type(s), "
            "%d angle type(s), %d dihedral type(s), %d improper type(s)",
            len(atom_types), len(nitems["_btype"]),
            len(nitems["_atype"]), len(nitems["_dtype"]), len(nitems["_itype"]),
        )

    def get_sorted_attribute(self, attribute: str) -> list:
        attr_map = {
            "bonds":     "_btype",
            "angles":    "_atype",
            "dihedrals": "_dtype",
            "impropers": "_itype",
            "atoms":     "_stype",
        }
        key = attr_map.get(attribute.lower(), attribute)

        indexes = []
        attributes = []
        for specie in self._all_species:
            for attr in specie.__getattribute__(key):
                indexes.append(attr.id)
                attributes.append(attr)

        return [attributes[ii] for ii in np.argsort(indexes)]
