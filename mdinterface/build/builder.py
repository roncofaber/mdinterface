#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fluent builder API for assembling multi-layer simulation boxes.
"""

from typing import List, Optional, Union, Tuple, Any

from mdinterface.utils.auxiliary import find_smallest_missing, label_to_element
from mdinterface.io.lammpswriter import DATAWriter, write_lammps_coefficients
from mdinterface.build.box import make_interface_slab, make_solvent_box, add_component

import ase
import MDAnalysis as mda
import numpy as np
import shutil

#%%

class BoxBuilder:
    """
    Fluent builder for assembling layered simulation boxes.

    Layers are stacked in the order they are added. Slabs (solid interfaces),
    solvent regions, and vacuum gaps can be freely combined.

    Example
    -------
    >>> system = (
    ...     BoxBuilder(xysize=[15, 15])
    ...         .add_slab(gold, nlayers=3)
    ...         .add_solvent(water, ions=[na, cl], nions=[5, 5], zdim=25, density=1.0)
    ...         .add_slab(gold, nlayers=3)
    ...         .add_vacuum(zdim=5)
    ...         .build()
    ...         .write_lammps("data.lammps")
    ... )
    """

    def __init__(self, xysize: Union[List[float], Tuple[float, float]]) -> None:
        if len(xysize) != 2:
            raise ValueError(f"xysize must have exactly 2 elements, got {len(xysize)}")
        self._xsize = float(xysize[0])
        self._ysize = float(xysize[1])
        self._layers: List[dict] = []
        self._all_species: List[Any] = []   # ordered list of Specie objects
        self._universe: Optional[mda.Universe] = None

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
            Number of layers to stack along Z.
        """
        slab_sp = species.copy()
        slab_idx = sum(1 for lay in self._layers if lay["type"] == "slab")
        suffix = f"_s{slab_idx}"
        for atom in slab_sp._stype:
            atom.set_label(atom.label + suffix)
        self._register(slab_sp)
        self._layers.append({"type": "slab", "species": slab_sp, "nlayers": nlayers})
        return self

    def add_solvent(
        self,
        solvent: Any,
        ions: Optional[List[Any]] = None,
        nions: Optional[Union[int, List[int]]] = None,
        zdim: Optional[float] = None,
        density: Optional[float] = None,
        nsolvent: Optional[int] = None,
        concentration: Optional[float] = None,
        conmodel: Optional[dict] = None,
        ion_pos: Optional[str] = None,
    ) -> "BoxBuilder":
        """
        Add a solvent (liquid) layer, optionally with dissolved ions.

        Parameters
        ----------
        solvent : Specie
            Solvent molecule species.
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
            Ion placement strategy passed to PACKMOL (``"random"``, ``"center"`` …).
        """
        if zdim is None:
            raise ValueError("'zdim' is required for add_solvent()")
        solv_copy = solvent.copy()
        ions_copies = [ion.copy() for ion in (ions or [])]
        self._register(solv_copy, *ions_copies)
        self._layers.append({
            "type":          "solvent",
            "solvent":       solv_copy,
            "ions":          ions_copies,
            "nions":         nions,
            "zdim":          zdim,
            "density":       density,
            "nsolvent":      nsolvent,
            "concentration": concentration,
            "conmodel":      conmodel,
            "ion_pos":       ion_pos,
        })
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

        Returns
        -------
        BoxBuilder
            Returns *self* so you can chain ``.write_lammps()``.
        """
        self._update_topology_indexes()

        xsize, ysize = self._xsize, self._ysize

        # First pass: build slab supercells and determine actual cell size.
        for layer in self._layers:
            if layer["type"] == "slab":
                tslab = make_interface_slab(
                    layer["species"], xsize, ysize, layers=layer["nlayers"]
                )
                layer["_slab"] = tslab
                if tslab is not None:
                    xi = np.dot([1, 0, 0], tslab.atoms.cell @ [1, 0, 0])
                    yi = np.dot([0, 1, 0], tslab.atoms.cell @ [0, 1, 0])
                    xsize = max(xsize, xi)
                    ysize = max(ysize, yi)

        # All species as universes for PACKMOL residue-name lookup.
        all_sp_univs = [sp.to_universe() for sp in self._all_species]

        # Second pass: stack layers.
        system = None
        zdim = 0.0

        for layer in self._layers:
            ltype = layer["type"]

            if ltype == "slab":
                slab_u = layer["_slab"].to_universe(
                    layered=layered, match_cell=match_cell, xydim=[xsize, ysize]
                )
                system, zdim = add_component(system, slab_u, zdim, padding=padding)

            elif ltype == "solvent":
                ions_list = layer["ions"] if layer["ions"] else None
                solv_box = make_solvent_box(
                    all_sp_univs,
                    layer["solvent"].to_universe(),
                    ions_list,
                    [xsize, ysize, layer["zdim"]],
                    layer["density"],
                    layer["nions"],
                    layer["concentration"],
                    layer["conmodel"],
                    layer["ion_pos"],
                    layer["nsolvent"],
                )
                system, zdim = add_component(system, solv_box, zdim, padding=padding)

            elif ltype == "vacuum":
                zdim += layer["zdim"]

        system.dimensions = [xsize, ysize, zdim] + [90, 90, 90]

        if center:
            system.atoms.translate([0, 0, zdim / 2])
            _ = system.atoms.wrap()

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
    # Internal helpers (adapted from SimulationBox)
    # ------------------------------------------------------------------

    def _register(self, *species: Any) -> None:
        """Add species to the global registry if not already present."""
        for sp in species:
            if sp is not None and not any(s is sp for s in self._all_species):
                self._all_species.append(sp)

    def _update_topology_indexes(self) -> None:
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
