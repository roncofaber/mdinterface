#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SimCell: layer-by-layer builder for assembling multi-layer MD simulation boxes.

Add slabs, solvent regions, and vacuum gaps one step at a time, then call
:meth:`SimCell.build` to pack and assemble the full system.
"""

import logging
from collections import Counter
from typing import List, Optional, Union, Tuple, Any

from mdinterface.utils.auxiliary import find_smallest_missing, label_to_element
from mdinterface.io.lammpswriter import DATAWriter, write_lammps_coefficients
from mdinterface.build.box import make_interface_slab, add_component
from mdinterface.build.solvent import make_solvent_box

import ase
import MDAnalysis as mda
import numpy as np
import shutil

from mdinterface.utils.logger import set_verbosity, log_banner, log_header, log_subheader

def _get_mdi_version() -> str:
    """Return the mdinterface version, read at call time to avoid stale metadata."""
    import sys
    mdi = sys.modules.get("mdinterface")
    return getattr(mdi, "__version__", "unknown")

#%%

logger = logging.getLogger(__name__)


def _log_layer_result(n_atoms, dims, zdim, layer_zdim, extra_lines=()):
    """
    Emit the standard size / z lines after a layer is assembled.

    Parameters
    ----------
    n_atoms : int or None
        Atom count; omit the size line if None.
    dims : tuple (x, y, z) or None
        Cell dimensions in Å; omit the size line if None.
    zdim : float
        Running total Z height after this layer.
    layer_zdim : float
        Z contributed by this layer.
    extra_lines : iterable of str
        Additional lines (without leading ``>`` or indent) emitted between the
        size line and the total-z line.  Each is prefixed with ``    > ``.
    """
    if n_atoms is not None and dims is not None:
        x, y, z = dims
        logger.info("  ├> %d atoms  |  %.3f x %.3f x %.3f Å", n_atoms, x, y, z)
    for line in extra_lines:
        logger.info("  ├> %s", line)
    logger.info("  └─> total z: %.2f Å  (+%.2f Å)", zdim, layer_zdim)


class SimCell:
    """
    Fluent builder for assembling layered simulation boxes.

    Layers are stacked in the order they are added. Slabs (solid interfaces),
    solvent regions, and vacuum gaps can be freely combined.

    Example
    -------
    >>> simbox = SimCell(xysize=[15, 15])
    >>> simbox.add_slab(gold, nlayers=3)
    >>> simbox.add_solvent(water, solute=[na, cl], nsolute=[5, 5], zdim=25, density=1.0)
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
        """
        Parameters
        ----------
        xysize : list or tuple of float
            XY dimensions of the simulation box in Å, e.g. ``[15.0, 15.0]``.
        verbose : None, bool, int, or str, optional
            Controls package-wide log verbosity via :func:`set_verbosity`.
            ``None`` (default) leaves the current level unchanged.

            Integer scale:

            - ``0`` / ``False`` — WARNING (quiet)
            - ``1`` / ``True``  — INFO  (normal)
            - ``2``             — DEBUG (detailed)
            - ``3``, ``4``, …  — DEBUG (same as 2 for values < 10)

            Integers ≥ 10 are treated as raw ``logging`` level constants.
            Strings are resolved by name (``"DEBUG"``, ``"INFO"``, …).
        """
        if verbose is not None:
            set_verbosity(verbose)
        
        xsize, ysize = self._validate_xysize(xysize)
        self._xsize = xsize
        self._ysize = ysize
        self._layers: List[dict] = []
        self._all_species: List[Any] = []
        self._universe: Optional[mda.Universe] = None
        
        log_banner(logger, "mdinterface :: SimCell", f"version {_get_mdi_version()}")
        logger.info("  └─> xysize: %.2f x %.2f Å", xsize, ysize)

    # ------------------------------------------------------------------
    # Layer-adding methods (chainable)
    # ------------------------------------------------------------------

    def add_slab(self, species: Any, nlayers: int = 1) -> None:
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
        self._layers.append({"type": "slab", "label": "slab",
                              "species": slab_sp, "nlayers": nlayers})
        logger.info("  + slab     %s,  %d layer(s)",
                    getattr(slab_sp, "resname", "?"), nlayers)

    def add_prebuilt(self, species: Any, nlayers: int = 1) -> None:
        """
        Add a pre-built layer whose atomic positions come from a prior MD run.

        Semantically equivalent to :meth:`add_slab` but intended for species
        whose positions have already been set (e.g. via
        :meth:`~mdinterface.core.specie.Specie.update_positions`).  No XY
        tiling is performed; the species cell is used as-is.

        Any pre-processing of the positions (centering, trimming, …) should
        be done on the species before passing it here.

        Parameters
        ----------
        species : Specie or Polymer
            Species with positions already set from a prior trajectory.
        nlayers : int
            Number of repeat units to stack along Z (usually 1).
        """
        slab_sp = species.copy()
        slab_idx = sum(1 for lay in self._layers if lay["type"] == "slab")
        suffix = f"_s{slab_idx}"
        for atom in slab_sp._stype:
            atom.set_label(atom.label + suffix)
        self._register(slab_sp)
        self._layers.append({"type": "slab", "label": "prebuilt",
                              "species": slab_sp, "nlayers": nlayers})
        logger.info("  + prebuilt %s,  %d layer(s)",
                    getattr(slab_sp, "resname", "?"), nlayers)

    def add_solvent(
        self,
        solvent: Optional[Any] = None,
        solute: Optional[List[Any]] = None,
        nsolute: Optional[Union[int, List[int]]] = None,
        zdim: Optional[float] = None,
        density: Optional[float] = None,
        nsolvent: Optional[Union[int, List[int]]] = None,
        concentration: Optional[float] = None,
        conmodel: Optional[dict] = None,
        solute_pos: Optional[str] = None,
        dilate: float = 1.0,
        packmol_tolerance: float = 2.0,
        ratio: Optional[List[float]] = None,
    ) -> None:
        """
        Add a solvent (liquid) layer, optionally with dissolved species.

        Parameters
        ----------
        solvent : Specie or list of Specie or None
            Solvent molecule(s). Pass ``None`` for a solute-only region.
        solute : list of Specie, optional
            Species to dissolve (ions, neutral molecules, …).
        nsolute : int or list of int, optional
            Number of each solute species (alternative to *concentration*).
        zdim : float
            Thickness of this region in Angstroms.
        density : float, optional
            Solvent density in g/cm³ (alternative to *nsolvent*).
        nsolvent : int or list of int, optional
            Explicit number of solvent molecules (alternative to *density*).
        concentration : float, optional
            Solute concentration in Molar (alternative to *nsolute*).
        conmodel : dict, optional
            Spatially varying concentration model.
        solute_pos : str, optional
            Solute placement strategy:

            - ``None`` / ``"packmol"`` — PACKMOL random placement in the full
              box (default).
            - ``"left"``   — PACKMOL random placement in the left half
              (z ≤ zdim/2).
            - ``"right"``  — PACKMOL random placement in the right half
              (z ≥ zdim/2).
            - ``"center"`` — each molecule fixed at the box centre.
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

        # Normalise solvent to a list; copy each species.
        if solvent is None:
            solv_copies = []
        elif isinstance(solvent, (list, tuple)):
            solv_copies = [s.copy() for s in solvent]
        else:
            solv_copies = [solvent.copy()]

        solute_copies = [sp.copy() for sp in (solute or [])]
        self._register(*solv_copies, *solute_copies)

        self._layers.append({
            "type":               "solvent",
            "solvent":            solv_copies,
            "solute":             solute_copies,
            "nsolute":            nsolute,
            "zdim":               zdim,
            "density":            density,
            "nsolvent":           nsolvent,
            "concentration":      concentration,
            "conmodel":           conmodel,
            "solute_pos":         solute_pos,
            "dilate":             dilate,
            "packmol_tolerance":  packmol_tolerance,
            "ratio":              ratio,
        })
        solv_str = "+".join(getattr(s, "resname", "?") for s in solv_copies) or "ions"
        rho_str  = (f"ρ={density:.2f} g/cm³" if density is not None
                    else f"nsolvent={nsolvent}" if nsolvent is not None else "density=?")
        solute_info = ""
        if solute_copies:
            sol_names = "+".join(getattr(s, "resname", "?") for s in solute_copies)
            solute_info = f",  solute: {sol_names} (n={nsolute})"
        logger.info("  + solvent  %s,  zdim=%.1f Å,  %s%s",
                    solv_str, zdim, rho_str, solute_info)

    def add_vacuum(self, zdim: float = 0.0) -> None:
        """
        Add an empty vacuum gap.

        Parameters
        ----------
        zdim : float
            Thickness of the vacuum gap in Angstroms.
        """
        self._layers.append({"type": "vacuum", "zdim": zdim})
        logger.info("  + vacuum   zdim=%.1f Å", zdim)

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    def build(
        self,
        padding: float = 0.5,
        center: bool = False,
        layered: bool = False,
        match_cell: Union[bool, Any] = True,
        hijack: Optional[ase.Atoms] = None,
        stack_axis: str = "z",
    ) -> None:
        """
        Assemble all layers into a simulation box.

        Parameters
        ----------
        padding : float
            Spacing (Å) inserted between adjacent layers.
        center : bool
            If True, shift the system so the center of the first layer falls
            on the periodic boundary (z=0).  This is the standard convention
            for electrode/electrolyte slabs where one electrode straddles the
            cell edge.
        layered : bool
            Assign distinct molecule indices to each slab layer for LAMMPS.
        match_cell : bool or Specie
            Controls XY cell matching:

            - ``True`` — scale all slabs to the largest fitted XY cell
              (default). Ensures the solid/liquid interface is well-defined
              with no XY mismatch between layers.
            - ``False`` — no scaling; each slab keeps its natural tiled XY.
            - *Specie* — lock XY to that species' cell and stretch all other
              slabs to match. The reference species is left unscaled. Useful
              when a polymer or pre-relaxed slab should define the cell and
              electrodes should conform to it.
        hijack : ase.Atoms, optional
            After assembly, override both positions and cell dimensions with
            this external ``ase.Atoms`` object. Useful when you have a
            pre-relaxed structure you want to map the topology onto.
        stack_axis : str, optional
            Axis along which layers are stacked in the output file.
            Assembly always happens along Z; a final coordinate permutation
            reorients the box. Accepted values: ``"x"``, ``"y"``, ``"z"``
            (default ``"z"``).

        """
        stack_axis = stack_axis.lower()
        if stack_axis not in ("x", "y", "z"):
            raise ValueError(
                f"stack_axis must be 'x', 'y', or 'z', got {stack_axis!r}"
            )

        log_header(logger, "Build")
        logger.info("  ├> %d layers, padding=%.2f Å",
                    len(self._layers), padding)

        self._update_topology_indexes()

        xsize, ysize, cell_ref, do_match = self._resolve_xy(match_cell)
        xsize, ysize = self._fit_slabs(xsize, ysize, cell_ref)

        all_sp_univs = [sp.to_universe() for sp in self._all_species]
        system, zdim, first_layer_zdim = self._stack_layers(
            xsize, ysize, all_sp_univs, padding, layered, do_match
        )

        system.dimensions = [xsize, ysize, zdim] + [90, 90, 90]
        res_counts = Counter(res.resname for res in system.residues)

        log_header(logger, "Done")
        logger.info("  ├> %d atoms  |  %.3f x %.3f x %.3f Å",
                    len(system.atoms), xsize, ysize, zdim)
        logger.info("  ├> %d layers  |  %d residues",
                    len(self._layers), len(system.residues))
        parts = [f"{name} ({n} mol)" for name, n in res_counts.items()]
        for i in range(0, len(parts), 3):
            logger.info("  ├> %s", ",  ".join(parts[i:i + 3]))

        if center:
            shift = zdim - first_layer_zdim / 2
            system.atoms.translate([0, 0, shift])
            _ = system.atoms.wrap()
            logger.info("  └─> system centered on first layer (shift %.2f Å)", shift)

        if hijack is not None:
            system.dimensions = hijack.get_cell_lengths_and_angles()
            system.atoms.positions = hijack.get_positions()
            logger.info("  └─> positions and cell overridden by hijack ase.Atoms")

        if stack_axis != "z":
            system, xsize, ysize = self._apply_stack_axis(system, xsize, ysize, zdim, stack_axis)

        self._universe = system
        self._xsize = xsize
        self._ysize = ysize
        
        return

    # Output
    def write_lammps(
        self,
        filename: str = "data.lammps",
        atom_style: str = "full",
        write_coeff: bool = True,
    ) -> None:
        """
        Write a LAMMPS data file (and optional force-field coefficients).

        .. note::
            This method is only needed for classical MD with LAMMPS.  If you
            are using the assembled structure for AIMD, ML-MD, or any other
            workflow, use :meth:`to_ase` or :attr:`universe` instead — no
            force-field parameters are required for those paths.

        Parameters
        ----------
        filename : str
            Output file path.
        atom_style : str
            LAMMPS atom style (``"full"`` or ``"atomic"``).
        write_coeff : bool
            Whether to write force-field coefficient blocks.

        """
        if self._universe is None:
            raise RuntimeError("Call build() before write_lammps().")

        log_header(logger, "Output")
        logger.info("  ├> %s  (style=%s,  coeff=%s)", filename, atom_style, write_coeff)
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
        logger.info("  └─> %d atoms,  %d bonds  written", len(system.atoms), nbonds)
        return

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
    
    def to_universe(self):
        return self.universe

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

    def _resolve_xy(self, match_cell):
        """
        Determine the starting XY cell dimensions and match-cell settings.

        Returns
        -------
        xsize, ysize : float
            Starting XY dimensions in Å.
        cell_ref : Specie or None
            Reference species that locks the XY cell, or ``None``.
        do_match : bool
            Whether slab conversion should apply XY scaling.
        """
        cell_ref, do_match = self._resolve_match_cell(match_cell)

        if cell_ref is not None:
            # Pre-fit through make_interface_slab to get *tiled* dimensions.
            # For a polymer (xrep=yrep=1) this is a no-op; for a crystalline
            # slab it tiles to xysize and returns the correct fitted cell.
            ref_fitted = make_interface_slab(cell_ref, self._xsize, self._ysize)
            if ref_fitted is not None:
                xsize = np.dot([1, 0, 0], ref_fitted.atoms.cell @ [1, 0, 0])
                ysize = np.dot([0, 1, 0], ref_fitted.atoms.cell @ [0, 1, 0])
            else:
                xsize = cell_ref.atoms.get_cell()[0][0]
                ysize = cell_ref.atoms.get_cell()[1][1]
            logger.info("  ├> Reference cell: %s:  %.3f x %.3f Å",
                        getattr(cell_ref, "resname", "?"), xsize, ysize)
            logger.info("  └─> All slabs will be stretched to match it")
        else:
            xsize, ysize = self._xsize, self._ysize
            if do_match:
                logger.info("  └─> No reference cell: using largest XY slab")

        return xsize, ysize, cell_ref, do_match

    def _fit_slabs(self, xsize: float, ysize: float, cell_ref) -> Tuple[float, float]:
        """
        First pass: tile every slab layer and finalise the XY cell size.

        Tiled slabs are cached in ``layer["_slab"]``.  When no *cell_ref* is
        set, *xsize*/*ysize* grow to accommodate the largest slab footprint.

        Returns
        -------
        xsize, ysize : float
            Final XY dimensions in Å.
        """
        slab_layers = [l for l in self._layers if l["type"] == "slab"]
        if slab_layers:
            log_subheader(logger, "Slab tiling")

        slab_dims = []
        for layer in self._layers:
            if layer["type"] != "slab":
                continue
            tslab = make_interface_slab(
                layer["species"], xsize, ysize, layers=layer["nlayers"]
            )
            layer["_slab"] = tslab
            if tslab is None:
                continue
            xi = np.dot([1, 0, 0], tslab.atoms.cell @ [1, 0, 0])
            yi = np.dot([0, 1, 0], tslab.atoms.cell @ [0, 1, 0])
            zi = np.dot([0, 0, 1], tslab.atoms.cell @ [0, 0, 1])
            slab_dims.append((xi, yi))
            layer["_native_xy"] = (xi, yi)
            logger.info("  ├> %s:  %.3f x %.3f x %.3f Å,  %d atoms",
                        getattr(layer["species"], "resname", "?"),
                        xi, yi, zi, len(tslab.atoms))

        if slab_dims:
            if cell_ref is None:
                # Cell XY is defined by the largest fitted slab, not the
                # requested xysize (which was only used as a tiling target).
                xsize = max(d[0] for d in slab_dims)
                ysize = max(d[1] for d in slab_dims)
            logger.info("  └─> Final XY size: %.3f x %.3f Å", xsize, ysize)
            self._check_xy_mismatch(slab_dims, xsize, ysize)

        return xsize, ysize

    def _stack_layers(
        self,
        xsize: float,
        ysize: float,
        all_sp_univs: list,
        padding: float,
        layered: bool,
        do_match: bool,
    ) -> Tuple[mda.Universe, float]:
        """
        Second pass: assemble all layers into a single Universe.

        Returns
        -------
        system : mda.Universe
        zdim : float
            Total Z height in Å (excluding final cell padding).
        first_layer_zdim : float
            Z height after the first layer only, used for centering.
        """
        system = None
        zdim = 0.0
        first_layer_zdim = 0.0
        n_layers = len(self._layers)

        for ii, layer in enumerate(self._layers):
            ltype = layer["type"]
            tag = f"[{ii + 1}/{n_layers}]"

            if ltype == "slab":
                name  = getattr(layer["species"], "resname", "?")
                label = layer.get("label", "slab")
                log_subheader(logger, f"Layer {tag}")
                logger.info("  ├> %s: %s,  %d layer(s)", label, name, layer["nlayers"])
                slab_u = layer["_slab"].to_universe(
                    layered=layered, match_cell=do_match, xydim=[xsize, ysize]
                )
                zdim_before = zdim
                system, zdim = add_component(system, slab_u, zdim, padding=padding)
                sx, sy, sz = slab_u.dimensions[:3]
                native_xi, native_yi = layer.get("_native_xy", (sx, sy))
                layer_zdim = zdim - zdim_before
                stretched = (do_match and not
                             (np.isclose(native_xi, sx, rtol=1e-2) and
                              np.isclose(native_yi, sy, rtol=1e-2)))
                extras = ([f"native XY: {native_xi:.3f} x {native_yi:.3f} Å"]
                          if stretched else [])
                _log_layer_result(len(slab_u.atoms), (sx, sy, sz),
                                  zdim, layer_zdim, extra_lines=extras)

            elif ltype == "solvent":
                solv_names = " + ".join(
                    getattr(s, "resname", "?") for s in layer["solvent"]
                ) or "ions"
                log_subheader(logger, f"Layer {tag}")
                logger.info("  ├> solvent: %s", solv_names)
                if layer["density"] is not None:
                    logger.info("  ├> density: %.2f g/cm³", layer["density"])
                solv_box = self._build_solvent_layer(layer, xsize, ysize, all_sp_univs)
                if solv_box is not None:
                    zdim_before = zdim
                    system, zdim = add_component(system, solv_box, zdim, padding=padding)
                    layer_zdim = zdim - zdim_before
                    svx, svy, svz = solv_box.dimensions[:3]
                    mol_counts = Counter(res.resname for res in solv_box.residues)
                    # one line listing every species (solvent + solute) with counts
                    all_sp = list(layer["solvent"]) + list(layer["solute"])
                    mol_line = ",  ".join(
                        f"{getattr(s, 'resname', '?')} ({mol_counts.get(getattr(s, 'resname', '?'), 0)} mol)"
                        for s in all_sp
                    )
                    if mol_line:
                        logger.info("  ├> %s", mol_line)
                    _log_layer_result(len(solv_box.atoms), (svx, svy, svz),
                                      zdim, layer_zdim)
                else:
                    logger.warning("  └─> empty; check packmol.log")
                    system, zdim = add_component(system, solv_box, zdim, padding=padding)

            elif ltype == "vacuum":
                layer_zdim = layer["zdim"]
                zdim += layer_zdim
                log_subheader(logger, f"Layer {tag}")
                logger.info("  ├> vacuum: %.1f Å", layer_zdim)
                _log_layer_result(None, None, zdim, layer_zdim)

            if ii == 0:
                first_layer_zdim = zdim

        return system, zdim, first_layer_zdim

    def _build_solvent_layer(
        self,
        layer: dict,
        xsize: float,
        ysize: float,
        all_sp_univs: list,
    ) -> Optional[mda.Universe]:
        """
        Build one solvent layer via PACKMOL.

        Applies dilation if requested and delegates to :func:`make_solvent_box`.
        """
        dilate   = layer["dilate"]
        eff_zdim = layer["zdim"] * dilate
        eff_rho  = layer["density"] / dilate if layer["density"] is not None else None

        if dilate != 1.0:
            logger.info(
                "  ├> dilation x%.2f: packing %.1f Å at %.3f g/cm³  (target: %.1f Å)",
                dilate, eff_zdim,
                eff_rho if eff_rho is not None else float("nan"),
                layer["zdim"],
            )

        return make_solvent_box(
            species=all_sp_univs,
            solvent=layer["solvent"] or None,
            solute=layer["solute"] or None,
            volume=[xsize, ysize, eff_zdim],
            density=eff_rho,
            nsolute=layer["nsolute"],
            concentration=layer["concentration"],
            conmodel=layer["conmodel"],
            solute_pos=layer["solute_pos"],
            nsolvent=layer["nsolvent"],
            tolerance=layer["packmol_tolerance"],
            ratio=layer["ratio"],
        )

    @staticmethod
    def _resolve_match_cell(match_cell):
        """
        Parse the *match_cell* parameter used in :meth:`build`.

        Returns
        -------
        cell_ref : Specie or None
            The reference species whose XY cell locks *xsize/ysize*, or
            ``None`` when a plain bool was supplied.
        do_match : bool
            Whether ``to_universe`` should apply XY scaling.
        """
        if match_cell is False or match_cell is True:
            return None, bool(match_cell)
        return match_cell, True

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

    @staticmethod
    def _apply_stack_axis(
        system: mda.Universe,
        xsize: float,
        ysize: float,
        zdim: float,
        stack_axis: str,
    ):
        """
        Permute atomic coordinates so that the stacking direction is *stack_axis*.

        Assembly always builds along Z.  This method applies a lossless
        coordinate permutation and updates the cell dimensions accordingly:

        - ``"x"`` : (x, y, z) → (z, y, x)  |  cell [zdim, ysize, xsize]
        - ``"y"`` : (x, y, z) → (x, z, y)  |  cell [xsize, zdim, ysize]

        Returns the updated (system, new_xsize, new_ysize).
        """
        pos = system.atoms.positions
        if stack_axis == "x":
            system.atoms.positions = pos[:, [2, 1, 0]]
            system.dimensions = [zdim, ysize, xsize, 90, 90, 90]
            new_xsize, new_ysize = zdim, ysize
        else:  # "y"
            system.atoms.positions = pos[:, [0, 2, 1]]
            system.dimensions = [xsize, zdim, ysize, 90, 90, 90]
            new_xsize, new_ysize = xsize, zdim
        logger.info("  └─> axis permuted Z -> %s", stack_axis.upper())
        return system, new_xsize, new_ysize

    def _register(self, *species: Any) -> None:
        """Add species to the global registry if not already present."""
        for sp in species:
            if sp is not None and not any(s is sp for s in self._all_species):
                self._all_species.append(sp)

    def _update_topology_indexes(self) -> None:
        n_species = len(self._all_species)
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
            "  └─> %d species,  %d atom types,  %d bond,  %d angle,  %d dihedral,  %d improper",
            n_species, len(atom_types), len(nitems["_btype"]),
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


# ---------------------------------------------------------------------------
# Backwards-compatibility alias
# ---------------------------------------------------------------------------

class BoxBuilder(SimCell):
    """Deprecated alias for :class:`SimCell`.  Will be removed in a future version."""

    def __init__(self, *args, **kwargs):
        import warnings
        warnings.warn(
            "BoxBuilder is deprecated and will be removed in a future version. "
            "Use SimCell instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)
