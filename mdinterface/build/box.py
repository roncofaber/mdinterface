#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
box.py — low-level box assembly utilities.

PACKMOL orchestration (``populate_box``), slab tiling (``make_interface_slab``),
and layer stacking (``add_component``).  Solvent-specific logic lives in
``solvent.py``.
"""

from typing import List, Optional, Union, Tuple, Any
from mdinterface.io.packmol import header, box_place, fix_place, structure_place

import logging
import os
import shutil
import tempfile
import MDAnalysis as mda
import ase.io
from ase import Atoms

import numpy as np
import subprocess

logger = logging.getLogger(__name__)


class PackmolError(RuntimeError):
    """PACKMOL execution or output-processing failure."""

    def __init__(self, message, tempdir, log_path, returncode=None):
        self.tempdir = tempdir
        self.log_path = log_path
        self.returncode = returncode
        details = [message]
        if returncode is not None:
            details.append(f"return code: {returncode}")
        details.extend([f"temporary files: {tempdir}", f"log: {log_path}"])
        super().__init__("; ".join(details))


def _indent_constraints(lines):
    """Join PACKMOL constraint lines, each indented for the .in file."""
    return "\n".join(f"    {line}" for line in lines)


def _write_packmol_template(specie, path):
    """Write a Specie's positions/resname to a PDB template for a PACKMOL structure block."""
    atoms = specie.atoms.copy()
    atoms.set_array("residuenames", np.array([specie.resname] * len(atoms), dtype=object))
    ase.io.write(path, atoms, format="proteindatabank")


def populate_box(
    volume: List[float],
    instructions: List[Tuple[Any, Union[int, List[float]], str]],
    tolerance: float = 2.0
) -> Optional[Atoms]:
    """
    Run PACKMOL to pack molecules into a box and return the result as an ASE Atoms object.

    Each instruction tuple describes one molecule type to place:

    - ``(specie, count, "box")`` -- pack *count* copies anywhere in the box.
    - ``(specie, count, "box", bounds)`` -- restrict placement to a
      sub-region ``[xmin, ymin, zmin, xmax, ymax, zmax]``.
    - ``(specie, count, "box", None, extra_constraints)`` -- pack *count* copies
      in the full box with additional PACKMOL constraint lines (e.g., to exclude a region).
      *extra_constraints* is a list of raw PACKMOL constraint strings.
    - ``(specie, count, "region", region_obj)`` -- pack *count* copies constrained to
      a Region object's interior (using its ``packmol_line("inside")`` method).
    - ``(specie, count, "region", region_obj, extra_constraints)`` -- as above, plus
      additional raw PACKMOL constraint lines (e.g. to exclude a nested sub-region).
    - ``(specie, coords, "fixed")`` -- place a single molecule at fixed
      fractional coordinates *coords*.
    - ``(specie, coords, "zfixed")`` -- place a single molecule at fixed
      fractional XY coordinates with Z placement in a thin bin near the center of mass.

    Parameters
    ----------
    volume : list of float
        Box dimensions ``[x, y, z]`` in Angstroms.
    instructions : list of tuple
        Packing instructions (see above).  An empty list returns ``None``
        immediately without calling PACKMOL.
    tolerance : float, default 2.0
        Minimum intermolecular distance in Angstroms passed to PACKMOL.

    Returns
    -------
    ase.Atoms or None
        The packed system with ``residuenames``/``residuenumbers`` arrays set,
        or ``None`` if the instructions list is empty.

    Raises
    ------
    TypeError
        If *volume* is not a numeric three-element sequence.
    ValueError
        If *volume* contains non-finite or non-positive dimensions.
    PackmolError
        If PACKMOL cannot be started, exits unsuccessfully, or produces
        missing or unreadable output. Temporary files are retained and their
        location is included in the exception.
    """
    if not instructions:
        return None

    if not isinstance(volume, (list, tuple, np.ndarray)):
        raise TypeError("volume must be a list, tuple, or numpy array with three dimensions.")
    if len(volume) != 3:
        raise ValueError(f"volume must contain exactly three dimensions, got {len(volume)}.")
    try:
        volume = np.asarray(volume, dtype=float)
    except (TypeError, ValueError) as exc:
        raise TypeError("volume dimensions must be numeric.") from exc
    if not np.all(np.isfinite(volume)) or np.any(volume <= 0):
        raise ValueError(f"volume dimensions must be finite and positive, got {volume.tolist()}.")

    # generate box boundaries with 1 AA padding
    box = np.concatenate(([1,1,1], volume-1)).tolist()

    # all PACKMOL files go in a temp dir; kept on failure for inspection
    tmpdir = tempfile.mkdtemp(prefix="packmol_")
    input_file  = os.path.join(tmpdir, "input_packmol.in")
    output_file = os.path.join(tmpdir, "system.pdb")
    log_file    = os.path.join(tmpdir, "packmol.log")

    try:
        with open(input_file, "w") as fout:
            fout.write(header.format(tolerance, output_file, np.random.randint(100000)))

            for cc, instruction in enumerate(instructions):

                # unpack instructions
                mol = instruction[0]
                rep = instruction[1]
                typ = instruction[2]
                custom_bounds = instruction[3] if len(instruction) > 3 else None

                if isinstance(rep, int):
                    if not rep:
                        continue

                if typ == "box":  # normal add, optionally with extra "outside" constraints
                    bounds = custom_bounds if custom_bounds is not None else box
                    extra_constraints = instruction[4] if len(instruction) > 4 else None
                    lines = ["inside box " + " ".join(map(str, bounds))]
                    if extra_constraints:
                        lines.extend(extra_constraints)
                    fout.write(structure_place.format(cc, rep, _indent_constraints(lines)))

                elif typ == "region":  # placed via a Region's own "inside" constraint,
                                       # optionally excluding nested sub-regions
                    region_obj = instruction[3]
                    extra_constraints = instruction[4] if len(instruction) > 4 else None
                    lines = [region_obj.packmol_line("inside")]
                    if extra_constraints:
                        lines.extend(extra_constraints)
                    fout.write(structure_place.format(cc, rep, _indent_constraints(lines)))

                elif typ == "fixed":  # coordinate -> fixed point
                    fout.write(fix_place.format(cc, *rep))

                elif typ == "zfixed":  # small bin to use in CM
                    tbox = box.copy()
                    tbox[2] = tbox[2] - mol.estimate_specie_radius()
                    tbox[5] = tbox[5] + mol.estimate_specie_radius()
                    lines = ["inside box " + " ".join(map(str, tbox))]
                    fout.write(structure_place.format(cc, 1, _indent_constraints(lines)))

                else:
                    raise ValueError("Wrong instructions")

                _write_packmol_template(mol, os.path.join(tmpdir, "mol_{}.pdb".format(cc)))

        # run packmol
        logger.debug("  >> Running PACKMOL (tolerance=%.1f Å, %d molecule type(s))",
                     tolerance, len(instructions))
        with open(input_file, "r") as stdin_f, open(log_file, "w") as stdout_f:
            try:
                result = subprocess.run(
                    ["packmol"],
                    stdin=stdin_f,
                    stdout=stdout_f,
                    stderr=subprocess.STDOUT,
                    cwd=tmpdir,
                )
            except FileNotFoundError as exc:
                raise PackmolError(
                    "PACKMOL executable was not found on PATH",
                    tmpdir,
                    log_file,
                ) from exc

        if result.returncode != 0:
            raise PackmolError(
                "PACKMOL exited unsuccessfully",
                tmpdir,
                log_file,
                returncode=result.returncode,
            )

        logger.debug("  >> PACKMOL converged")

        if not os.path.isfile(output_file):
            raise PackmolError(
                f"PACKMOL did not create the expected output file {output_file}",
                tmpdir,
                log_file,
                returncode=result.returncode,
            )

        try:
            packed = ase.io.read(output_file, format="proteindatabank")
            logger.debug("  >> PACKMOL output: %d atoms", len(packed))
        except Exception as exc:
            raise PackmolError(
                f"PACKMOL output could not be read from {output_file}",
                tmpdir,
                log_file,
                returncode=result.returncode,
            ) from exc

        # success -- clean up
        shutil.rmtree(tmpdir, ignore_errors=True)
        return packed

    except Exception:
        logger.error("  >> PACKMOL failed; temporary files kept at: %s", tmpdir)
        raise

def make_interface_slab(interface_uc, xsize, ysize, layers=1):
    """
    Tile a unit-cell Specie into a surface slab of the requested XY size.

    The unit cell is replicated in X and Y using the smallest integer repeat
    counts that cover *xsize* x *ysize* Angstroms, then repeated *layers*
    times along Z.

    Parameters
    ----------
    interface_uc : Specie
        Unit-cell species to tile.
    xsize : float
        Target X dimension in Angstroms.
    ysize : float
        Target Y dimension in Angstroms.
    layers : int, default 1
        Number of Z repetitions (i.e. atomic layers along the stacking axis).

    Returns
    -------
    Specie or None
        Tiled slab Specie, or ``None`` if *layers* is 0 or *interface_uc*
        is ``None``.
    """
    if layers == 0 or interface_uc is None:
        return None

    if xsize <= 0 or ysize <= 0:
        raise ValueError(f"Slab target dimensions must be positive, got [{xsize}, {ysize}].")

    cell = interface_uc.atoms.get_cell()
    xperiod = float(cell[0][0])
    yperiod = float(cell[1][1])
    if xperiod <= 0 or yperiod <= 0:
        raise ValueError("Slab cell must have positive X and Y lattice-vector projections.")
    if not np.allclose(cell, np.diag(np.diag(cell))):
        logger.warning("Non-orthogonal slab cell will be converted to an orthogonal cell during tiling.")

    xrep = max(1, int(np.ceil(xsize / xperiod)))
    yrep = max(1, int(np.ceil(ysize / yperiod)))
    
    slab = interface_uc.copy()
    
    if not (xrep, yrep, 1) == (1,1,1):
        slab.repeat((xrep, yrep, 1), make_cubic=True)
    
    if layers > 1: # helps with indexing
        slab.repeat([1,1,layers])
    
    slab.atoms.center()
    # slab.atoms.rattle()
    
    return slab

def add_component(system, component, zdim, padding=0):
    """
    Translate and merge a new component Universe into a growing system Universe.

    The component is shifted so that its bottom face sits at the current Z
    cursor (*zdim* + *padding*), then merged with *system*.  The Z cursor is
    advanced by the component's Z extent.

    Parameters
    ----------
    system : mda.Universe or None
        The assembly built so far.  ``None`` means the component becomes the
        first element.
    component : mda.Universe or None
        The layer to append.  ``None`` is a no-op (returns *system* and
        *zdim* unchanged).
    zdim : float
        Current Z height of the assembly in Angstroms (bottom of the new layer).
    padding : float, default 0
        Extra gap in Angstroms inserted below the new component.

    Returns
    -------
    system : mda.Universe
        Updated assembly with the component merged in.
    zdim : float
        New Z height after appending the component.
    """
    # nothing to add here
    if component is None:
        return system, zdim
    
    # ohh, let's lego the shit out of this
    component = component.copy()
    
    # component: "look at me, I am the system now."
    if system is None:
        component.atoms.translate([0, 0, zdim])
        system = component
        zdim += component.dimensions[2]
    
    # make space and add it to the pile
    else:
        component.atoms.translate([0, 0, zdim + padding])
        system = mda.Merge(system.atoms, component.atoms)
        zdim += component.dimensions[2] + padding
    return system, zdim
