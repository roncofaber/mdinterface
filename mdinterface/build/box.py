#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
box.py — low-level box assembly utilities.

PACKMOL orchestration (``populate_box``), slab tiling (``make_interface_slab``),
and layer stacking (``add_component``).  Solvent-specific logic lives in
``solvent.py``.
"""

from typing import List, Optional, Union, Tuple, Any
from mdinterface.io.packmol import header, box_place, fix_place

import logging
import os
import shutil
import tempfile
import MDAnalysis as mda

import numpy as np
import subprocess

logger = logging.getLogger(__name__)


def populate_box(
    volume: List[float],
    instructions: List[Tuple[Any, Union[int, List[float]], str]],
    tolerance: float = 2.0
) -> Optional[mda.Universe]:
    """
    Run PACKMOL to pack molecules into a box and return the result as a Universe.

    Each instruction tuple describes one molecule type to place:

    - ``(universe, count, "box")`` -- pack *count* copies anywhere in the box.
    - ``(universe, count, "box", bounds)`` -- restrict placement to a
      sub-region ``[xmin, ymin, zmin, xmax, ymax, zmax]``.
    - ``(universe, coords, "fixed")`` -- place a single molecule at fixed
      fractional coordinates *coords*.

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
    mda.Universe or None
        The packed system as an MDAnalysis Universe, or ``None`` if the
        instructions list is empty or PACKMOL fails to write an output file.
    """
    if not instructions:
        return None

    # check volume
    assert len(volume) == 3, "Check volume!"

    # generate box boundaries with 1 AA padding
    box = np.concatenate(([1,1,1], np.asarray(volume)-1)).tolist()

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

                if typ == "box":  # normal add
                    bounds = custom_bounds if custom_bounds is not None else box
                    fout.write(box_place.format(cc, rep, " ".join(map(str, bounds))))

                elif typ == "fixed":  # coordinate -> fixed point
                    fout.write(fix_place.format(cc, *rep))

                elif typ == "zfixed":  # small bin to use in CM
                    tbox = box.copy()
                    tbox[2] = tbox[2] - mol.estimate_specie_radius()
                    tbox[5] = tbox[5] + mol.estimate_specie_radius()
                    fout.write(box_place.format(cc, 1, " ".join(map(str, tbox))))
                    mol = mol.to_universe()

                else:
                    raise ValueError("Wrong instructions")

                mol.atoms.write(os.path.join(tmpdir, "mol_{}.pdb".format(cc)))

        # run packmol
        logger.debug("  ├> Running PACKMOL (tolerance=%.1f Å, %d molecule type(s))",
                     tolerance, len(instructions))
        with open(input_file, "r") as stdin_f, open(log_file, "w") as stdout_f:
            result = subprocess.run(["packmol"], stdin=stdin_f, stdout=stdout_f, cwd=tmpdir)

        if result.returncode != 0:
            logger.warning("  ├> PACKMOL may not have converged - check %s", tmpdir)
        else:
            logger.debug("  ├> PACKMOL converged")

        try:
            universe = mda.Universe(output_file)
            logger.debug("  └─> PACKMOL output: %d atoms", len(universe.atoms))
        except Exception:
            logger.warning("  └─> Could not load PACKMOL output; temp files kept at: %s", tmpdir)
            return None

        # success -- clean up
        shutil.rmtree(tmpdir, ignore_errors=True)
        return universe

    except Exception:
        logger.warning("  └─> PACKMOL failed; temp files kept at: %s", tmpdir)
        raise

def make_interface_slab(interface_uc, xsize, ysize, layers=1):
    """
    Tile a unit-cell Specie into a surface slab of the requested XY size.

    The unit cell is replicated in X and Y using the nearest integer repeat
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
    
    xrep = int(np.round(xsize/interface_uc.atoms.get_cell()[0][0]))
    yrep = int(np.round(ysize/interface_uc.atoms.get_cell()[1][1]))
    
    slab = interface_uc.copy()
    
    if not np.isclose(np.dot(slab.atoms.cell[0], [1,0,0]), slab.atoms.cell[0][0]):
        xrep += 1
        logger.warning("Non-orthogonal cell along X ; check interface pattern")

    if not np.isclose(np.dot(slab.atoms.cell[1], [0,1,0]), slab.atoms.cell[1][1]):
        yrep += 1
        logger.warning("Non-orthogonal cell along Y ; check interface pattern")
    
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
