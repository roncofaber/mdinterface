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
import MDAnalysis as mda

import numpy as np
import subprocess

logger = logging.getLogger("mdinterface.box")


def populate_box(
    volume: List[float],
    instructions: List[Tuple[Any, Union[int, List[float]], str]],
    input_file: str = "input_packmol.in",
    output_file: str = "system.pdb",
    tolerance: float = 2.0
) -> Optional[mda.Universe]:

    if not instructions:
        return None

    # check volume
    assert len(volume) == 3, "Check volume!"
    
    # generate box boundaries with 1 AA padding
    box = np.concatenate(([1,1,1], np.asarray(volume)-1)).tolist()
    
    tmp_files = ["packmol.log", "input_packmol.in", "system.pdb"]
    with open(input_file, "w") as fout:
        
        fout.write(header.format(tolerance, output_file, np.random.randint(100000)))
        
        for cc, instruction in enumerate(instructions):
            
            # unpack instructions
            mol = instruction[0]
            rep = instruction[1]
            typ = instruction[2]
            
            if isinstance(rep, int):
                if not rep:
                    continue
            
            if typ == "box": # normal add
                fout.write(box_place.format(cc, rep, " ".join(map(str, box))))
            
            elif typ == "fixed": # coordinate -> fixed point
                fout.write(fix_place.format(cc, *rep))
                
            elif typ == "zfixed": # small bin to use in CM
                tbox = box.copy()
                tbox[2] = tbox[2] - mol.estimate_specie_radius()
                tbox[5] = tbox[5] + mol.estimate_specie_radius()
                fout.write(box_place.format(cc, 1, " ".join(map(str, tbox))))
                mol = mol.to_universe()
            
            else:
                raise ValueError("Wrong instructions")
            
            # write tmp pdb file and store info
            mol.atoms.write("mol_{}.pdb".format(cc))
            tmp_files.append("mol_{}.pdb".format(cc))
            
    # run packmol
    logger.info("  Running PACKMOL (tolerance=%.1f A, %d molecule type(s))",
                tolerance, len(instructions))
    try:
        with open(input_file, 'r') as stdin_f, open('packmol.log', 'w') as stdout_f:
            subprocess.run(['packmol'], stdin=stdin_f, stdout=stdout_f, check=True)
        logger.debug("    PACKMOL converged")

    except subprocess.CalledProcessError:
        logger.warning("  PACKMOL may not have converged - check packmol.log")

    try:
        universe = mda.Universe(output_file)
        logger.debug("    PACKMOL output: %d atoms", len(universe.atoms))
    except Exception:
        logger.warning("  Could not load PACKMOL output file '%s'", output_file)
        universe = None

    # remove temp mol files and packmol files
    subprocess.call(['rm'] + tmp_files)
    
    return universe

# generate a slab from a unit cell
def make_interface_slab(interface_uc, xsize, ysize, layers=1):
    
    if layers == 0 or interface_uc is None:
        return None
    
    xrep = int(np.round(xsize/interface_uc.atoms.get_cell()[0][0]))
    yrep = int(np.round(ysize/interface_uc.atoms.get_cell()[1][1]))
    
    slab = interface_uc.copy()
    
    if not np.isclose(np.dot(slab.atoms.cell[0], [1,0,0]), slab.atoms.cell[0][0]):
        xrep +=1
        print("WARNING: check interface if pattern matches")
    
    if not np.isclose(np.dot(slab.atoms.cell[1], [0,1,0]), slab.atoms.cell[1][1]):
        yrep +=1
        print("WARNING: check interface if pattern matches")
    
    if not (xrep, yrep, 1) == (1,1,1):
        slab.repeat((xrep, yrep, 1), make_cubic=True)
    
    if layers > 1: # helps with indexing
        slab.repeat([1,1,layers])
    
    slab.atoms.center()
    # slab.atoms.rattle()
    
    return slab

# add a component to the system
def add_component(system, component, zdim, padding=0):
    
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
