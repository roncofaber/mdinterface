#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script demonstrating how to use the mdinterface package to build a polymer.

Author: roncofaber
"""

import numpy as np
from mdinterface import SimulationBox
from mdinterface.core.polymer import Polymer
from mdinterface.database import Water, Hydronium
from mdinterface.io.read import read_lammps_data_file

#%% make polymer

# load an existing monomer (in this case, we already have force field params,
# so we take advantage and read them)
monomer, atoms, bonds, angles, dihedrals, impropers = read_lammps_data_file("data/nafion_n1_dep.data.lammps")

# add nominal charge (important for charged molecules (e.g. ionomer))
monomer.set_array("nominal_charge", np.array(len(monomer)*[0]))
monomer.arrays["nominal_charge"][69] = -1

# define number of repeating units
nrep = 5

# make the polymer specie. You need LigParGen to be configured to make this work!
# if not, set refine_charges=False
nafn = Polymer(atoms=monomer, name="NAFN", bonds=bonds, angles=angles,
             dihedrals=dihedrals, fix_missing=True, start_end_idxs=[34,52],
             substitute="F", nrep=nrep, refine_charges=True, offset=True)

nafn.round_charges(7)

#%% load some other species

mol = Water()
hyd = Hydronium()

#%%
simobj = SimulationBox(
    # solvent   = mol,
    solute    = [nafn, hyd, mol],#O2
    # interface = Graphene(),
    # enderface = enderface
    )

#%% Define simulation box parameters

xysize = [250, 250] # XY cross section

# add subsequent layers of building blocks
layering = [
    # add a solvent slab, nions can be a list with the number of ions as provided above in solute
    {"type": "solvent", "zdim": 250, "nions": [10, nrep*10, 200], "ion_pos": 'box'},
    ]

#%% Build the simulation box

# main function to build a box given the instructions defined above
system = simobj.make_simulation_box(
    xysize, # XY cross section
    layering, # layering info
    padding = 0.5, # pad this amount between any layer (Angstrom)
    to_ase = True, # return an ase.Atoms object, otherwishe mda.Universe
    write_data = False, # write LAMMPS data file
    filename = "data.lammps", # name of the LAMMPS data file
    center_electrode = False, # shifts everything by 50% along Z to place the first interface in the center (it's better to layer instead)
    layered = False, # assign different molecule indexes to each layer of interface/enderface type slabs for LAMMPS
    hijack = None, # give an ase.Atoms to override the positions (ATTENTION: should have exactly same order, use with CARE (or just don't...))
    match_cell = False, # if trying to mix slabs with different cross section, use this to deform them so they match XY, use with care
    remove_charges = False # do not write charges in LAMMPS data file
    )
