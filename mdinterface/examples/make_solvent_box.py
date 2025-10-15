#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script demonstrating how to use the mdinterface package to build a simulation box.

Author: roncofaber
"""

from mdinterface import SimulationBox
from mdinterface.core.specie import Specie
from mdinterface.database import Water

#%% Set up simulation box

# Solvent setup (e.g., TIP3P water model with Ewald correction)
wat = Water(model="ewald")

# add extra specie (e.g. ammonia) with LigParGen parameters
amm = Specie("NH3", ligpargen=True)

#%% Create a SimulationBox instance

# here we set the water as solvent
simobj = SimulationBox(
    solvent   = wat,
    solute    = [amm],
)

#%% Define simulation box parameters

xysize = [20, 20] # XY cross section

# add subsequent layers of building blocks
layering = [
    # add a solvent layer
    # nspecies (preferred) or nions (deprecated) can be:
    # - int: same number for all ion types, e.g., nspecies: 1
    # - list: different numbers for each ion type, e.g., nspecies: [1, 2]
    # - None: no ions (solvent only)
    {"type": "solvent", "rho": 1.0, "zdim": 20, "nspecies": 15},
    ]
# water density is set to 1 g/cm3, yet it should probably be lower if you have solute

# Alternative examples:
# For concentration-based approach (in molar):
# {"type": "solvent", "rho": 1.0, "zdim": 25, "concentration": 0.1}

# For specifying exact number of solvent molecules:
# {"type": "solvent", "nsolvent": 500, "zdim": 25, "nspecies": 10}

#%% Build the simulation box

# main function to build a box given the instructions defined above
system = simobj.make_simulation_box(
    xysize, # XY cross section
    layering, # layering info
    padding = 0.5, # pad this amount between any layer (Angstrom)
    to_ase = True, # return an ase.Atoms object, otherwishe mda.Universe
    write_data = True, # write LAMMPS data file
    filename = "data.lammps", # name of the LAMMPS data file
    center_electrode = False, # shifts everything by 50% along Z to place the first interface in the center (it's better to layer instead)
    layered = False, # assign different molecule indexes to each layer of interface/enderface type slabs for LAMMPS
    hijack = None, # give an ase.Atoms to override the positions (ATTENTION: should have exactly same order, use with CARE (or just don't...))
    match_cell = False, # if trying to mix slabs with different cross section, use this to deform them so they match XY, use with care
    atom_style = "full", # write LAMMPS data file with the format 'full' or 'atomic'
    write_coeff = True   # write or not coefficients in the LAMMPS data file
    )