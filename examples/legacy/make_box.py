#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script demonstrating how to use the mdinterface package to build a simulation box.

Author: roncofaber
"""

from mdinterface import SimulationBox
from mdinterface.core.specie import Specie
from mdinterface.core.topology import Bond, Angle
from mdinterface.database import Water, Metal111

#%% Set up simulation box

# Solvent setup (e.g., TIP3P water model)
# https://docs.lammps.org/Howto_tip3p.html (Ewald model, check Price paper)
b1 = Bond("O", "H", kr=450, r0=0.9572)
a1 = Angle("H", "O", "H", kr=55, theta0=104.52)

mol = Specie("H2O", charges=[-0.83, 0.415, 0.415], bonds=b1, angles=a1,
             lj={"O": [0.102, 3.188], "H": [0.0, 0.0]})

# alternatively you can load water from the database as:
wat = Water()

# Ions setup
# https://www.sciencedirect.com/science/article/pii/S0167732216322760
Na = Specie("Na", charges=+0.8, lj={"Na": [0.00280, 3.3304]})
Cl = Specie("Cl", charges=-0.8, lj={"Cl": [0.1178, 4.41720]})

# Interface setup
interface = Metal111("Au")

#%% Create a SimulationBox instance

simobj = SimulationBox(
    solvent   = mol,
    solute    = [Na, Cl],
    interface = interface,
    enderface = interface
)

#%% Define simulation box parameters

xysize = [15, 15] # XY cross section

# add subsequent layers of building blocks
layering = [
    # add the interface type slab
    {"type": "interface", "nlayers" : 1 },
    # add a solvent slab, nions can be a list with the number of ions as provided above in solute
    {"type": "solvent", "rho": 1.0, "zdim": 25, "nions": 1},
    # add an enderface type slab
    {"type": "enderface", "nlayers" : 1 },
    # add vacuum, if you want
    {"type": "vacuum", "zdim": 5},
    ]

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
