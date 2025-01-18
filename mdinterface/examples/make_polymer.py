#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script demonstrating how to use the mdinterface package to build a polymer.

Author: roncofaber
"""

from mdinterface import SimulationBox
from mdinterface.core.specie import Specie
from mdinterface.database import Water, Graphene
from mdinterface.build.polymer import build_polymer
from mdinterface.io.read import read_lammps_data_file

#%% make polymer

# load an existing monomer (in this case, we already have force field params,
# so we take advantage and read them)
monomer, atoms, bonds, angles, dihedrals, impropers = read_lammps_data_file("data/nafion_n1.data.lammps")

# define the ends that will be connected when making a polymer
# you can also set start_end_idxs=[34, 52] as build_polymer option
symbols = monomer.get_chemical_symbols()
symbols[34] = "X"
symbols[52] = "X"
monomer.set_chemical_symbols(symbols)

# build a polymer with 5 repeating units
nafion = build_polymer(monomer, "F", 5)

#%% create species with topology

nas = Specie(atoms=nafion, name="NAFI", bonds=bonds, angles=angles,
             dihedrals=dihedrals, fix_missing=True)
mol = Water()


#%% make simulation box

simobj = SimulationBox(
    # solvent   = mol,
    solute    = [nas],#O2
    # interface = Graphene(),
    # enderface = enderface
    )

#%% populate system.
# Here the box is initially huge cause we need to relax the polymer carefully...

solvent_vol = [100, 100, 100] # size of solvent \AA
solvent_rho = 0.15            # density g/cm3
nions       = 12              # number of ions
nlayers     = 2               # number of interface layers

system = simobj.make_simulation_box(solvent_vol,
                                    solvent_rho,
                                    nions = nions,
                                    # concentration = 0.35, #molar conc.
                                    # conmodel = cm,
                                    # layers=nlayers, 
                                    to_ase = True,
                                    write_data = True,
                                    ion_pos="box",
                                    # padding=2.0,
                                    # center_electrode = True,
                                    # vacuum=20,\
                                    )
