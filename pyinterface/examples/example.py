#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script demonstrating how to use the pyinterface package to build a simulation box.

Author: Roncofaber
"""

import ase
from pyinterface import SimulationBox
from pyinterface.core.specie import Specie
from pyinterface.core.topology import Bond, Angle

#%% Set up simulation box

# Solvent setup (e.g., TIP3P water model)
# https://docs.lammps.org/Howto_tip3p.html (Ewald model, check Price paper)
b1 = Bond("O", "H", kr=450, r0=0.9572)
a1 = Angle("H", "O", "H", kr=55, theta0=104.52)

mol = Specie("H2O", charges=[-0.83, 0.415, 0.415], bonds=b1, angles=a1,
             lj={"O": [0.102, 3.188], "H": [0.0, 0.0]})

# Ions setup
# https://www.sciencedirect.com/science/article/pii/S0167732216322760
Na = Specie("Na", charges=+0.8, lj={"Na": [0.00280, 3.3304]})
Cl = Specie("Cl", charges=-0.8, lj={"Cl": [0.1178, 4.41720]})

# Interface setup
interface = Specie("gold.pdb", charges=0.0, lj={"Au": [5.29, 2.62904212]})

#%% Create a SimulationBox instance

simobj = SimulationBox(
    solvent   = mol,
    solute    = [Na, Cl],
    interface = interface,
    # enderface = interface
)

#%% Define simulation box parameters

solvent_vol = [15, 15, 25]  # Size of solvent box in Å
solvent_rho = 0.99713       # Density in g/cm³
nions       = 1             # Number of ions
nlayers     = 1             # Number of interface layers

#%% Build the simulation box

system = simobj.make_simulation_box(
    solvent_vol      = solvent_vol,
    solvent_rho      = solvent_rho,
    nions            = nions,
    # concentration  = 0.35,  # Molar concentration
    # conmodel       = cm,
    layers           = nlayers,
    to_ase           = True,
    write_data       = True,
    # padding       = -1.5
    center_electrode = False,
    ion_pos          = "left"
)
