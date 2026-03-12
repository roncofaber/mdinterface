#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SimCell example: pure solvent box with dissolved species.

Equivalent to make_solvent_box.py but using the fluent SimCell API.

Author: roncofaber
"""

from mdinterface import SimCell
from mdinterface.database import Water
from mdinterface.core.specie import Specie

#%% Define species

# Solvent: TIP3P water (Ewald model)
water = Water(model="ewald")

# Solute: ammonia with LigParGen parameters
amm = Specie("NH3", ligpargen=True)

#%% Set up simulation box

simbox = SimCell(xysize=[20, 20])

#%% Add layers

simbox.add_solvent(
    water,
    solute=[amm],
    nsolute=15,      # 15 ammonia molecules
    zdim=20,       # 20 Å thick region
    density=1.0,   # water density in g/cm³
)

#%% Build

simbox.build(padding=0.5)

#%% Output

# Write LAMMPS data file (with force-field coefficients)
simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True)

# Convert to ASE Atoms (e.g. for visualisation or further manipulation)
atoms = simbox.to_ase()

# Access the raw MDAnalysis Universe if needed
universe = simbox.universe
