#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SimCell example: a CO2 gas pocket embedded in a bulk water solvent.

Demonstrates the Region-based `regions=` parameter on add_solvent(): a
Sphere region carved out of the bulk water layer and filled with a
different species and density.

Author: roncofaber
"""

from mdinterface import SimCell
from mdinterface.database import Water
from mdinterface.core.specie import Specie
from mdinterface.build.regions import Sphere

#%% Define species

water = Water()
co2 = Specie("CO2")

#%% Set up simulation box

simbox = SimCell(xysize=[25, 25])

#%% Define the pocket: a CO2-filled sphere embedded in the water layer

bubble = Sphere(center=(12.5, 12.5, 15.0), radius=6.0)

simbox.add_solvent(
    water,
    zdim=30,
    density=1.0,
    regions=[bubble.fill(co2, density=0.8)],
)

#%% Build

simbox.build(padding=0.5)

#%% Output

atoms = simbox.to_ase()
universe = simbox.universe
