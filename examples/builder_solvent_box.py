#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BoxBuilder example: pure solvent box with dissolved species.

Equivalent to make_solvent_box.py but using the fluent BoxBuilder API.

Author: roncofaber
"""

from mdinterface import BoxBuilder
from mdinterface.database import Water
from mdinterface.core.specie import Specie

#%% Define species

# Solvent: TIP3P water (Ewald model)
water = Water(model="ewald")

# Solute: ammonia with LigParGen parameters
amm = Specie("NH3", ligpargen=True)

#%% Build with BoxBuilder

system = (
    BoxBuilder(xysize=[20, 20])
        .add_solvent(
            water,
            ions=[amm],
            nions=15,      # 15 ammonia molecules
            zdim=20,       # 20 Å thick region
            density=1.0,   # water density in g/cm³
        )
        .build(padding=0.5)
        .write_lammps("data.lammps", atom_style="full", write_coeff=True)
)

# The MDAnalysis Universe is available if you need it for further processing
universe = system.universe
