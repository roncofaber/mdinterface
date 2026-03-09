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

# ---------------------------------------------------------------------------
# 1. Define species
# ---------------------------------------------------------------------------

# Solvent: TIP3P water (Ewald model)
water = Water(model="ewald")

# Solute: ammonia with LigParGen parameters
amm = Specie("NH3", ligpargen=True)

# ---------------------------------------------------------------------------
# 2. Assemble layers
# ---------------------------------------------------------------------------

builder = (
    BoxBuilder(xysize=[20, 20])
        .add_solvent(
            water,
            ions=[amm],
            nions=15,      # 15 ammonia molecules
            zdim=20,       # 20 Å thick region
            density=1.0,   # water density in g/cm³
        )
)

# ---------------------------------------------------------------------------
# 3. Build
# ---------------------------------------------------------------------------

builder.build(padding=0.5)

# ---------------------------------------------------------------------------
# 4. Output
# ---------------------------------------------------------------------------

# Write LAMMPS data file (with force-field coefficients)
builder.write_lammps("data.lammps", atom_style="full", write_coeff=True)

# Convert to ASE Atoms (e.g. for visualisation or further manipulation)
atoms = builder.to_ase()

# Access the raw MDAnalysis Universe if needed
universe = builder.universe
