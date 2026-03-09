#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BoxBuilder example: Au(111) / NaCl electrolyte / Au(111) sandwich.

Equivalent to make_box.py but using the fluent BoxBuilder API.
Layers are stacked in the order the methods are called, so the
geometry is immediately readable from the code.

Author: roncofaber
"""

from mdinterface import BoxBuilder
from mdinterface.database import Water, Ion, Metal111

# ---------------------------------------------------------------------------
# 1. Define species
# ---------------------------------------------------------------------------

water = Water(model="ewald")

na = Ion("Na", ffield="Cheatham")
cl = Ion("Cl", ffield="Cheatham")

gold = Metal111("Au")

# ---------------------------------------------------------------------------
# 2. Assemble layers  (bottom → top)
# ---------------------------------------------------------------------------

builder = (
    BoxBuilder(xysize=[15, 15])
        .add_slab(gold, nlayers=1)
        .add_solvent(
            water,
            ions=[na, cl],
            nions=[5, 5],  # 5 Na⁺ and 5 Cl⁻
            zdim=25,
            density=1.0,
        )
        .add_slab(gold, nlayers=1)
        .add_vacuum(zdim=5)
)

# ---------------------------------------------------------------------------
# 3. Build
# ---------------------------------------------------------------------------

builder.build(
    padding=0.5,
    center=False,   # set True to center the first slab in the box
    layered=False,  # set True to tag each slab layer with a unique mol-id
)

# ---------------------------------------------------------------------------
# 4. Output
# ---------------------------------------------------------------------------

# Write LAMMPS data file (with force-field coefficients)
builder.write_lammps("data.lammps", atom_style="full", write_coeff=True)

# Convert to ASE Atoms (e.g. for visualisation or further manipulation)
atoms = builder.to_ase()

# Access the raw MDAnalysis Universe if needed
universe = builder.universe
