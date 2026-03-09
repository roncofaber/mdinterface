#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BoxBuilder example: multi-layer cell with heterogeneous electrolyte regions.

This setup is impractical to express with the old SimulationBox API (which
only supports one interface, one enderface, and one miderface slot). With
BoxBuilder you can stack as many slabs and solvent regions as needed.

Layout (bottom to top):
    Au(111) slab        — bottom electrode
    NaCl solution       — first electrolyte region
    Pt(111) slab        — middle electrode
    KF solution         — second electrolyte region
    Pt(111) slab        — top electrode
    vacuum gap

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
k  = Ion("K",  ffield="Cheatham")
f  = Ion("F",  ffield="Dang")

gold     = Metal111("Au")
platinum = Metal111("Pt")

# ---------------------------------------------------------------------------
# 2. Assemble layers  (bottom → top)
# ---------------------------------------------------------------------------

builder = (
    BoxBuilder(xysize=[20, 20])
        .add_slab(gold,     nlayers=3)                                  # bottom electrode
        .add_solvent(water, ions=[na, cl], nions=[4, 4], zdim=20, density=1.0)  # NaCl region
        .add_slab(platinum, nlayers=2)                                  # middle electrode
        .add_solvent(water, ions=[k,  f],  nions=[4, 4], zdim=20, density=1.0)  # KF region
        .add_slab(platinum, nlayers=2)                                  # top electrode
        .add_vacuum(zdim=10)                                            # gap at periodic boundary
)

# ---------------------------------------------------------------------------
# 3. Build
# ---------------------------------------------------------------------------

builder.build(padding=0.5)

# ---------------------------------------------------------------------------
# 4. Output
# ---------------------------------------------------------------------------

# Write LAMMPS data file (with force-field coefficients)
builder.write_lammps("data_multilayer.lammps", atom_style="full", write_coeff=True)

# Convert to ASE Atoms (e.g. for visualisation or further manipulation)
atoms = builder.to_ase()

# Access the raw MDAnalysis Universe if needed
universe = builder.universe
