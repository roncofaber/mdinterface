#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BoxBuilder example: Au(111) / NaCl electrolyte / Au(111) sandwich.

Equivalent to make_box.py but using the fluent BoxBuilder API.
Layers are added in the order they appear in the cell (bottom → top).

Author: roncofaber
"""

from mdinterface import BoxBuilder
from mdinterface.database import Water, Ion, Metal111

#%% Define species

water = Water(model="ewald")

na = Ion("Na", ffield="Cheatham")
cl = Ion("Cl", ffield="Cheatham")

gold = Metal111("Au")

#%% Set up simulation box

simbox = BoxBuilder(xysize=[15, 15])

#%% Add layers  (bottom → top)

simbox.add_slab(gold, nlayers=1)

simbox.add_solvent(
    water,
    ions=[na, cl],
    nions=[5, 5],           # 5 Na⁺ and 5 Cl⁻
    zdim=25,
    density=1.0,
    # dilate=1.25,          # for concentrated systems: packs into a box 25%
    #                       # taller so PACKMOL has breathing room; NpT MD
    #                       # will compress back to the correct density
    # packmol_tolerance=1.5,# lower from default 2.0 Å if PACKMOL fails to pack
)

simbox.add_slab(gold, nlayers=1)

simbox.add_vacuum(zdim=5)

#%% Build

simbox.build(
    padding=0.5,
    center=False,   # set True to center the first slab in the box
    layered=False,  # set True to tag each slab layer with a unique mol-id
    # hijack=atoms, # override positions + cell from an external ase.Atoms
    #               # (e.g. a pre-relaxed slab); topology is preserved
)

#%% Output

# Write LAMMPS data file (with force-field coefficients)
simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True)

# Convert to ASE Atoms (e.g. for visualisation or further manipulation)
atoms = simbox.to_ase()

# Access the raw MDAnalysis Universe if needed
universe = simbox.universe
