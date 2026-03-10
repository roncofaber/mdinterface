#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BoxBuilder example: mixed-solvent box (water + methanol).

Demonstrates the three ways to specify solvent composition when multiple
solvent species are present.  A dissolved NaCl salt is added on top.

Mixing modes (choose one):

  A) ratio + density   -- most physical: volumes/masses set the count
  B) ratio + nsolvent  -- fix total molecule count, split by ratio
  C) nsolvent (list)   -- full manual control, one count per species

Author: roncofaber
"""

from mdinterface import BoxBuilder
from mdinterface.database import Water, Ion
from mdinterface.core.specie import Specie

#%% Define species

# Solvents
water    = Water(model="ewald")
methanol = Specie("CH3OH", ligpargen=True)

# Dissolved ions (optional — comment out if not needed)
na = Ion("Na", ffield="Cheatham")
cl = Ion("Cl", ffield="Cheatham")

#%% Set up simulation box

simbox = BoxBuilder(
    xysize=[25, 25],
    verbose=True,
)

#%% Add mixed-solvent layer

# --- Mode A: ratio + density (recommended for liquid mixtures) ---
# 3 water molecules for every 1 methanol, total density 0.95 g/cm³.
# Molecule counts are derived automatically from the molar masses.
simbox.add_solvent(
    [water, methanol],
    ratio=[3, 1],           # molar mixing ratio (water : methanol)
    density=0.95,           # mixture density in g/cm³
    zdim=30,                # region thickness in Å
    ions=[na, cl],
    nions=[5, 5],           # 5 Na+ and 5 Cl-
)

# --- Mode B: ratio + fixed total count ---
# Uncomment to use instead of Mode A.
# simbox.add_solvent(
#     [water, methanol],
#     ratio=[3, 1],
#     nsolvent=200,          # 150 water + 50 methanol (split proportionally)
#     zdim=30,
#     ions=[na, cl],
#     nions=[5, 5],
# )

# --- Mode C: explicit per-species counts ---
# Uncomment to use instead of Mode A.
# simbox.add_solvent(
#     [water, methanol],
#     nsolvent=[150, 50],    # exactly 150 water and 50 methanol
#     zdim=30,
#     ions=[na, cl],
#     nions=[5, 5],
# )

#%% Build

simbox.build(padding=0.5)

#%% Output

simbox.write_lammps("data_multisolvent.lammps", atom_style="full", write_coeff=True)

atoms    = simbox.to_ase()
universe = simbox.universe
