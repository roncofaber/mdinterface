#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Polymer example: piperion anion-exchange membrane.

Two monomers are loaded from LAMMPS data files, marked with polymerization
connection points, and assembled into a co-polymer chain.  The chain is then
packed with water and chloride ions into a membrane box, using an explicit
water count (hydration number lambda) rather than a density target.

RDKit generates and minimizes the whole-chain starting geometry. Production
work still requires MD equilibration appropriate to the target system.

Author: roncofaber
"""

from pathlib import Path

import numpy as np

from mdinterface import Polymer, SimCell
from mdinterface.database import Water, Ion
from mdinterface.io.read import read_lammps_data_file

DATA = Path(__file__).parent / "data"

#%% Load monomers from LAMMPS data files

mon1, atoms1, bonds1, angles1, dihedrals1, impropers1 = \
    read_lammps_data_file(DATA / "mon1" / "monomer_1.lammps.lmp")

# ato_start_idx offsets the atom-type indices of mon2 so they do not clash
# with those of mon1 when the two topology lists are combined below.
mon2, atoms2, bonds2, angles2, dihedrals2, impropers2 = \
    read_lammps_data_file(DATA / "mon2" / "monomer_2.lammps.lmp",
                          ato_start_idx=len(atoms1))

#%% Mark polymerization connection points and formal charges

# polymerize array: 0 = not a junction, 1 = head, 2 = tail.
# The marked atoms are the leaving atoms (H in this case) at each chain end.
# They are deleted during assembly; the bond forms between the heavy atoms they
# were attached to.

mon1.new_array("polymerize", np.zeros(len(mon1), dtype=int))
mon1.arrays["polymerize"][26] = 1   # head H atom
mon1.arrays["polymerize"][47] = 2   # tail H atom

mon2.new_array("polymerize", np.zeros(len(mon2), dtype=int))
mon2.arrays["polymerize"][32] = 1   # head H atom
mon2.arrays["polymerize"][39] = 2   # tail H atom

# nominal_charge: formal integer charge per atom.
# mon1 carries one piperidinium cation; mon2 is neutral.
mon1.set_array("nominal_charge", np.zeros(len(mon1), dtype=int))
nitrogens = [atom.index for atom in mon1 if atom.symbol == "N"]
if len(nitrogens) != 1:
    raise ValueError("Expected one piperidinium nitrogen in monomer 1.")
mon1.arrays["nominal_charge"][nitrogens[0]] = 1

mon2.set_array("nominal_charge", np.zeros(len(mon2), dtype=int))

#%% Build the co-polymer chain

# Sequence of monomer indices: 0 -> mon1 (charged), 1 -> mon2 (neutral).
# Here: 3 charged units and 1 neutral linker for a short demonstration.
# Typical production chains use 17-20 units (e.g. 17×mon1 + 3×mon2).
sequence = [0, 1, 0, 0]   # 3 × mon1, 1 × mon2

chain = Polymer(
    monomers  = [mon1, mon2],
    sequence  = sequence,
    bonds     = bonds1 + bonds2,
    angles    = angles1 + angles2,
    dihedrals = dihedrals1 + dihedrals2,
    impropers = impropers1 + impropers2,
    name      = "PI00",
)

#%% Whole-chain geometry from the molecular graph

energy = chain.generate_conformer(minimize=True, seed=42, max_iterations=5000)
print(f"Initial geometry: MMFF94 energy {energy:.6f} kcal/mol")

#%% Topology refinement at junctions (requires LigParGen)
#
# Refines OPLS-AA atom types and partial charges at each inter-monomer
# junction using LigParGen on small local snippets.
#
charge_report = chain.refine_junctions(snippet_radius=12, charge_correction="uniform", cap_element="H")
print("Junction charge audit:", charge_report)
chain.round_charges(7)
if not np.isclose(chain.charges.sum(), sequence.count(0), atol=1e-6):
    raise ValueError("Refined polymer charge does not match its ionic sites.")

#%% MD equilibration
#
# A minimized isolated-chain conformer is not an equilibrated membrane.
# Equilibrate with the chosen simulation force field and thermodynamic conditions.
# Once equilibrated, save the chain with dill for easy reuse:
#
# import dill
# with open("chain.pkl", "wb") as f:
#     dill.dump(chain, f)

#%% Build the membrane box

water = Water()
cl    = Ion("Cl", ffield="opls-aa")

# Water content is set by the hydration number lambda (water molecules per
# ionic site) rather than a density.  The box starts large and NPT MD shrinks
# it to the correct density during equilibration.
n_chains = 5
n_mon1   = sequence.count(0)   # charged monomers per chain
lam      = 10                  # water molecules per ionic site
n_cl     = n_chains * n_mon1
n_water  = n_chains * n_mon1 * lam

simbox = SimCell(xysize=[200, 200], verbose=True)

simbox.add_solvent(
    water,
    zdim     = 200,
    nsolvent = n_water,
    solute   = [chain, cl],
    nsolute  = [n_chains, n_cl],
)

simbox.build(padding=0.5)
if not np.isclose(simbox.universe.atoms.charges.sum(), 0.0, atol=1e-6):
    raise ValueError("Polymer membrane is not charge neutral.")

#%% Output

simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True, expected_charge=0.0)

atoms    = simbox.to_ase()
universe = simbox.universe
