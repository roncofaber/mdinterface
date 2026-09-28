"""Build a charged polymer from mapped SMILES with LigParGen parameters."""

import numpy as np

from mdinterface import Polymer, Specie, SimCell
from mdinterface.database import Ion

monomer = Specie(smiles="[CH3:1][N+](C)(C)C[CH3:2]", name="MON")
monomer.parameterize()
monomer.mark_attachment_sites(head_map=1, tail_map=2)
chain = Polymer(monomer, nrep=3, name="POL")
chain.generate_conformer(seed=42)
energy = chain.minimize_geometry(max_iterations=5000)
report = chain.refine_junctions(snippet_radius=6, charge_correction="uniform")
assert len(chain.junction_bonds) == 2
assert np.isclose(report["final_charge"], 3.0)
print(f"MMFF94 geometry energy: {energy:.6f} kcal/mol")
print("Junction bonds:", chain.junction_bonds)
print("Charge audit:", report)
chain.atoms.write("polymer.xyz", format="xyz")
box = SimCell(xysize=[40, 40])
box.add_solvent(chain, nsolvent=1, zdim=40, solute=[Ion("Cl", ffield="opls-aa")], nsolute=[3])
box.build()
assert np.isclose(box.universe.atoms.charges.sum(), 0.0)
box.write_lammps("polymer.data", expected_charge=0.0)
