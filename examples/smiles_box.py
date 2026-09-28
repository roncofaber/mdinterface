"""Build an unparameterized ethanol box from SMILES for ASE-based workflows."""

from mdinterface import SimCell, Specie

ethanol = Specie(smiles="CCO", name="ETOH", seed=42)
box = SimCell(xysize=[25, 25])
box.add_solvent(ethanol, nsolvent=30, zdim=25)
box.build()
box.to_ase().write("ethanol.xyz")
print(ethanol.to_smiles())
