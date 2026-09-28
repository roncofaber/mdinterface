# Species

A `Specie` is the fundamental unit in `mdinterface`. It holds the geometry of a single molecule or repeating unit, and optionally its force-field parameters (charges, bonds, angles, …) for LAMMPS output.

## Loading from the database

The easiest way to get a `Specie` is from the built-in database:

```python
from mdinterface.database import Water, Ion, Metal111, Graphene, NobleGas

water   = Water(model="ewald")           # modified TIP3P for Ewald electrostatics
na      = Ion("Na", ffield="Cheatham")
cl      = Ion("Cl", ffield="Cheatham")
gold    = Metal111("Au")                 # FCC (111) gold surface
grap    = Graphene()
argon   = NobleGas("Ar")
```

See the [Database guide](database.md) for a full list of available entries.

## SMILES and RDKit molecules

RDKit is a core dependency. Use the explicit `smiles=` keyword to preserve chemical connectivity, bond orders, formal charges, and stereochemistry:

```python
from mdinterface import Specie
from rdkit import Chem

ethanol = Specie(smiles="CCO", name="ETOH", seed=42)
cation = Specie(smiles="C[N+](C)(C)C", name="TMA")
from_rdkit = Specie(Chem.MolFromSmiles("C[C@H](O)C(=O)O"), seed=42)

mol = cation.to_rdkit()
smiles = ethanol.to_smiles()
mapped_smiles = ethanol.to_smiles(mapped=True)
```

Explicit mapped hydrogens in SMILES are preserved. RDKit adds any remaining hydrogens and generates a conformer with ETKDG when coordinates are absent. Existing RDKit coordinates are retained. Inputs must describe one connected species; create separate species for disconnected salts. `Specie("CO")` retains its ASE molecule-name meaning, while `Specie(smiles="CO")` creates methanol. `smiles` cannot be combined with `atoms` or `lammps_data`.

The chemical graph stays fixed when coordinates change. `Specie.graph` remains a NetworkX graph and includes bond orders for RDKit-backed structures. `to_rdkit()` returns an independent molecule in ASE atom order with current coordinates. For ASE inputs without a stored chemical graph, it attempts RDKit bond perception using the specified total charge and raises if the chemistry cannot be determined. Periodic slabs retain the existing geometric graph and do not require a molecular RDKit representation.

`nominal_charge` contains integer formal-charge annotations; `atom_map` contains unique persistent atom identifiers. These are separate from force-field partial charges. A SMILES input alone does not provide OPLS parameters or partial charges. Supply `ligpargen=True` to parameterize, or provide force-field data explicitly. Conflicting `tot_charge` or modified formal-charge annotations raise an error rather than silently changing the chemistry.

RDKit-backed inputs are passed to LigParGen as MOL files, preserving bonds and formal charges. This path requires the mdinterface-compatible LigParGen fork at commit [`ad78036`](https://github.com/roncofaber/ligpargen/commit/ad78036842318f166531be41cfcbc3563d7c5476) pinned in the [installation instructions](../installation.md); older revisions reconstruct the molecule from PDB coordinates and can lose chemical information. ASE-only inputs retain the XYZ/Open Babel path.

`repeat()` scales the stored total charge with the number of copies and retains the first copy's atom maps while assigning distinct maps to additional copies.

## Preparing molecular geometry

```python
ethanol.generate_conformer(seed=42)
energy = ethanol.minimize_geometry(max_iterations=1000)
# Or embed and minimize together, applying coordinates only if both succeed:
energy = ethanol.generate_conformer(seed=42, minimize=True)
```

`generate_conformer()` embeds a new ETKDGv3 conformer; `minimize_geometry()` applies MMFF94 to existing coordinates. Both preserve chemistry, atom maps and simulation parameters and leave coordinates unchanged on failure. Periodic, disconnected or constrained structures are rejected. MMFF94 energies are in kcal/mol and are not OPLS-AA simulation energies. See the [polymer workflow](polymer.md) for junction parameterization and charge auditing.

## Charge estimation

`estimate_charges(method="ligpargen")`, `estimate_charges(method="obabel")`, `estimate_charges(method="resp")`, and `estimate_OPLSAA_parameters()` use the species' stored total charge. An explicit conflicting `charge` raises an error. Set `tot_charge` when constructing a coordinate-only species; SMILES-backed species derive it from the chemical graph. Open Babel receives RDKit-backed chemistry through MOL data, and RESP receives the total charge explicitly when constructing its electronic structure calculation.

## Defining a Specie manually

For molecules not in the database, create a `Specie` from an ASE `Atoms` object and set the parameters manually:

```python
from ase import Atoms
from mdinterface import Specie

mol = Atoms("OCO", positions=[[0,0,0],[1.16,0,0],[-1.16,0,0]])
co2 = Specie(mol, charges=[-0.3298, 0.6596, -0.3298])
# ... set atom types, bonds, angles as needed
```

## Generating OPLS-AA parameters with LigParGen

For organic molecules, `mdinterface` can call LigParGen to obtain OPLS-AA force-field parameters automatically:

```python
from mdinterface import Specie

methanol = Specie("CH3OH", ligpargen=True)
```

This requires a working LigParGen installation and the `BOSSdir` configured in `config.ini` (see [Installation](../installation.md)).

Missing LigParGen, Open Babel, or BOSS configuration and failed LigParGen runs raise `LigParGenError` with actionable setup guidance. Failed runs retain their temporary directory and `ligpargen.log` for inspection.

For an existing species, `parameterize()` assigns a complete LigParGen result atomically:

```python
cation = Specie(smiles="C[N+](C)(C)C")
report = cation.parameterize(charge_correction="none")
cation.validate_force_field()
```

Coordinates stay unchanged. The charge audit reports the initial, raw refined and final totals, formal target, residual, correction per atom and number of junctions. Charge correction is always explicit: `"none"` is the default, and `"uniform"` distributes the residual over all atoms. `ligpargen=True` uses this same parameterization path. Charge-audit totals and residuals within `1e-12 e` of zero are reported as `0.0`, without applying a numerical-noise correction or changing individual atomic partial charges.

**Large molecules (>200 atoms):** `parameterize()` splits the molecule into chemically valid capped segments and refines the junctions. `segment_size` counts caps and cannot exceed 200. All segments and expanded junction snippets are checked before any LigParGen calculation. If the molecule cannot be split safely or a junction remains too large, parameterization fails with guidance rather than sending an oversized structure to BOSS. `snippet_radius` controls junction extent, and `cap_element` selects a neutral monovalent cap. All calculations and validation must succeed before parameters are applied. `refine_large_topology()` delegates to the same workflow with the same options and report.

`validate_force_field()` checks finite charges, assigned pair coefficients, and molecular bond, angle and proper-torsion coverage. Explicit zero-valued coefficients are valid. Assigned improper coefficients are checked, but missing impropers cannot be inferred universally from connectivity. This validation does not select constraints or establish compatibility between different force-field conventions.

## OpenFF status

OpenFF-to-LAMMPS export and re-import into `Specie` have been validated, but OpenFF is not yet a supported parameterization backend. The LAMMPS data file preserves numeric coefficients but does not encode every required convention, including Fourier torsion style, mixing rules, switching behavior, and 1-4 scaling. Until mdinterface carries this metadata through system assembly and output, do not treat OpenFF-generated coefficients as OPLS coefficients or mix OpenFF and LigParGen molecular species without an independent force-field compatibility analysis.

## Polymers

`Polymer` extends `Specie` to build linear chains from one or more monomer units, including co-polymers with arbitrary sequences. See the dedicated [Polymer guide](polymer.md) for the full workflow.

## Inspecting a Specie

```python
specie.atoms          # ase.Atoms
specie.universe       # mda.Universe
print(specie)         # summary: formula, n atoms, charges, ...
```
