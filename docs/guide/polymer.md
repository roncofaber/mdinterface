# Polymer Builder

`Polymer` is a subclass of `Specie` that assembles a linear chain from one or more monomer `Specie` objects. Assembly preserves monomer chemistry and records explicit junction bonds. Parameterized `Specie` monomers carry their parameters into the chain automatically, with separate atom labels for every repeat. Refine junction parameters before classical MD export. ASE-only monomers still require explicit bonded parameter lists.

```python
from mdinterface import Polymer
```

## Preparing a monomer

Each monomer must have two arrays set on its ASE `Atoms` object before polymerization:

| Array | Values | Role |
|-------|--------|------|
| `polymerize` | `0` = not a junction, `1` = head, `2` = tail | Marks the leaving atom (e.g. H or F) at each chain end that will be **removed** to form the inter-monomer bond |
| `nominal_charge` | integer per atom | Tracks the formal charge contribution of each atom; used to keep the total chain charge consistent after topology refinement |

```python
import numpy as np

# load monomer from a LAMMPS data file (or build from scratch)
from mdinterface.io.read import read_lammps_data_file
mon, atoms, bonds, angles, dihedrals, impropers = read_lammps_data_file("monomer.lammps.lmp")

# mark the leaving atoms at the two chain ends (e.g. H or F)
mon.new_array("polymerize", np.zeros(len(mon), dtype=int))
mon.arrays["polymerize"][head_idx] = 1   # leaving atom at the head end
mon.arrays["polymerize"][tail_idx] = 2   # leaving atom at the tail end

# set formal charges (0 for neutral atoms, integer for charged sites)
mon.set_array("nominal_charge", np.zeros(len(mon), dtype=int))
mon.arrays["nominal_charge"][charged_atom_idx] = 1   # e.g. +1 on N
```

> The `polymerize`-marked atoms are the **leaving atoms** at each chain end (e.g. H or F depending on the monomer chemistry). They are deleted during assembly, and the bond is formed between the heavy atoms they were attached to.

Monomers may be ASE `Atoms` or `Specie` objects. For SMILES-backed species, formal charges and bond orders come directly from RDKit; mark the leaving hydrogens in `monomer.atoms.arrays["polymerize"]`. Atom maps remain unique in the assembled chain. ASE-only monomers undergo RDKit bond perception, and inconsistent nonzero formal-charge annotations are rejected. A monomer must be one connected molecule with exactly one head and one tail, each marking a neutral, singly bonded leaving atom.

### Attachment sites from mapped SMILES

Use atom maps to identify attachment sites without depending on input atom order:

```python
from mdinterface import Specie, Polymer

monomer = Specie(smiles="[CH3:1][N+](C)(C)C[CH3:2]", name="MON")
monomer.parameterize()
monomer.mark_attachment_sites(head_map=1, tail_map=2)
chain = Polymer(monomer, nrep=3, name="POL")
```

`mark_attachment_sites()` accepts mapped leaving atoms or mapped anchors with terminal leaving atoms of `leaving_element` (default `"H"`). Missing maps, charged or nonterminal leaving atoms, duplicate sites, and ambiguous non-hydrogen choices raise an error without replacing existing annotations. Multiple equivalent hydrogens on a center without a specified chiral tag are accepted. Explicitly map the desired leaving atom when stereochemistry matters.

`Polymer` automatically inherits charges, atom types, and bonded parameters from `Specie` monomers when no explicit topology overrides are supplied. Each repeat receives its own labels, so identical labels in different monomers do not collide and local junction updates remain local. Source monomers are not modified. If you supply explicit topology arguments, the legacy explicit-data path is used instead of automatic bonded-parameter inheritance.

The runnable `examples/polymer/polymer_from_smiles.py` continues through geometry preparation, junction refinement, packing with chloride counterions, and validated LAMMPS export. It requires LigParGen and BOSS.

## Building the chain

### Homopolymer

Pass a parameterized `Specie` monomer and `nrep` to repeat it:

```python
from mdinterface import Polymer

chain = Polymer(monomer, nrep=20, name="POLY")
```

### Co-polymer with an explicit sequence

Pass a list of parameterized `Specie` monomers and a `sequence` of monomer indices:

```python
import random

seq = 17 * [0] + 3 * [1]   # 17 × mon_A and 3 × mon_B
random.shuffle(seq)

chain = Polymer(
    monomers  = [mon_A, mon_B],
    sequence  = seq,
    name      = "COPOL",
)
```

## Recommended preparation order

Assemble the chemical graph, generate and minimize the initial geometry, refine junction parameters, inspect the charge audit, and only then round partial charges if needed. Pack the prepared chains with solvent and counterions, validate the exported system, and equilibrate under the intended simulation conditions.

### Whole-chain geometry with RDKit

After assembling the chemical graph, generate the whole-chain coordinates with RDKit before junction parameterization:

```python
chain.generate_conformer(seed=42)
energy = chain.minimize_geometry(max_iterations=5000)
print(f"MMFF94 energy: {energy:.6f} kcal/mol")
```

`generate_conformer()` replaces the assembled coordinates using [RDKit ETKDGv3](https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html) with optional MMFF94 minimization via `minimize=True`. Call `minimize_geometry()` to relax existing coordinates without embedding another conformer. Both methods are defined on `Specie` and inherited by `Polymer`. They preserve atom order, junction bonds, formal charges, partial charges, and simulation topology. A fixed seed makes generation reproducible within the same RDKit version. Missing MMFF parameters, embedding failure, or unconverged minimization raises an error without changing the original coordinates. Periodic structures and ASE constraints are not supported by this operation.

MMFF94 is used only to prepare the geometry; it does not supply the OPLS-AA parameters exported to LAMMPS. This is one minimized isolated-chain conformer, not an equilibrated polymer ensemble. Long chains may need another seed, more minimization iterations, or a different preparation method. The existing assembled coordinates remain available if you skip this optional step.

ML relaxation is optional. If used, choose a model supporting the chain's chemistry and charge state and supply its actual formal charge, never a hard-coded zero for an ionomer. Production equilibration should use the chosen simulation force field and target thermodynamic conditions.

### Topology refinement with LigParGen

Junction atoms sit in a chemical environment that neither monomer alone can describe correctly. `refine_junctions` calls LigParGen on small snippets around each junction to obtain accurate OPLS-AA atom types and partial charges:

```python
report = chain.refine_junctions(
    snippet_radius=12,      # graph distance in bonds around the junction
    charge_correction="uniform",  # explicitly correct the chain total
    cap_element="H",   # element used to cap dangling bonds in each snippet
)
print(report)
```

> Requires a working LigParGen installation (see [Installation](../installation.md)). Only needed for classical MD with LAMMPS.

Refinement stages topology and charge changes and applies them only after every junction succeeds, leaving the original chain unchanged on failure. All capped snippets are checked against LigParGen's 200-atom limit before any backend calculation; preserved rings and charged groups can expand a snippet beyond its requested radius. `chain.junction_bonds` exposes inter-monomer bonds as pairs of current ASE atom indices.

Refinement requires a connected graph and `snippet_radius >= 4`. Snippets preserve touched rings, formal-charge sites and their immediate neighbors, and multiple bonds; artificial caps are neutral. Charges and topology are recalculated for every junction without the former element-string cache. RDKit-backed snippets use the MOL handoff described in the [species guide](specie.md).

Assembly removes leaving-atom partial charges, so an unrefined chain need not have its intended total charge. `round_charges()` only rounds the existing total. Refinement returns a charge audit containing the initial total, formal target, refined total, residual (refined minus formal), per-atom correction, final total, and junction count. It also logs the residual before applying the optional uniform `charge_correction="uniform"` correction; this correction enforces the chain total, not an integer charge on each repeat unit. Check local charge changes as well as the total. The piperion example runs refinement and checks both the chain charge and final membrane neutrality before export. Its whole-chain RDKit geometry step replaces the former optional ML-relaxation recipe.

### Serialization (optional)

Polymer objects can take a while to prepare. Saving them to disk lets you reuse a pre-built chain without repeating the refinement steps. `dill` is recommended over `pickle` because polymer objects can contain non-picklable internals:

```python
import dill

# save
with open("chain.pkl", "wb") as f:
    dill.dump(chain, f)

# reload in a later session
with open("chain.pkl", "rb") as f:
    chain = dill.load(f)
```

## Using a polymer in SimCell

A `Polymer` is a `Specie`, so it slots directly into the normal `SimCell` workflow. Pass pre-built chains as solutes alongside solvent(s) and other solute(s).

!!! tip Setting solvent content explicitly
    For polymer membranes, **solvent density is not a meaningful input**: the starting simulation box is likely much larger than the equilibrated cell, and NPT MD will shrink it to the correct density. The solvent content should instead be set explicitly. For example, the hydration number λ (number of water molecules per ionic site in ionomers) and converted to an explicit molecule count.

```python
import dill
from mdinterface import SimCell
from mdinterface.database import Water, Ion

with open("chain.pkl", "rb") as f:
    chain = dill.load(f)

water = Water()
cl    = Ion("Cl", ffield="opls-aa")

n_chains        = 15
n_sites         = 17          # ionic sites per chain (e.g. 17 ammonium groups)
lam             = 20          # lambda: water molecules per ionic site
n_cl            = n_chains * n_sites
n_water         = n_chains * n_sites * lam

# use a large initial box -- NPT MD will equilibrate it to the correct density
simbox = SimCell(xysize=[400, 400])
simbox.add_solvent(
    water,
    zdim     = 400,
    nsolvent = n_water,
    solute   = [chain, cl],
    nsolute  = [n_chains, n_cl],
)
simbox.build(padding=0.5)
simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True)
```
