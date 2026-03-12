# Polymer Builder

`Polymer` is a subclass of `Specie` that assembles a linear chain from one or more monomer `Specie` objects. All the geometry, bonding, and force-field data of the individual monomers are carried through to the assembled chain.

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

## Building the chain

### Homopolymer

Pass a single monomer and `nrep` to repeat it:

```python
from mdinterface import Polymer

chain = Polymer(monomer, nrep=20, bonds=bonds, angles=angles,
                dihedrals=dihedrals, impropers=impropers, name="POLY")
```

### Co-polymer with an explicit sequence

Pass a list of monomers and a `sequence` of monomer indices:

```python
import random

seq = 17 * [0] + 3 * [1]   # 17 × mon_A and 3 × mon_B
random.shuffle(seq)

chain = Polymer(
    monomers  = [mon_A, mon_B],
    sequence  = seq,
    bonds     = bonds_A + bonds_B,
    angles    = angles_A + angles_B,
    dihedrals = dihedrals_A + dihedrals_B,
    impropers = impropers_A + impropers_B,
    name      = "COPOL",
)
chain.round_charges(7)
```

## Best practices

### Geometry refinement

The geometry at inter-monomer junctions is approximate after assembly. A **quick geometry relaxation** is normally needed to fix obviously bad bond lengths.

#### Geometry relaxation with FAIRChem

Hookean constraints at the new bonds prevent the chain from collapsing during the first optimization steps:

```python
import ase.constraints
from ase.data import atomic_numbers, covalent_radii
from fairchem.core import FAIRChemCalculator

# add Hookean springs at each junction
constraints = []
for pair in chain._get_connection_elements():
    d1 = covalent_radii[atomic_numbers[chain.atoms[pair[0]].symbol]]
    d2 = covalent_radii[atomic_numbers[chain.atoms[pair[1]].symbol]]
    constraints.append(
        ase.constraints.Hookean(int(pair[0]), int(pair[1]), 8, rt=d1 + d2)
    )
chain.atoms.set_constraint(constraints)

# attach a machine-learning potential
calc = FAIRChemCalculator.from_model_checkpoint(
    "uma-s-1.pt", device="cuda", task_name="omol"
)
chain.atoms.calc = calc
chain.atoms.info["charge"] = 0
chain.atoms.info["spin"]   = 0

chain.relax_structure(trajectory="relax.traj", optimizer="FIRE", steps=1000)

# remove constraints before topology refinement
chain.atoms.set_constraint()
```

### Topology refinement with LigParGen

Junction atoms sit in a chemical environment that neither monomer alone can describe correctly. `refine_polymer_topology` calls LigParGen on small snippets around each junction to obtain accurate OPLS-AA atom types and partial charges:

```python
chain.refine_polymer_topology(
    Nmax=12,      # atoms on each side of the junction to include in the snippet
    offset=True,  # shift charges so the total matches nominal_charge.sum()
    ending="H",   # element used to cap dangling bonds in each snippet
)
```

> Requires a working LigParGen installation (see [Installation](../installation.md)). Only needed for classical MD with LAMMPS.

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
    For polymer membranes, **solvent density is not a meaningful input**: the starting simulation box is likely much larger than the equilibrated cell, and NPT MD will shrink it to the correct density. The solvent content should instead be setexplicitly. For example, the hydration number λ (number of water molecules per ionic site in ionomers) and converted to an explicit molecule count.

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
