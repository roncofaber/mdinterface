# Species

A `Specie` is the fundamental unit in `mdinterface`. It holds the geometry of a single molecule or repeating unit, and optionally its force-field parameters (charges, bonds, angles, …) for LAMMPS output.

## Loading from the database

The easiest way to get a `Specie` is from the built-in database:

```python
from mdinterface.database import Water, Ion, Metal111, Graphene, NobleGas

water   = Water(model="ewald")           # SPC/E water
na      = Ion("Na", ffield="Cheatham")
cl      = Ion("Cl", ffield="Cheatham")
gold    = Metal111("Au")                 # FCC (111) gold surface
grap    = Graphene()
argon   = NobleGas("Ar")
```

See the [Database guide](database.md) for a full list of available entries.

## Defining a Specie manually

For molecules not in the database, create a `Specie` from an ASE `Atoms` object and set the parameters manually:

```python
from ase import Atoms
from mdinterface import Specie

mol = Atoms("CO2", positions=[[0,0,0],[1.16,0,0],[-1.16,0,0]])
co2 = Specie(mol)
co2.set_charges([-0.3298, 0.6596, -0.3298])
# ... set atom types, bonds, angles as needed
```

## Generating OPLS-AA parameters with LigParGen

For organic molecules, `mdinterface` can call LigParGen to obtain OPLS-AA force-field parameters automatically:

```python
from mdinterface import Specie

methanol = Specie("CH3OH", ligpargen=True)
```

This requires a working LigParGen installation and the `BOSSdir` configured in `config.ini` (see [Installation](../installation.md)).

**Large molecules (>200 atoms):** LigParGen's input limit is 200 atoms. When `ligpargen=True` and the molecule exceeds this limit, `mdinterface` automatically splits it into segments along clean backbone bonds, runs LigParGen independently on each segment, and then refines the parameters at every junction using a local snippet — all transparently. No extra configuration is needed.

## Polymers

`Polymer` extends `Specie` to build linear chains from one or more monomer units — including co-polymers with arbitrary sequences. See the dedicated [Polymer guide](polymer.md) for the full workflow.

## Inspecting a Specie

```python
specie.atoms          # ase.Atoms
specie.universe       # mda.Universe
print(specie)         # summary: formula, n atoms, charges, ...
```
