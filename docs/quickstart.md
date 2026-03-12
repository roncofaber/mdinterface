# Quick Start

## Simple workflow

```python
from mdinterface import SimCell
from mdinterface.database import Water, Metal111

# use database structure for convenient pre-built geometry
water = Water()          
gold  = Metal111("Au")

# define builder object
simbox = SimCell(xysize=[15, 15])

# stack your layers
simbox.add_slab(gold, nlayers=3)
simbox.add_solvent(water, zdim=20, density=1.0)
simbox.build()

#convert to ase.Atoms ready for AIMD, ML-MD, or any other tool
atoms = simbox.to_ase()    
```

## LAMMPS workflow

To output a LAMMPS data file, force field information is required. The examples below use the built-in database with pre-parameterised force-field data for writing LAMMPS input files.

## Defining species

Species are the building blocks. They can come from the built-in database or be defined from scratch:

```python
from mdinterface import SimCell
from mdinterface.database import Water, Ion, Metal111

water = Water(model="ewald")          # SPC/E water with Ewald-compatible charges
na    = Ion("Na", ffield="Cheatham")
cl    = Ion("Cl", ffield="Cheatham")
gold  = Metal111("Au")                # FCC Au (111) surface
```

## Electrode / electrolyte sandwich

```python
simbox = SimCell(xysize=[15, 15], verbose=True)

simbox.add_slab(gold, nlayers=3)
simbox.add_solvent(water, solute=[na, cl], nsolute=[5, 5], zdim=25, density=1.0)
simbox.add_slab(gold, nlayers=3)
simbox.add_vacuum(zdim=5)

simbox.build(padding=0.5)
simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True)
```

## Mixed-solvent electrolyte

```python
from mdinterface.core.specie import Specie

methanol = Specie("CH3OH", ligpargen=True)

simbox = SimCell(xysize=[25, 25])
simbox.add_solvent(
    [water, methanol],
    ratio=[3, 1],       # 3 water : 1 methanol by mole
    density=0.95,
    zdim=30,
    solute=[na, cl],
    nsolute=[5, 5],
)
simbox.build(padding=0.5)
simbox.write_lammps("data_mixture.lammps", atom_style="full", write_coeff=True)
```

## Multi-layer system

```python
simbox = SimCell(xysize=[20, 20], verbose=True)

simbox.add_slab(gold, nlayers=4)
simbox.add_vacuum(zdim=3)
simbox.add_solvent(water, zdim=20, density=1.0)
simbox.add_vacuum(zdim=3)
simbox.add_slab(gold, nlayers=4)

simbox.build(padding=0.5, center=True)
```

## Converting the output

```python
atoms    = simbox.to_ase()    # ase.Atoms with cell and PBC
universe = simbox.universe    # mda.Universe
```

## More examples

Complete runnable scripts are in the [`examples/`](https://github.com/roncofaber/mdinterface/tree/main/examples) directory:

| Script | What it shows |
|--------|--------------|
| `electrode_interface.py` | Au / NaCl electrolyte / Au sandwich |
| `solvent_box.py` | Pure solvent + dissolved ions |
| `multisolvent_box.py` | Mixed-solvent box with ratio/density/count modes |
| `multilayer.py` | Five-layer multi-slab system |
| `sandwich_from_traj.py` | Load a relaxed structure via `hijack` |
