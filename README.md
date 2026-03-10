<div style="display: flex; align-items: center;">
  <img src="./assets/mdinterface.png" alt="Logo" width="80"  style="margin-right: 10px;">
  <div style="display: flex; flex-direction: column;">
    <h1 style="margin: 0;">mdinterface: Build Interface Systems for Molecular Dynamics Simulations</h1>
  </div>
</div>

[![PyPI version](https://badge.fury.io/py/mdinterface.svg?icon=si%3Apython)](https://pypi.org/project/mdinterface/) [![GitHub version](https://badge.fury.io/gh/roncofaber%2Fmdinterface.svg?icon=si%3Agithub)](https://github.com/roncofaber/mdinterface)

`mdinterface` is a Python package for building systems for Molecular Dynamics (MD) simulations. Initially developed for electrolyte/electrode solid-liquid interfaces, it is equally suited for pure solvent boxes, mixed-solvent electrolytes, and polymer networks.

## Features

- **Fluent `BoxBuilder` API** -- stack slabs, solvent regions, and vacuum gaps layer by layer; call `.build()` when done.
- **Multi-solvent support** -- mix solvents by molar ratio + density, ratio + total count, or explicit per-species molecule counts.
- **Ion placement** -- dissolve ions by count, molar concentration, or a spatially-varying concentration profile.
- **PACKMOL integration** -- handles molecular packing automatically; tolerance and dilation are tunable per layer.
- **Force-field database** -- pre-defined parameters for common metals, noble gases, water models, and ions; or generate OPLS-AA parameters on the fly with [LigParGen](https://github.com/Isra3l/ligpargen).
- **Polymer builder** -- generate chains of arbitrary length from a monomer `Specie`.
- **RESP charges** -- estimate partial charges with [PySCF](https://github.com/pyscf/pyscf) / [gpu4pyscf](https://github.com/pyscf/gpu4pyscf) (optional).
- **AIMD with FAIRChem** -- run ML-potential dynamics via [FAIRChem](https://github.com/facebookresearch/fairchem) (optional).
- **LAMMPS output** -- writes data files and force-field coefficient blocks ready to run.
- **MDAnalysis integration** -- every object converts to `mda.Universe` with a single call.

## Requirements

Check [requirements.txt](requirements.txt) for mandatory dependencies. `pip install mdinterface` handles them automatically.

You also need `packmol` installed and on your `PATH`:

```bash
conda install -c conda-forge packmol
```

### Optional packages

#### LigParGen (automatic OPLS-AA parameters)

Follow the instructions on the [LigParGen GitHub](https://github.com/Isra3l/ligpargen) (or try [this fork](https://github.com/roncofaber/ligpargen) if you hit installation issues). Point `mdinterface` to your BOSS directory via `config.ini`:

```ini
# ~/.config/mdinterface/config.ini  (path is OS-dependent)
[settings]
BOSSdir = /path/to/boss
```

#### RESP charges with PySCF

Install [PySCF](https://github.com/pyscf/pyscf) and [PyMBXAS](https://gitlab.com/roncofaber/pymbxas). RESP fitting currently requires [gpu4pyscf](https://github.com/pyscf/gpu4pyscf).

#### AIMD with FAIRChem

```bash
pip install fairchem-core
```

## Installation

- **Python** 3.8+
- **PACKMOL** (see above)

```bash
# Stable release
pip install mdinterface

# Development version
git clone https://github.com/roncofaber/mdinterface.git
cd mdinterface
pip install -e .
```

Optional extras:

```bash
pip install mdinterface[resp]   # RESP charge analysis
pip install mdinterface[aimd]   # FAIRChem AIMD
pip install mdinterface[all]    # everything
```

## Quick start

### Define species

```python
from mdinterface import BoxBuilder
from mdinterface.database import Water, Ion, Metal111

water = Water(model="ewald")
na    = Ion("Na", ffield="Cheatham")
cl    = Ion("Cl", ffield="Cheatham")
gold  = Metal111("Au")
```

### Build a gold / NaCl electrolyte / gold sandwich

```python
simbox = BoxBuilder(xysize=[15, 15], verbose=True)

simbox.add_slab(gold, nlayers=3)
simbox.add_solvent(water, ions=[na, cl], nions=[5, 5], zdim=25, density=1.0)
simbox.add_slab(gold, nlayers=3)
simbox.add_vacuum(zdim=5)

simbox.build(padding=0.5)
simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True)
```

### Mixed-solvent box (water + methanol)

```python
from mdinterface.core.specie import Specie

methanol = Specie("CH3OH", ligpargen=True)

simbox = BoxBuilder(xysize=[25, 25])
simbox.add_solvent(
    [water, methanol],
    ratio=[3, 1],      # 3 water : 1 methanol by mole
    density=0.95,
    zdim=30,
    ions=[na, cl],
    nions=[5, 5],
)
simbox.build(padding=0.5)
simbox.write_lammps("data_mixture.lammps", atom_style="full", write_coeff=True)
```

### Convert to ASE or MDAnalysis

```python
atoms    = simbox.to_ase()       # ase.Atoms with cell and PBC
universe = simbox.universe       # mda.Universe
```

See the [examples/box_builder/](examples/box_builder/) directory for more complete scripts:

| Script | What it shows |
|--------|--------------|
| `builder_electrode_interface.py` | Au / NaCl electrolyte / Au sandwich |
| `builder_solvent_box.py` | Pure solvent + dissolved species |
| `builder_multisolvent_box.py` | Mixed-solvent box with ratio/density/count modes |
| `builder_multilayer.py` | Five-layer multi-slab system |
| `builder_sandwich_from_traj.py` | Load a relaxed structure via `hijack` |

The legacy `SimulationBox` API is still available and unchanged; see [examples/simulation_box/](examples/simulation_box/).

## Roadmap

Since the original idea was to make a package to build MD boxes layer by layer, I am strongly debating renaming everything as "Workflow for Easy Molecular DYnamics Simulations", aka WEMDYS.

## Questions & Issues

Sir, this is a WEMDY'S. Please contact me or open an issue, glad to talk about ideas and improvements!
