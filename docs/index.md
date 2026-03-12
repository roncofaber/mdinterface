<div style="display: flex; align-items: center;">
  <img src="assets/mdinterface.png" alt="Logo" width="80" style="margin-right: 14px;">
  <div>
    <h1 style="margin: 0;">mdinterface</h1>
  </div>
</div>

`mdinterface` is a Python package for building systems for Molecular Dynamics (MD) simulations. Initially developed for electrolyte/electrode solid-liquid interfaces, it is equally suited for pure solvent boxes, mixed-solvent electrolytes, and polymer networks.

## Features

- **Layer-by-layer `SimCell` builder**: add slabs, solvent regions, and vacuum gaps one step at a time; call `.build()` when done.
- **ASE & MDAnalysis integration**: the assembled box converts to `ase.Atoms` or `mda.Universe` with a single call, ready for any downstream tool.
- **Multi-solvent support**: mix solvents by molar ratio + density, ratio + total count, or explicit per-species molecule counts.
- **Ion placement**: dissolve ions by count, molar concentration, or a spatially-varying concentration profile.
- **PACKMOL integration**: handles molecular packing automatically; tolerance and dilation are tunable per layer.
- **Configurable stacking axis**: build along Z (default) and permute to X or Y at the end.
- **Polymer builder**: generate chains of arbitrary length from a monomer `Specie`.
- **AIMD with FAIRChem**: run ML-potential dynamics via FAIRChem (optional).
- **RESP charges**: estimate partial charges with PySCF / gpu4pyscf (optional).
- **Force-field database**: pre-defined parameters for common metals, noble gases, water models, and ions; or generate OPLS-AA parameters on the fly with LigParGen.
- **LAMMPS output**: writes data files and force-field coefficient blocks ready to run.

!!! warning
    The documentation is still under construction and *might* contain mistakes.

## Quick example

```python
from mdinterface import SimCell
from mdinterface.database import Water, Ion, Metal111

water = Water(model="ewald")
na    = Ion("Na", ffield="Cheatham")
cl    = Ion("Cl", ffield="Cheatham")
gold  = Metal111("Au")

simbox = SimCell(xysize=[15, 15], verbose=True)
simbox.add_slab(gold, nlayers=3)
simbox.add_solvent(water, solute=[na, cl], nsolute=[5, 5], zdim=25, density=1.0)
simbox.add_slab(gold, nlayers=3)
simbox.build(padding=0.5)
simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True)
```

See the [Quick Start](quickstart.md) page for more examples, or jump straight to the [User Guide](guide/simcell.md).
