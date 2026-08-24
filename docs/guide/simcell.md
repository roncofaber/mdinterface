# SimCell

`SimCell` is the main class for building MD simulation boxes. Add layers in the order they should appear along the z axis (you can rotate it later), then call `.build()` to assemble and pack everything.

```python
from mdinterface import SimCell
```

## Creating a SimCell

```python
simbox = SimCell(xysize=[15, 15], verbose=True)
```

| Parameter | Description |
|-----------|-------------|
| `xysize` | Target XY dimensions in Angstroms. Used as a tiling guide for slabs; the actual cell XY is set by the largest fitted slab. |
| `verbose` | Verbosity: `False`/`0` = quiet, `True`/`1` = INFO, `2` = DEBUG, or a string like `"DEBUG"`. |

## Adding layers

You can add any type of layer to your system with the following methods. Call them in the order you want them along the stacking axis.

### `add_slab(species, nlayers)`

Adds a periodic solid layer by tiling a unit-cell `Specie` in XY. The number of repeats is chosen to cover the requested `xysize`; the actual XY cell is then set by the tiled footprint (an integer multiple of the unit cell). `nlayers` controls how many unit-cell layers are stacked along Z.

```python
from mdinterface.database import Metal111
gold = Metal111("Au")
simbox.add_slab(gold, nlayers=4)
```

When multiple slabs are present, each gets distinct atom-type labels (`_s0`, `_s1`, …) so their force-field parameters remain independent.

### `add_prebuilt(species, nlayers)`

Like `add_slab`, but skips XY tiling — the species positions are used exactly as provided. Intended for structures taken from a previous MD run or geometry optimisation (e.g. a relaxed electrode loaded with `hijack` and then saved back as a `Specie`).

```python
simbox.add_prebuilt(relaxed_electrode, nlayers=1)
```

Any centering or trimming should be done on the species before passing it here.

### `add_solvent(solvent, zdim, ...)`

Adds a PACKMOL-packed solvent (liquid) layer. Accepts a single `Specie`, a list of species for mixtures, or `None` for a solute-only region.

```python
from mdinterface.database import Water, Ion
water = Water(model="ewald")
na    = Ion("Na", ffield="Cheatham")

simbox.add_solvent(
    water,
    zdim=25,
    density=1.0,
    solute=[na],
    nsolute=5,
    solute_pos="right",   # "left", "right", "center", or None (full box)
)
```

**Mixing solvents:**

```python
simbox.add_solvent(
    [water, methanol],
    ratio=[3, 1],     # molar ratio
    density=0.95,
    zdim=30,
)
```

**Spatially heterogeneous regions (pockets):**

Carve a `Sphere`, `Box`, or `Cylinder` out of the layer and fill it with different content - e.g.
an argon gas pocket embedded in bulk water:

```python
from mdinterface.build.regions import Sphere
from mdinterface.database import Argon

argon  = Argon()
bubble = Sphere(center=(7.5, 7.5, 15.0), radius=5.0)

simbox.add_solvent(
    water,
    zdim=30,
    density=1.0,
    regions=[bubble.fill(argon, density=0.8)],
)
```

A `FilledRegion` can itself carry `regions=[...]` for further nested sub-regions, to arbitrary depth - e.g. a neon bubble nested inside the argon pocket:

```python
from mdinterface.build.regions import Sphere
from mdinterface.database import Argon, Neon

argon  = Argon()
neon   = Neon()
pocket = Sphere(center=(7.5, 7.5, 15.0), radius=5.0)
bubble = Sphere(center=(7.5, 7.5, 15.0), radius=1.0)

simbox.add_solvent(
    water,
    zdim=30,
    density=1.0,
    regions=[pocket.fill(argon, density=0.8, regions=[bubble.fill(neon, density=0.1)])],
)
```

`regions` isn't limited to one nested pocket - multiple independent, non-nested regions carved
out of the same layer are just as common, e.g. two separate gas pockets on opposite sides of the
box:

```python
from mdinterface.build.regions import Sphere
from mdinterface.database import Argon, Krypton

argon   = Argon()
krypton = Krypton()

simbox.add_solvent(
    water,
    zdim=30,
    density=1.0,
    regions=[
        Sphere(center=(4.0, 7.5, 15.0), radius=3.0).fill(argon, density=0.8),
        Sphere(center=(11.0, 7.5, 15.0), radius=3.0).fill(krypton, density=0.8),
    ],
)
```

A region can also confine a solute rather than swap the solvent - e.g. pinning a handful of Na+
ions inside a box-shaped subvolume near one edge instead of letting them roam the whole layer:

```python
from mdinterface.build.regions import Box

corner = Box.from_bounds(0.0, 0.0, 5.0, 6.0, 6.0, 15.0)

simbox.add_solvent(
    water,
    zdim=30,
    density=1.0,
    regions=[corner.fill(solute=[na], nsolute=5)],
)
```

`Region.fill()` accepts the same content parameters as `add_solvent` itself (`solvent`, `solute`,
`nsolute`, `density`, `nsolvent`, `concentration`, `ratio`), minus `zdim`/`xysize` since the
region defines its own extent. `conmodel` is not yet supported inside a region. The bulk solvent's
molecule count automatically accounts for the volume carved out by each region.

A region's `center` can be `"random"` instead of a coordinate tuple, to have a non-overlapping
placement chosen automatically within its parent volume (rejection-sampled against the parent's
bounds and any sibling regions). Pass `seed` to `add_solvent` to make the placement reproducible:

```python
simbox.add_solvent(
    water,
    zdim=30,
    density=1.0,
    regions=[Sphere(center="random", radius=3.0).fill(argon, density=0.8)],
    seed=42,
)
```

**Key parameters:**

| Parameter | Description |
|-----------|-------------|
| `solvent` | `Specie`, list of `Specie`, or `None` (solute-only) |
| `zdim` | Z thickness of the layer in Angstroms |
| `density` | Target density in g/cm³ |
| `nsolvent` | Explicit molecule count (int or list) |
| `ratio` | Molar ratio for multi-solvent mixtures |
| `solute` | `Specie` or list of `Specie` to dissolve |
| `nsolute` | Molecule count per solute species |
| `solute_pos` | Placement region: `None`, `"center"`, or a `Region` instance. `"left"`/`"right"` strings are deprecated - pass an equivalent `Box` instead (e.g. `Box.from_bounds(...)`) |
| `regions` | List of `FilledRegion` (from `SomeRegion(...).fill(...)`) for spatially heterogeneous sub-regions within the layer |
| `seed` | RNG seed for resolving `center="random"` regions, for reproducible placement |
| `dilate` | Expand the PACKMOL box by this fraction to help convergence |
| `packmol_tolerance` | PACKMOL distance tolerance in Angstroms (default 2.0) |

### `add_vacuum(zdim)`

Adds an empty vacuum gap.

```python
simbox.add_vacuum(zdim=10)
```

## Building

```python
simbox.build(padding=0.5, center=False, stack_axis="z")
```

| Parameter | Description |
|-----------|-------------|
| `padding` | Extra space (Angstroms) added above/below each PACKMOL region |
| `center` | Shift the system so the center of the first layer falls in the middle of the box, keeping it intact. Whatever ends up opposite it (typically a vacuum gap, if present) absorbs the periodic seam instead. |
| `layered` | Keep per-layer residue numbering instead of merging |
| `match_cell` | `True` (default): stretch all slabs to the largest XY cell, ensuring a consistent solid/liquid interface. `False`: each slab keeps its natural tiled XY. A `Specie`: lock XY to that species' cell and stretch everything else to match — useful when a pre-relaxed slab or polymer defines the cell. |
| `hijack` | Replace positions and cell with a prebuilt `ase.Atoms` object |
| `stack_axis` | Stacking direction: `"z"` (default), `"x"`, or `"y"` |

### Stacking axis

By default the system is assembled along Z. Use `stack_axis` to permute coordinates at the end:

```python
simbox.build(stack_axis="x")   # stacking direction becomes X
```

This is a coordinate permutation applied after assembly, but the build process itself always runs along Z.

## Output

```python
# LAMMPS data file
simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True)

# ASE Atoms object
atoms = simbox.to_ase()

# MDAnalysis Universe
universe = simbox.universe
```

### GROMACS output *(experimental)*

!!! warning
    GROMACS output is experimental.  Unit conversions and dihedral mappings
    have been validated for OPLS-AA molecules generated by LigParGen, but
    other force fields or atom styles may need adjustments.  Verify the output
    against a reference before production use.

`write_gromacs()` produces three files in the directory of your choice:

- `{prefix}.gro` — atomic coordinates and box vectors
- `{resname}.itp` — per-species include topology (one per unique species)
- `{prefix}.top` — system topology with molecule counts matching the GRO file

```python
simbox.write_gromacs(prefix="system", outdir="gromacs_input/")
```

Individual ITP files can also be written directly from a `Specie`:

```python
water.write_gromacs_itp("water.itp")
```

## Verbosity

```python
import mdinterface
mdinterface.set_verbosity("DEBUG")   # package-wide
mdinterface.set_verbosity(1)         # INFO
mdinterface.set_verbosity(0)         # quiet
```

See the [Logging guide](logging.md) for details.
