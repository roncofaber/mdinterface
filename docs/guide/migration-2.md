# Migrating to mdinterface 2.0

Version 2.0.0 removes the deprecated `SimulationBox` implementation, the `mdinterface.simulationbox` module, and the `BoxBuilder` alias. `SimCell` is the supported system builder. `Specie`, `Polymer`, and the database force-field models remain available. Python 3.10 or newer is required.

## Molecular graphs and SMILES

RDKit is now a required dependency. Use `Specie(smiles="CCO")` for explicit SMILES input or pass an RDKit molecule through `atoms`. Positional strings retain their ASE interpretation. SMILES formal charges define the molecular charge, but simulation partial charges still require force-field parameterization. Polymer assembly now validates monomer chemistry and preserves explicit junction bonds. RDKit-backed LigParGen parameterization requires the compatible fork with chemistry-preserving atom reordering; see the [Specie guide](specie.md).

## Replace the old builders

For `BoxBuilder`, replace the import and constructor name with `SimCell`; its layer-building methods already match the current API.

For `SimulationBox`, construct `SimCell(xysize=...)`, add each layer explicitly, then build and export the system:

```python
from mdinterface import SimCell
from mdinterface.database import Water, Ion, Metal111

water = Water()
gold = Metal111("Au")
box = SimCell(xysize=[15, 15])
box.add_slab(gold, nlayers=1)
box.add_solvent(
    water, zdim=25, density=1.0,
    solute=[Ion("Na"), Ion("Cl")], nsolute=[5, 5],
)
box.add_slab(gold, nlayers=1)
box.add_vacuum(zdim=5)
box.build(padding=0.5, match_cell=False)
box.write_lammps("data.lammps")
atoms = box.to_ase()
universe = box.universe
```

| Previous workflow | Current workflow |
| --- | --- |
| `interface`, `miderface`, or `enderface` layer dictionaries | `add_slab(species, nlayers=...)` in the desired order |
| Solvent layer dictionaries | `add_solvent(solvent, solute=..., nsolute=..., zdim=...)` |
| Layer keys `nions` and `ion_pos` | `nsolute` and `solute_pos` |
| Vacuum layer dictionary | `add_vacuum(zdim=...)` |
| `make_simulation_box(...)` | `build(...)`, followed by `universe` or `to_ase()` |
| `write_data=True` or `write_lammps_file(...)` | `write_lammps(...)` after `build()` |
| `center_electrode=True` | Review `build(center=True)` and the coordinate change below |

`SimCell.build()` returns `None`. Obtain the assembled structure through `box.universe` or `box.to_ase()`. Set `padding` and `match_cell` explicitly when migrating: the legacy builder defaulted to `padding=1.5` and `match_cell=False`, whereas `SimCell` defaults to `0.5` and `True`.

Legacy examples have been removed. The current examples cover solvent boxes, mixed solvents, electrode interfaces, multiple layers, spatial pockets, trajectory-based membranes, and polymer assembly.

## Review coordinate-dependent settings

In 1.x, `SimCell.build(center=True)` centered the first layer on the periodic boundary. In 2.0 it centers that layer at the box midpoint, using its allocated thickness, then wraps individual atoms. For a 100 Å box with a 20 Å first layer, that layer moves from 90-100 Å and 0-10 Å to 40-60 Å.

Update spatial selections, restraints, and analysis bins that depend on absolute coordinates. The default `center=False` is unchanged. The periodic seam need not fall in vacuum, and other molecules can be split across the boundary. See the [SimCell guide](simcell.md) for `stack_axis` and `hijack` behavior.

## Handle external-tool failures

PACKMOL is installed with the core package. Failed execution or unreadable output now raises `PackmolError` instead of returning `None`. LigParGen failures raise `LigParGenError`; both exceptions provide diagnostic information. LigParGen and a licensed BOSS installation remain optional requirements for automatic parameterization.

The low-level `populate_box()` API accepts `Specie` templates and returns `ase.Atoms`, replacing its previous MDAnalysis-based interface. Use `SimCell.universe` when you need the assembled MDAnalysis topology.

## Validate region content

Nested regions must fit inside their parent's actual shape. Both bulk solvent and bulk solute exclude pockets, and density and concentration calculations use the remaining volume. Count lists must match the species lists; mixtures require explicit per-species counts or a mixing ratio. `regions` cannot be combined with `conmodel` or `solute_pos="center"`; use a region fill for confined solutes.

## Polymer preparation API

Replace `chain.refine_polymer_topology(Nmax=12, offset=True, ending="H")` with `chain.refine_junctions(snippet_radius=12, charge_correction="uniform", cap_element="H")`. The default correction is `"none"`. The method returns a charge audit and applies changes only after all junction calculations succeed. Constructor-driven refinement still uses `refine_polymer=True`, with the renamed `charge_correction` and `cap_element` options. Use the public `chain.junction_bonds` property instead of `_get_connection_elements()`.

Geometry preparation is available on both `Specie` and `Polymer`: `generate_conformer(seed=42)` embeds a new conformer, `minimize_geometry()` minimizes existing coordinates, and `generate_conformer(seed=42, minimize=True)` combines both operations. Charge-estimation helpers now default to the stored molecular charge and reject conflicting overrides instead of assuming zero.

## Parameter transfer and validation

Use `Specie.parameterize()` to assign a validated LigParGen result without manually unpacking parameter lists. `Polymer` inherits parameters automatically from `Specie` monomers unless explicit topology overrides are provided. `mark_attachment_sites(head_map=..., tail_map=...)` replaces manual index-based attachment annotations for mapped molecules.

`refine_large_topology()` now uses `snippet_radius`, `cap_element`, and `charge_correction`, matching `refine_junctions()`. Automatic large-molecule parameterization no longer applies an implicit uniform charge correction. Request `charge_correction="uniform"` explicitly and inspect the returned audit. Segment sizes include caps, and oversize junctions are rejected before backend execution.

Classical exports with coefficients now reject incomplete or nonfinite parameters before opening the output file. Explicit zero coefficients remain valid and are retained in the topology. Use `expected_charge=0.0` to require neutrality; charged systems remain supported. Imported LAMMPS bonds now remain fixed when positions change.
