# mdinterface project instructions

`mdinterface` is a Python package for building molecular dynamics systems, including electrolyte/electrode interfaces, solvent boxes, mixed-solvent electrolytes, and polymer networks. Its primary output is LAMMPS data, with experimental GROMACS support.

- Package: `mdinterface`
- Main branch: `main`
- Supported Python: 3.10-3.14
- License: Apache-2.0
- Core stack: ASE, MDAnalysis, NetworkX, NumPy, and PACKMOL
- Main pipeline: `SimCell` to PACKMOL to `ase.Atoms` or `MDAnalysis.Universe` to `write_lammps()` or `write_gromacs()`

## Working agreements

- Make focused changes and do not refactor unrelated code.
- Preserve user changes in a dirty working tree.
- Add or update tests for behavior changes and bug fixes.
- Update user documentation in the same change when modifying a documented public API or workflow.
- Add a concise `CHANGELOG.md` entry under `Unreleased` for every user-visible change.
- Use NumPy-style docstrings for public Python APIs.
- Do not hard-wrap Markdown paragraphs or list items. Keep each paragraph or item on one source line and let the renderer wrap it.
- Do not place temporary agent plans or scratch files under `docs/`. Use an ignored location such as `.agent-work/`.

## Instruction portability

For machine-readable LAMMPS handoffs, use optional `SimCell.write_lammps(metadata=...)` as documented in `docs/guide/simcell.md`. Record final exported IDs and file checksums, not assumed internal type numbers. Structural metadata must not silently select constraints, force-field styles, boundary conditions or simulation protocols for consumers.

- `AGENTS.md` is the canonical source for repository-wide agent guidance.
- Tool-specific instruction files may exist only as minimal adapters that import or direct the tool to `AGENTS.md`.
- Do not duplicate shared instructions in tool-specific files because duplicated guidance will drift.
- Keep project conventions and workflows in `AGENTS.md` or regular project documentation. Keep tool-specific files limited to capabilities or settings that have no portable equivalent.

## Development commands

Install the core package and test dependency:

```bash
pip install -e .[test]
```

Create the full parameterization development environment when working on LigParGen, BOSS, or OpenFF:

```bash
mamba env create -f environment-full.yml
mamba activate mdinterface-full
```

Install and preview the documentation:

```bash
pip install -r docs/requirements.txt
mkdocs serve
```

Verify a change with the relevant subset first, then run the complete checks when practical:

```bash
pytest -q
pytest -q -m integration
mkdocs build --strict
```

The editable install provides the upstream PACKMOL package and executable used by packing workflows and tests. Tests that execute PACKMOL use the `integration` marker. Optional RESP and AIMD dependencies are heavy and should only be installed for work involving those features.

## Documentation map

| File | Covers | Update when changing |
|---|---|---|
| `docs/guide/simcell.md` | `SimCell` layer construction and build behavior | `SimCell` public API or build semantics |
| `docs/guide/specie.md` | `Specie`, topology fields, LigParGen, and RESP | `Specie` public API |
| `docs/guide/polymer.md` | Polymer construction | `Polymer` public API |
| `docs/guide/database.md` | Built-in species | Database species or factories |
| `docs/guide/logging.md` | Verbosity behavior | Logging behavior |
| `docs/api/*.md` | Generated API reference | Public docstrings and exports |
| `CONTRIBUTING.md` | Development setup and contribution workflow | Development or verification workflow |
| `docs/development/releasing.md` | Versioning, changelog promotion, tags, and publishing | Release workflow |

## Architecture constraints

- `SimCell` in `mdinterface/build/builder.py` is the current fluent builder. `BoxBuilder` is a deprecated alias.
- `SimulationBox` in `mdinterface/simulationbox.py` is a separate deprecated implementation retained for compatibility. It duplicates some topology-indexing behavior, so fixes to `SimCell` do not automatically apply to it.
- `Specie` combines `ase.Atoms` with optional LAMMPS topology and force-field data. A bare `Specie(atoms)` is valid for workflows that do not write LAMMPS data.
- Database entries such as `Water`, `Ion`, and `Metal111` are pre-built `Specie` objects or factories, not a separate domain type.
- PACKMOL is the only packing backend. `mdinterface/build/box.py::populate_box` runs it in a temporary directory.
- LigParGen is the supported automatic OPLS-AA/CM1A parameterization backend. Preserve `ligpargen=True` compatibility when introducing a backend-neutral API.
- OpenFF-to-LAMMPS export and re-import have been validated, but OpenFF is not yet a supported `Specie` backend. The current topology and GROMACS paths assume OPLS-style torsions, while OpenFF LAMMPS export uses Fourier proper torsions and CVFF impropers.
- Any OpenFF integration must preserve bond, angle, dihedral, and improper styles together with mixing rules, switching behavior, and van der Waals and Coulomb 1-4 scaling. Do not silently combine molecular species with incompatible global LAMMPS conventions.
- LAMMPS atom types and bonded interaction types use independent numbering paths. Atom labels are renumbered in `DATAWriter._write_atoms`, while bonded type IDs originate in `_update_topology_indexes`. Output code must renumber against the interactions present in the assembled universe rather than assume stored IDs are contiguous LAMMPS type numbers.

## Testing risks

Tests roughly mirror module names rather than package paths. Coverage of `mdinterface/io/lammpswriter.py`, `mdinterface/build/box.py`, and `mdinterface/simulationbox.py` is comparatively thin. For LAMMPS output changes, assert file structure and consistency between header counts and section rows. When available, loading the result with LAMMPS is stronger verification than checking only that `write_lammps()` did not raise.

## Changelog and releases

`CHANGELOG.md` is written for package users, not as a complete repository activity log. Include public features, behavior changes, bug fixes, deprecations, installation changes, and compatibility changes. Omit internal refactoring, agent instructions, contributor-only documentation, CI maintenance, and release-process changes unless they alter what users install or receive.

Group entries under `Added`, `Changed`, `Fixed`, or `Deprecated`. Use one concise sentence per bullet and cover one observable change per bullet. Lead with the affected API when practical. Include implementation details only when they explain compatibility or user-visible behavior. Do not combine unrelated changes with semicolons, and do not hard-wrap bullets in the source file.

Versioned changelog sections require a one- or two-sentence summary paragraph before their first category heading. Release CI enforces this with `scripts/extract_changelog_summary.py`.

Follow `docs/development/releasing.md` for releases. Do not commit, push, tag, publish packages, or publish a GitHub release unless the user explicitly requests those actions.

## Optional integrations

- LigParGen requires a BOSS backend configured through the platform-specific `mdinterface` configuration directory.
- `environment-full.yml` provides the reproducible Python 3.12 and NumPy 1.x development stack required by the current AmberTools dependencies for LigParGen and OpenFF work. The core package separately supports Python 3.10-3.14. The environment installs the public LigParGen fork, AmberTools, Open Babel, OpenFF Toolkit and Interchange, and CPU-only NAGL and PyTorch, but never BOSS.
- The `boss-container` repository distributes build recipes only. User-built images contain the user's licensed BOSS files and must not be treated as redistributable project artifacts.
- OpenFF packages are currently distributed through conda-forge rather than PyPI. Keep OpenFF optional and do not add it to the core pip dependency set.
- RESP work requires PySCF, PyMBXAS, and gpu4pyscf. The GPU package is platform-specific and is not installed by the `resp` extra.
- AIMD work requires `fairchem-core` through the `aimd` extra.
