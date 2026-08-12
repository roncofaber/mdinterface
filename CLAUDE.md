# mdinterface - Project Instructions

Python package for building Molecular Dynamics simulation systems - electrolyte/electrode
interfaces, solvent boxes, mixed-solvent electrolytes, polymer networks - with LAMMPS output
(and experimental GROMACS output).

- **Package** `mdinterface` · **Main branch** `main` · Python 3.8+ · License Apache-2.0
- **Stack** ASE, MDAnalysis, NetworkX, NumPy, PACKMOL (external binary); optional PySCF/gpu4pyscf
  (RESP charges), fairchem-core (AIMD), LigParGen + BOSS (OPLS-AA parameter generation)
- **Pipeline** `SimCell` (fluent layer-by-layer builder) → PACKMOL → `ase.Atoms` /
  `MDAnalysis.Universe` → `write_lammps()` / `write_gromacs()`

## Documentation map

This file is the fast-reference layer; `docs/` (mkdocs-material, published to GitHub Pages) is
the deep reference. **Change something a `docs/` file describes, update that file in the same
change.**

| File | Covers | Update when you... |
|---|---|---|
| `docs/guide/simcell.md` | `SimCell` layer-by-layer workflow: `add_slab`/`add_solvent`/`add_vacuum`, `dilate`, `match_cell`, `stack_axis` | Change `SimCell`'s public API or `build()` semantics |
| `docs/guide/specie.md` | `Specie` construction, topology fields, LigParGen/RESP | Change `Specie`'s public API |
| `docs/guide/polymer.md` | Polymer chain building | Change `Polymer`'s public API |
| `docs/guide/database.md` | Built-in species (`Water`, `Ion`, `Metal111`, ...) | Add/rename a database species |
| `docs/guide/logging.md` | Verbosity levels (`verbose=` scale) | Change logging behavior |
| `docs/api/*.md` | mkdocstrings-generated reference (numpy docstring style, `show_source: true`) | Regenerated from docstrings - keep docstrings accurate instead |
| `CHANGELOG.md` | User-facing history | Every user-visible change, same commit |
| `.claude/skills/release` | Version bump, changelog promotion, tag/push, release workflow stages | Change the release process itself |

## Build & test

```bash
pip install -e .[all]      # editable install, all optional extras (resp, aimd)
pytest                     # tests/ (pyproject.toml: testpaths = ["tests"])
mkdocs serve                # preview docs locally
```

PACKMOL must be on `PATH` (`conda install -c conda-forge packmol`). Most of `build/` and every
solvent-layer test shells out to it; failures leave the temp working dir behind for inspection
(logged as a warning) instead of being cleaned up.

## Versioning

`mdinterface/__init__.py`'s `__version__` is the single source of truth (`pyproject.toml` reads
it via `dynamic = ["version"]` / `attr = "mdinterface.__version__"`). `.github/workflows/release.yml`
runs `preflight` → `test` → `build` → `publish` (PyPI) → `release` (draft GitHub release) on a
`v*.*.*` tag push or manual `workflow_dispatch`; `preflight` hard-fails if `v{tag}` doesn't match
`__version__`, or if the changelog summary check (see Changelog discipline) fails - bump the
string and write the changelog entry before tagging, not after. See `.claude/skills/release` for
the full procedure. `docs.yml` redeploys GitHub Pages on push to `main` when `docs/`,
`mdinterface/`, `mkdocs.yml`, or `pyproject.toml` change.

## Changelog discipline

Add a `CHANGELOG.md` entry under `## [Unreleased]` **in the same change** that fixes a bug or
changes behavior, grouped `### Added`/`### Changed`/`### Fixed`. Unlike a consumer-app changelog,
entries here name the actual function/class and the mechanism - this is a library other code
depends on, and "fixed a bug" doesn't tell someone whether an upgrade affects them. One line per
fix, e.g.:

- `NameError` in `Specie.estimate_charges(assign=True)` for non-RESP methods (`atoms` was only
  defined in the RESP branch)
- Wrong return-type annotation on `SimCell._stack_layers` (declared 2-tuple, returns 3-tuple)

At release time (see `.claude/skills/release`), every versioned section (`## [X.Y.Z]`) must also
get a one-to-two sentence summary paragraph before its first `### ` group -
`release.yml`'s `preflight` job runs `scripts/extract_changelog_summary.py` and **fails the
release, before tests or PyPI publish run, if it's missing.**

## Key architecture decisions

- **`SimCell`** (`build/builder.py`) is the current fluent builder: `add_slab`/`add_solvent`/
  `add_vacuum` queue layers, `.build()` assembles them. `BoxBuilder` is a deprecated alias.
- **`SimulationBox`** (`simulationbox.py`) is a separate, older implementation kept only for
  backward compatibility (single interface/enderface/midface slot). It is **not** a subclass of
  `SimCell` and duplicates some of its topology-indexing logic independently - a fix in one place
  does not automatically apply to the other (see gotchas).
- **`Specie`** (`core/specie.py`) bundles an `ase.Atoms` with LAMMPS topology (atom types,
  charges, bonds/angles/dihedrals/impropers, LJ params). Database species and `Polymer` are all
  `Specie`s. Every topology field is optional - a bare `Specie(atoms)` is valid for AIMD/ML-MD
  workflows that never call `write_lammps()`.
- **`database/`** entries (`Water`, `Ion`, `Metal111`, ...) are just pre-built `Specie`
  instances/factories, not a distinct type.
- **PACKMOL is the only packing backend.** `build/box.py::populate_box` shells out to it via a
  temp dir per call.
- **Two independent numbering schemes must agree in LAMMPS output.** Atom "type" is stored as a
  human-readable label string (`Specie.get_atom_types()` → `extended_label`), and
  `io/lammpswriter.py::DATAWriter._write_atoms` renumbers it to a contiguous `1..N` via
  `np.unique(..., return_inverse=True)` *at write time*. Bond/angle/dihedral/improper "type" is a
  raw numeric `.id` assigned once by `_update_topology_indexes` (`build/builder.py`, and a
  duplicate implementation in `simulationbox.py`). These two schemes are not naturally consistent
  with each other - `io/lammpswriter.py::write_lammps_coefficients` renumbers everything against
  what is actually present in the assembled `Universe` at write time rather than trusting `.id`
  directly; do not add a new output path that trusts `.id` as if it were already a clean LAMMPS
  type index.

## Testing

`tests/` roughly mirrors module names, not exact package paths (`test_builder.py`,
`test_solvent.py`, `test_specie.py`, `test_topology.py`, `test_database.py`,
`test_auxiliary.py`). There is currently **no dedicated test file for `io/lammpswriter.py` or
`build/box.py`**; coverage of the LAMMPS output format itself is thin - prefer adding assertions
on file structure (header counts vs. section row counts) over just "does it run" when writing
tests in this area. A clean `write_lammps()` log does not prove the file is correct: `lmp` itself
is the only reliable check for that output format, since `test_builder.py` only checks that
`write_lammps()` doesn't raise, never that the file loads or that section row counts match the
header.

## Optional dependencies

- **LigParGen** (OPLS-AA auto-parameterization) needs a BOSS backend pointed to via
  `~/.config/mdinterface/config.ini` (`[settings] BOSSdir = ...`, path is OS-dependent per
  `platformdirs`) - native dir, Apptainer/Singularity `.sif`, or Docker image tag.
- **RESP charges** need PySCF + PyMBXAS + gpu4pyscf (`pip install mdinterface[resp]`).
- **AIMD** needs `fairchem-core` (`pip install mdinterface[aimd]`).
