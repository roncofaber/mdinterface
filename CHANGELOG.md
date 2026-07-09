# Changelog

All notable changes to mdinterface are documented here.

## [Unreleased]

---

## [1.5.3] — 2026-07-09

### Fixed
- `NameError` in `DATAWriter` when `convert_units=False` (`coordinates` and `triv` were only assigned inside the conditional branch)
- Improper type label in LAMMPS coefficient output was silently truncated to the first atom; now writes the full four-atom type string
- `NameError` in `Specie.estimate_charges(assign=True)` for non-RESP methods (`atoms` was only defined in the RESP branch)
- Wrong return-type annotation on `SimCell._stack_layers` (declared 2-tuple, returns 3-tuple)
- `map_impropers(None)` returned a 2-tuple instead of the consistent 3-tuple returned by all other `map_*` functions
- Dead unreachable error-checking code removed from `generate_missing_interactions`

### Changed
- Canonical repository moved from GitLab to GitHub
- Deployment: tag/version consistency check and GitHub Release creation added to publish workflow
- Deployment: docs workflow now only rebuilds when source files change

---

## [1.5.2] — 2026-03-31

### Fixed
- Topology label mismatch for large molecules in LigParGen segmentation

---

## [1.5.1] — 2026-03-31

### Added
- GitHub Actions workflow for automated PyPI publishing on version tags
- GROMACS output: `write_gromacs_itp`, `write_gromacs_top`, `Specie.write_gro()`
- `refine_large_specie_topology` for LigParGen parameterisation of molecules with more than 200 atoms
- `Specie.write_gro()` for writing GROMACS structure files
- `BOSSdir` config key supporting local install, container, and direct-path modes

### Fixed
- Junction LJ type correction in segment and polymer refinement
- Docs CI: use `--no-deps` to avoid building compiled dependencies (libarvo)
- PACKMOL and LigParGen now use a temp directory; kept on failure for inspection

### Changed
- `libarvo` made optional; only required for volume/radius estimation
- Logging style updated throughout; improved SimCell and externals API docs

---

## [1.5.0] — 2026-02-01

### Added
- MkDocs documentation site with Material theme
- Polymer guide and reorganised examples
- `SimCell` fluent builder API (replaces `BoxBuilder`)
- Multi-solvent support with `ratio` and mixed `nsolvent` lists
- `solvent.py` extracted from `build.py` for clarity
- `dilate` and `packmol_tolerance` parameters on `add_solvent`
- `hijack`, `stack_axis`, `match_cell` options on `SimCell.build()`
- Logging infrastructure with `set_verbosity` and structured headers

### Changed
- `ions`/`nions` renamed to `solute`/`nsolute` throughout
- `populate_with_ions` renamed to `populate_solutes`
- `BoxBuilder` retained as a deprecated alias for `SimCell`
- Examples reorganised into `simulation_box/` and `box_builder/` subfolders

### Fixed
- Duplicate log output caused by missing `propagate=False`
- Mutable default arguments in several functions

---

## [1.4.0]

### Added
- Initial `BoxBuilder` fluent API for multi-layer simulation box assembly
- `to_ase()` output method
- Noble gases in database

### Fixed
- Shell injection and insecure temp file naming
- Bare `except` clauses replaced throughout
- Raise-string bugs fixed
- Improper detection corrected

---
