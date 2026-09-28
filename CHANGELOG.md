# Changelog

All notable changes to mdinterface are documented here.

## [Unreleased]

## [2.0.1] - 2026-09-28

Fixes GROMACS topology generation and validation for single- and multi-species systems, with LAMMPS/GROMACS energy and force regression coverage. Charge audits now report numerical noise near zero as `0.0` without changing atomic partial charges.

### Fixed
- Parameterization and polymer charge audits normalize numerical noise within `1e-12 e` of zero without altering atomic partial charges.
- GROMACS exports restore scaled 1-4 interactions using explicit pairs derived from bond connectivity, including zero-valued torsions.
- `SimCell.write_gromacs()` writes shared atom types before molecule definitions so mixed-species topologies can be preprocessed.
- GROMACS system topologies preserve consecutive molecule order to match the coordinate file.
- GROMACS exports reject incomplete or unsupported parameters and conflicting species definitions before writing files.
- `SimCell.write_gromacs()` rejects assembled residues whose atom types, charges, masses, or atom counts differ from their species definitions.

## [2.0.0] - 2026-09-28

Makes `SimCell` the sole system builder, removing the deprecated `SimulationBox` and `BoxBuilder` APIs and adding spatially constrained solvent regions and optional LAMMPS structural metadata. This major release also changes explicit centering behavior, requires Python 3.10 or newer, and adds RDKit molecular preparation, atomic force-field parameterization, export validation, and improved external-tool diagnostics.

### Added
- `Specie.parameterize()` applies complete LigParGen results atomically and returns a charge audit for direct or segmented parameterization.
- `Specie.mark_attachment_sites()` selects validated polymer leaving atoms using atom maps.
- `Specie.validate_force_field()` and coefficient-bearing LAMMPS exports check parameter completeness, with optional `expected_charge` validation at export.
- `Specie.generate_conformer()` and `minimize_geometry()`, inherited by `Polymer`, provide separate or combined RDKit embedding and MMFF94 minimization while preserving simulation parameters.
- A mapped-SMILES polymer example prepares and parameterizes a charged chain and exports it with counterions.
- `Specie(smiles=...)` and RDKit molecule inputs preserve chemical graphs, formal charges, stereochemistry, and atom maps, with `to_rdkit()` and `to_smiles()` conversion.
- `SimCell.write_lammps(metadata=...)` optionally exports versioned structural JSON with final atom/type IDs, species groups, connectivity, coefficients and a data-file checksum, without choosing simulation settings.
- `Sphere`, `Box`, and `Cylinder` regions for spatially constrained PACKMOL placement
- `Region.fill()` for assigning content to a region
- `SimCell.add_solvent(regions=...)` for spatially heterogeneous solvent layers
- Nested regions through `Region.fill(regions=...)`
- Reproducible automatic region placement through `center="random"` and `SimCell.add_solvent(seed=...)`

### Changed
- `Polymer` inherits force-field parameters from `Specie` monomers with distinct labels per repeat unless explicit topology overrides are supplied.
- Large-molecule parameterization uses explicit charge-correction options and checks capped segment and expanded junction sizes before running LigParGen.
- `Polymer.refine_junctions(snippet_radius=..., charge_correction=..., cap_element=...)` replaces `refine_polymer_topology()` and returns a charge audit, with public `junction_bonds` for connectivity.
- The piperion example prepares geometry with RDKit and prints the refinement charge audit instead of recommending a neutral-charge ML relaxation for an ionic chain.
- RDKit is now a required core dependency for molecular graph handling and SMILES support.
- `Polymer` preserves explicit junction bonds and validates monomer chemistry instead of deriving junction connectivity from coordinates.
- RDKit-backed LigParGen inputs use MOL files, and installation instructions and the full environment pin the compatible fork to commit `ad78036` to preserve chemistry during atom reordering.
- `SimulationBox`, its module `mdinterface.simulationbox`, and the `BoxBuilder` alias are removed; migrate to `SimCell` using the 2.0 migration guide.
- Legacy examples are removed from the repository and source distribution; current examples use `SimCell` and `Polymer`.
- `SimCell.build(center=True)` now places the first layer at the box midpoint instead of across the periodic boundary; the default remains `False`.
- Minimum supported Python version is now 3.10
- Supported Python versions are now tested through Python 3.14
- PACKMOL is now installed automatically from its upstream PyPI package
- Importing `mdinterface` no longer loads optional AIMD and plotting dependencies or reads user configuration
- PACKMOL execution and output failures now raise `PackmolError` with retained diagnostic file locations instead of returning `None`
- LigParGen setup, execution, and output failures now raise `LigParGenError` with actionable installation or configuration guidance and retained diagnostic file locations
- `build.box.populate_box()` now accepts `Specie` objects and returns `ase.Atoms`
- PACKMOL templates now use ASE instead of MDAnalysis, eliminating PDB-completeness warnings for temporary files

### Fixed
- Source distributions include the full parameterization environment and MkDocs configuration referenced by the bundled setup instructions.
- SMILES parsing preserves explicit hydrogen atom maps used to select polymer attachment sites.
- `Specie.repeat()` scales the stored molecular charge and preserves the original copy's atom maps.
- Junction refinement retains distinct improper assignments with identical coefficients.
- Database ions retain their formal molecular charge independently of scaled force-field partial charges.
- LAMMPS coefficient export uses unique temporary files instead of overwriting `tmp_data.lammps` in the working directory.
- Explicit zero-valued bonded parameters are retained instead of being discarded during type mapping.
- Imported LAMMPS bond connectivity survives trajectory coordinate updates, atom reordering, and `Specie.repeat()` without distance-based bond inference.
- Junction refinement leaves the original chain unchanged if any parameterization fails.
- Charge estimation uses the stored molecular charge, rejects conflicting overrides, and passes the charge through the Open Babel and RESP backends.
- `Specie(ligpargen=True)` now retains the partial charges returned by parameterization.
- Polymer refinement preserves charged groups, uses neutral snippet caps, and no longer reuses parameters based only on element strings.
- The piperion example locates its charged nitrogen correctly and validates charge neutrality after refinement.
- Nested regions now remain inside their parent's actual shape, including randomly placed regions.
- `Region.fill()` validates content parameters, and solvent mixtures reject ambiguous scalar counts without a mixing ratio instead of silently omitting species.
- Bulk solutes now exclude filled regions and use the remaining volume for concentration-based counts; incompatible fixed-center and concentration-profile placement raises an error.
- Water model documentation now describes `Water(model="ewald")` as modified TIP3P for Ewald electrostatics (its parameters were always TIP3P-Ewald, not SPC/E) and no longer swaps variable names in the database guide
- Slab tiling producing cells smaller than the requested XY dimensions when the nearest repeat count rounded down
- Spurious MDAnalysis topology-guessing warnings in `Specie.to_universe()` and `build.box.populate_box()`
- Missing `elements` topology data in universes created by `Specie.to_universe()`
- Region-constrained PACKMOL placements exceeding requested boundaries because of loose default solver precision
- Wheels including repository documentation, examples, tests, and local planning files because package discovery was not limited to `mdinterface`
- Missing `mdinterface/config.ini` in wheels, which broke the default configuration fallback after installation

### Deprecated
- `add_solvent(solute_pos="left"/"right")` - pass an equivalent `Region`, such as `Box.from_bounds(...)`, instead

---

## [1.5.4] - 2026-08-12

Automates PyPI releases with a test gate and a required changelog summary, and fixes `SimulationBox` crashing on its most common usage (density-based single solvent) while making its deprecation warning actually visible.

### Changed
- Release process: `publish.yml` replaced by `release.yml` - now runs the test suite and requires a changelog summary paragraph before publishing to PyPI, and creates a draft (not published) GitHub release with changelog-sourced notes instead of autogenerated commit notes
- `.github/workflows/release.yml` also accepts manual `workflow_dispatch` (with a `tag` input) for re-running a release without re-tagging

### Added
- `.claude/skills/release` documenting the version-bump/changelog/tag procedure
- `scripts/extract_changelog_summary.py` for extracting a version's changelog summary paragraph

### Deprecated
- `SimulationBox` now emits a `DeprecationWarning` on construction; use `SimCell` instead

### Fixed
- `simulationbox.py` called `warnings.filterwarnings('ignore')` at import time, silently suppressing all warnings process-wide for the rest of the program (including the existing `BoxBuilder` deprecation warning and this module's own warnings) - removed
- `SimulationBox` crashed with `AttributeError: AtomGroup has no attribute get_masses` when using a single solvent with density (no `ratio`/`nsolvent`) - its most common usage

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
