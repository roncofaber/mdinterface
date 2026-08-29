# Contributing to mdinterface

## Development setup

Create and activate a Python environment, then install the package with its test dependencies in editable mode. This also installs the upstream PACKMOL package and executable:

```bash
pip install -e .[test]
```

Install optional dependencies only when working on the corresponding feature:

```bash
pip install -e .[resp]
pip install -e .[aimd]
```

RESP fitting also requires the platform-specific `gpu4pyscf` package, which is not installed by the `resp` extra.

Install `libarvo` separately when working on molecular volume or radius estimation.

## Tests

Run the relevant test module while developing, then run the complete suite before submitting a change:

```bash
pytest -q tests/test_builder.py
pytest -q
```

Tests marked `integration` execute PACKMOL. Run only unit tests with `pytest -q -m "not integration"` or only integration tests with `pytest -q -m integration`. Changes to LAMMPS output should verify section structure and counts, and should load the generated data with LAMMPS when that executable is available.

Generate an informational coverage report with `pytest -q --cov=mdinterface --cov-report=term-missing`. The project does not currently enforce a coverage threshold.

## Documentation

Install the documentation dependencies and build the site with warnings treated strictly:

```bash
pip install -r docs/requirements.txt
mkdocs build --strict
```

Update the relevant guide when changing a documented public API or workflow. API reference pages are generated from NumPy-style docstrings, so keep public docstrings current.

## Changelog

Add every user-visible change to `CHANGELOG.md` under `Unreleased`, grouped under `Added`, `Changed`, `Fixed`, or `Deprecated`. Write one concise, unwrapped sentence per bullet and describe one observable change. Omit internal refactoring, agent instructions, contributor-only documentation, CI maintenance, and release-process changes unless they affect what package users install or receive.

See `docs/development/releasing.md` for the release process.
