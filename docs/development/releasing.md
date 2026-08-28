# Releasing mdinterface

Only maintainers with permission to publish to PyPI and create GitHub releases can complete this workflow.

## 1. Bump the version

`mdinterface/__init__.py` defines `__version__`, which `pyproject.toml` reads dynamically. Set it to `X.Y.Z` without a `v` prefix and update `__date__`.

## 2. Promote the changelog

The `Unreleased` section should already contain the user-visible changes accumulated during development.

- Rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD`.
- Add a one- or two-sentence summary paragraph below the version header and above the first category heading.
- Add a new empty `## [Unreleased]` section above the released version.

The summary is required. Release CI validates it with:

```bash
python scripts/extract_changelog_summary.py X.Y.Z
```

## 3. Verify the release commit

Run the tests, build the documentation, and build the distributions:

```bash
python -m pip install build
pytest -q
mkdocs build --strict
python -m build
python scripts/verify_distribution_contents.py dist
```

Confirm that the wheel and source distribution contain the expected files before publishing.

## 4. Commit and tag

Commit the version and changelog together, push the release commit, then create and push the matching tag:

```bash
git add mdinterface/__init__.py CHANGELOG.md
git commit -m "Bump version to X.Y.Z"
git push origin main
git tag vX.Y.Z
git push origin vX.Y.Z
```

Pushing the tag runs `.github/workflows/release.yml`. The workflow validates the version and changelog, runs tests, builds the distributions, publishes them to PyPI, and creates a draft GitHub release.

## 5. Publish the GitHub release

Review the draft release notes and artifacts on GitHub, then publish the draft manually. PyPI publishing completes before the draft GitHub release is created.

## Re-running a failed release

Run the Release workflow manually and provide the existing `vX.Y.Z` tag. The workflow checks out that tag before validation, testing, and building. Re-run publishing only when the previous attempt did not successfully upload that version to PyPI, because PyPI does not allow replacing an existing release file.
