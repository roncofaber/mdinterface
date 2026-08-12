---
name: release
description: Cut a tagged release of mdinterface - version bump, changelog promotion, tag and push. Use this whenever asked to cut, ship, or publish a release, or bump the version. Do NOT improvise these steps from memory; the changelog summary paragraph is enforced by CI and a missing one blocks the PyPI publish.
disable-model-invocation: true
---

# Releasing mdinterface

## 1. Bump the version

`mdinterface/__init__.py`'s `__version__` is the single source of truth (`pyproject.toml` reads it dynamically). Bump it to `X.Y.Z`, no `v` prefix.

## 2. Promote the changelog

`## [Unreleased]` should already be populated incrementally as changes land (see CLAUDE.md's changelog discipline - never backfilled from git history at release time).

- Rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD`
- Add a one- or two-sentence summary paragraph directly below the header, above the first `### ` group
- Add a fresh, empty `## [Unreleased]` above it

The summary paragraph is not optional: `.github/workflows/release.yml`'s `preflight` job runs `scripts/extract_changelog_summary.py` and **fails the entire release - before tests run, before PyPI publish - if it's missing.**

## 3. Commit, push, tag, push the tag

```bash
git add mdinterface/__init__.py CHANGELOG.md
git commit -m "Bump version to X.Y.Z"
git push
git tag vX.Y.Z
git push --tags
```

In that order. Pushing the tag triggers `.github/workflows/release.yml`: `preflight` → `test` → `build` → `publish` (PyPI, automatic) → `release` (draft GitHub release, notes = changelog summary + link to CHANGELOG.md).

## 4. Publish the draft GitHub release

The GitHub release is created as a **draft** - go to the repo's Releases page and publish it by hand once you've skimmed the auto-filled notes. PyPI publishing already happened in the `publish` job and does not wait for this.

## Re-running a failed release

`workflow_dispatch` accepts a `tag` input, for re-running `release.yml` against an existing tag without re-tagging (e.g. after a transient PyPI outage). It does not re-verify a fresh `git push` - the tag must already exist and match `__version__`.
