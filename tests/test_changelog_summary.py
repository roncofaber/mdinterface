import importlib.util
from pathlib import Path

import pytest

_SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "extract_changelog_summary.py"
_spec = importlib.util.spec_from_file_location("extract_changelog_summary", _SCRIPT_PATH)
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
extract_summary = _module.extract_summary


def test_extract_summary_found():
    changelog = """# Changelog

## [Unreleased]

---

## [1.5.4] - 2026-08-11

Fixes a silent LAMMPS output corruption bug.

### Fixed
- Something
"""
    assert extract_summary(changelog, "1.5.4") == "Fixes a silent LAMMPS output corruption bug."


def test_extract_summary_multiline():
    changelog = """## [1.5.4] - 2026-08-11

Fixes a silent LAMMPS output corruption bug affecting systems that reuse
the same Specie template across multiple layers.

### Fixed
- Something
"""
    expected = (
        "Fixes a silent LAMMPS output corruption bug affecting systems that reuse\n"
        "the same Specie template across multiple layers."
    )
    assert extract_summary(changelog, "1.5.4") == expected


def test_missing_version_header_raises():
    changelog = """## [1.5.3] - 2026-07-09

Some summary.

### Fixed
- Something
"""
    with pytest.raises(ValueError, match="No changelog header found"):
        extract_summary(changelog, "1.5.4")


def test_missing_summary_paragraph_raises():
    changelog = """## [1.5.4] - 2026-08-11

### Fixed
- Something
"""
    with pytest.raises(ValueError, match="No summary paragraph found"):
        extract_summary(changelog, "1.5.4")


def test_second_paragraph_not_included():
    changelog = """## [1.5.4] - 2026-08-11

First paragraph only.

Second paragraph should not be included.

### Fixed
- Something
"""
    assert extract_summary(changelog, "1.5.4") == "First paragraph only."


def test_header_prefix_does_not_match_longer_version():
    # "## [1.5.4]" must not match when looking for "1.5" or "1.5.40"
    changelog = """## [1.5.40] - 2026-08-11

Wrong version's summary.

### Fixed
- Something
"""
    with pytest.raises(ValueError, match="No changelog header found"):
        extract_summary(changelog, "1.5.4")
