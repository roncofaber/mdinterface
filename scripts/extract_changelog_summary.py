#!/usr/bin/env python3
"""Extract a version's summary paragraph from CHANGELOG.md.

Usage: python scripts/extract_changelog_summary.py VERSION

Prints the summary paragraph to stdout and exits 0 on success. Exits 1 with
an error on stderr if the version's changelog header, or its summary
paragraph, is missing.
"""
import sys
from pathlib import Path


def extract_summary(changelog_text: str, version: str) -> str:
    """Return the summary paragraph for `## [version]` in changelog_text.

    The summary is every non-blank line immediately following the header
    (after skipping any blank lines right after it), stopping at the next
    blank line or the next '## '/'### ' header — i.e. exactly one paragraph.
    """
    header = f"## [{version}]"
    lines = changelog_text.splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith(header):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(
            f"No changelog header found for version {version!r} "
            f"(expected a line starting with {header!r})"
        )

    i = header_idx + 1
    while i < len(lines) and lines[i].strip() == "":
        i += 1

    summary_lines = []
    while i < len(lines):
        line = lines[i]
        if line.strip() == "" or line.startswith("## ") or line.startswith("### "):
            break
        summary_lines.append(line)
        i += 1

    summary = "\n".join(summary_lines).strip()
    if not summary:
        raise ValueError(
            f"No summary paragraph found under {header!r} in CHANGELOG.md - "
            "add one before the first '### ' group"
        )

    return summary


def main(argv):
    if len(argv) != 2:
        print("Usage: extract_changelog_summary.py VERSION", file=sys.stderr)
        return 1

    version = argv[1]
    changelog_path = Path(__file__).resolve().parent.parent / "CHANGELOG.md"
    changelog_text = changelog_path.read_text()

    try:
        summary = extract_summary(changelog_text, version)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
