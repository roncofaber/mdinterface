#!/usr/bin/env python3
import sys
import tarfile
import zipfile
from pathlib import Path


def verify_wheel(path: Path) -> None:
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()

    unexpected = [
        name for name in names
        if not name.startswith("mdinterface/") and ".dist-info/" not in name
    ]
    if unexpected:
        raise ValueError(f"Unexpected wheel contents: {unexpected}")
    if "mdinterface/config.ini" not in names:
        raise ValueError("Wheel is missing mdinterface/config.ini")


def verify_sdist(path: Path) -> None:
    with tarfile.open(path, "r:gz") as archive:
        names = archive.getnames()

    leaked = [name for name in names if "/docs/superpowers/" in name]
    if leaked:
        raise ValueError(f"Source distribution contains local planning files: {leaked}")


def main(argv) -> int:
    if len(argv) != 2:
        print("Usage: verify_distribution_contents.py DIST_DIR", file=sys.stderr)
        return 1

    dist_dir = Path(argv[1])
    wheels = list(dist_dir.glob("*.whl"))
    sdists = list(dist_dir.glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        print("ERROR: expected exactly one wheel and one source distribution", file=sys.stderr)
        return 1

    try:
        verify_wheel(wheels[0])
        verify_sdist(sdists[0])
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Verified {wheels[0].name} and {sdists[0].name}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
