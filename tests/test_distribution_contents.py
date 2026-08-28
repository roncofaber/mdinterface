import io
import tarfile
import zipfile

import pytest

from scripts.verify_distribution_contents import verify_sdist, verify_wheel


def write_wheel(path, names):
    with zipfile.ZipFile(path, "w") as archive:
        for name in names:
            archive.writestr(name, "")


def write_sdist(path, names):
    with tarfile.open(path, "w:gz") as archive:
        for name in names:
            info = tarfile.TarInfo(name)
            archive.addfile(info, io.BytesIO())


def test_verify_wheel_accepts_package_and_metadata(tmp_path):
    path = tmp_path / "package.whl"
    write_wheel(path, ["mdinterface/__init__.py", "mdinterface/config.ini", "mdinterface-1.0.dist-info/METADATA"])
    verify_wheel(path)


def test_verify_wheel_rejects_unexpected_top_level_content(tmp_path):
    path = tmp_path / "package.whl"
    write_wheel(path, ["mdinterface/config.ini", "tests/test_package.py"])
    with pytest.raises(ValueError, match="Unexpected wheel contents"):
        verify_wheel(path)


def test_verify_wheel_requires_default_config(tmp_path):
    path = tmp_path / "package.whl"
    write_wheel(path, ["mdinterface/__init__.py"])
    with pytest.raises(ValueError, match="missing mdinterface/config.ini"):
        verify_wheel(path)


def test_verify_sdist_rejects_local_planning_files(tmp_path):
    path = tmp_path / "package.tar.gz"
    write_sdist(path, ["package/docs/superpowers/plans/private.md"])
    with pytest.raises(ValueError, match="local planning files"):
        verify_sdist(path)
