"""Unit tests for the LigParGen integration."""

from pathlib import Path
import subprocess

from ase import Atoms
import pytest

from mdinterface import LigParGenError
from mdinterface.externals.ligpargen import run_ligpargen


class TestRunLigParGenFailures:

    def test_missing_executable_raises_actionable_error(self, monkeypatch):
        monkeypatch.setattr("mdinterface.externals.ligpargen.shutil.which", lambda name: None)

        with pytest.raises(LigParGenError, match="executable was not found") as error:
            run_ligpargen(Atoms("NH3"))

        assert "python -m pip install" in str(error.value)
        assert "ligpargen -h" in str(error.value)
        assert error.value.tempdir is None

    def test_missing_bossdir_raises_actionable_error(self, monkeypatch):
        monkeypatch.setattr("mdinterface.externals.ligpargen.shutil.which", lambda name: f"/usr/bin/{name}")
        monkeypatch.setattr("mdinterface.config.load_config", lambda: None)
        monkeypatch.delenv("BOSSdir", raising=False)
        monkeypatch.setenv("MDINT_CONFIG_DIR", "/tmp/mdinterface-config")

        with pytest.raises(LigParGenError, match="BOSSdir is not configured") as error:
            run_ligpargen(Atoms("NH3"))

        assert "/tmp/mdinterface-config/config.ini" in str(error.value)

    def test_missing_openbabel_raises_actionable_error(self, monkeypatch):
        monkeypatch.setattr(
            "mdinterface.externals.ligpargen.shutil.which",
            lambda name: "/usr/bin/ligpargen" if name == "ligpargen" else None,
        )

        with pytest.raises(LigParGenError, match="Open Babel executable") as error:
            run_ligpargen(Atoms("NH3"))

        assert "conda install -c conda-forge openbabel" in str(error.value)
        assert "obabel -V" in str(error.value)

    def test_nonzero_exit_retains_diagnostics(self, monkeypatch, tmp_path):
        workdir = tmp_path / "ligpargen"
        workdir.mkdir()
        monkeypatch.setenv("BOSSdir", "boss-container:latest")
        monkeypatch.setattr("mdinterface.externals.ligpargen.shutil.which", lambda name: f"/usr/bin/{name}")
        monkeypatch.setattr("mdinterface.externals.ligpargen.tempfile.mkdtemp", lambda prefix: str(workdir))

        def unsuccessful(*args, **kwargs):
            raise subprocess.CalledProcessError(17, args[0], output="output", stderr="failure")

        monkeypatch.setattr("mdinterface.externals.ligpargen.subprocess.run", unsuccessful)

        with pytest.raises(LigParGenError, match="return code: 17") as error:
            run_ligpargen(Atoms("NH3"))

        assert error.value.tempdir == str(workdir)
        assert error.value.log_path == str(workdir / "ligpargen.log")
        assert Path(error.value.log_path).read_text() == "STDOUT:\noutput\nSTDERR:\nfailure\n"

    def test_missing_output_retains_diagnostics(self, monkeypatch, tmp_path):
        workdir = tmp_path / "ligpargen"
        workdir.mkdir()
        monkeypatch.setenv("BOSSdir", "boss-container:latest")
        monkeypatch.setattr("mdinterface.externals.ligpargen.shutil.which", lambda name: f"/usr/bin/{name}")
        monkeypatch.setattr("mdinterface.externals.ligpargen.tempfile.mkdtemp", lambda prefix: str(workdir))
        monkeypatch.setattr(
            "mdinterface.externals.ligpargen.subprocess.run",
            lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 0, stdout="output", stderr=""),
        )

        with pytest.raises(LigParGenError, match="did not create the expected output") as error:
            run_ligpargen(Atoms("NH3"))

        assert error.value.tempdir == str(workdir)
        assert "LigParGen output: output" in str(error.value)
        assert Path(error.value.log_path).exists()
