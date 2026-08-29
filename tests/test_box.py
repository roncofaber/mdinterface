#!/usr/bin/env python3
"""Unit tests for mdinterface.build.box."""

from pathlib import Path
import subprocess

from ase import Atoms
import numpy as np
import pytest

from mdinterface import PackmolError, Specie
from mdinterface.build.box import make_interface_slab, populate_box
from mdinterface.build.regions import Sphere
from mdinterface.database import Ion


@pytest.fixture
def na():
    return Ion("Na")


class TestPopulateBoxRegion:

    def test_region_instruction_packs_inside_sphere(self, na):
        sphere = Sphere(center=(10.0, 10.0, 10.0), radius=4.0)
        packed = populate_box(
            [20.0, 20.0, 20.0],
            [(na, 5, "region", sphere)],
        )
        assert packed is not None
        assert len(packed) == 5
        for pos in packed.positions:
            dist = ((pos[0] - 10.0) ** 2 + (pos[1] - 10.0) ** 2 + (pos[2] - 10.0) ** 2) ** 0.5
            assert dist <= 4.1

    def test_box_with_extra_outside_constraint(self, na):
        # Pack into the full box but exclude a sphere in the middle - no atom
        # should land inside the excluded sphere.
        sphere = Sphere(center=(10.0, 10.0, 10.0), radius=5.0)
        packed = populate_box(
            [20.0, 20.0, 20.0],
            [(na, 8, "box", None, [sphere.packmol_line("outside")])],
        )
        assert packed is not None
        assert len(packed) == 8
        for pos in packed.positions:
            dist = ((pos[0] - 10.0) ** 2 + (pos[1] - 10.0) ** 2 + (pos[2] - 10.0) ** 2) ** 0.5
            assert dist >= 4.9

    def test_box_without_extra_constraints_unchanged(self, na):
        # 3-tuple form (no bounds, no extra constraints) still works exactly as before.
        packed = populate_box([20.0, 20.0, 20.0], [(na, 3, "box")])
        assert packed is not None
        assert len(packed) == 3

    def test_unknown_instruction_type_raises(self, na):
        with pytest.raises(ValueError, match="Wrong instructions"):
            populate_box([20.0, 20.0, 20.0], [(na, 1, "diagonal")])


class TestPopulateBoxFailures:

    @staticmethod
    def use_tempdir(monkeypatch, tmp_path):
        workdir = tmp_path / "packmol"
        workdir.mkdir()
        monkeypatch.setattr(
            "mdinterface.build.box.tempfile.mkdtemp",
            lambda prefix: str(workdir),
        )
        return workdir

    def test_missing_executable_raises_packmol_error(self, na, monkeypatch, tmp_path):
        workdir = self.use_tempdir(monkeypatch, tmp_path)

        def executable_missing(*args, **kwargs):
            raise FileNotFoundError("packmol")

        monkeypatch.setattr("mdinterface.build.box.subprocess.run", executable_missing)

        with pytest.raises(PackmolError, match="executable was not found") as error:
            populate_box([20.0, 20.0, 20.0], [(na, 1, "box")])

        assert error.value.tempdir == str(workdir)
        assert error.value.log_path == str(workdir / "packmol.log")
        assert error.value.returncode is None
        assert workdir.exists()

    def test_nonzero_exit_raises_packmol_error(self, na, monkeypatch, tmp_path):
        workdir = self.use_tempdir(monkeypatch, tmp_path)
        monkeypatch.setattr(
            "mdinterface.build.box.subprocess.run",
            lambda *args, **kwargs: subprocess.CompletedProcess(["packmol"], 17),
        )

        with pytest.raises(PackmolError, match="return code: 17") as error:
            populate_box([20.0, 20.0, 20.0], [(na, 1, "box")])

        assert error.value.returncode == 17
        assert workdir.exists()

    def test_missing_output_raises_packmol_error(self, na, monkeypatch, tmp_path):
        workdir = self.use_tempdir(monkeypatch, tmp_path)
        monkeypatch.setattr(
            "mdinterface.build.box.subprocess.run",
            lambda *args, **kwargs: subprocess.CompletedProcess(["packmol"], 0),
        )

        with pytest.raises(PackmolError, match="did not create the expected output"):
            populate_box([20.0, 20.0, 20.0], [(na, 1, "box")])

        assert workdir.exists()

    def test_unreadable_output_raises_packmol_error(self, na, monkeypatch, tmp_path):
        workdir = self.use_tempdir(monkeypatch, tmp_path)

        def create_output(*args, **kwargs):
            Path(kwargs["cwd"], "system.pdb").write_text("invalid")
            return subprocess.CompletedProcess(["packmol"], 0)

        monkeypatch.setattr("mdinterface.build.box.subprocess.run", create_output)
        monkeypatch.setattr(
            "mdinterface.build.box.ase.io.read",
            lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("invalid PDB")),
        )

        with pytest.raises(PackmolError, match="output could not be read"):
            populate_box([20.0, 20.0, 20.0], [(na, 1, "box")])

        assert workdir.exists()

    def test_success_removes_temporary_files(self, na, monkeypatch, tmp_path):
        workdir = self.use_tempdir(monkeypatch, tmp_path)

        def create_output(*args, **kwargs):
            Path(kwargs["cwd"], "system.pdb").write_text("output")
            return subprocess.CompletedProcess(["packmol"], 0)

        monkeypatch.setattr("mdinterface.build.box.subprocess.run", create_output)
        monkeypatch.setattr("mdinterface.build.box.ase.io.read", lambda *args, **kwargs: Atoms("Na"))

        packed = populate_box([20.0, 20.0, 20.0], [(na, 1, "box")])

        assert len(packed) == 1
        assert not workdir.exists()

    @pytest.mark.parametrize(
        ("volume", "error", "message"),
        [
            ("20,20,20", TypeError, "list, tuple, or numpy array"),
            ([20.0, 20.0], ValueError, "exactly three"),
            ([20.0, "wide", 20.0], TypeError, "must be numeric"),
            ([20.0, 0.0, 20.0], ValueError, "finite and positive"),
            ([20.0, np.nan, 20.0], ValueError, "finite and positive"),
        ],
    )
    def test_invalid_volume_raises_before_creating_tempdir(self, na, monkeypatch, volume, error, message):
        monkeypatch.setattr(
            "mdinterface.build.box.tempfile.mkdtemp",
            lambda *args, **kwargs: pytest.fail("temporary directory should not be created"),
        )

        with pytest.raises(error, match=message):
            populate_box(volume, [(na, 1, "box")])


class TestMakeInterfaceSlab:

    @staticmethod
    def specie_with_cell(cell):
        return Specie(Atoms("He", positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=True))

    @pytest.mark.parametrize(
        ("target", "expected_cell", "expected_atoms"),
        [
            ((2.0, 1.0), (4.0, 3.0), 1),
            ((8.0, 6.0), (8.0, 6.0), 4),
            ((8.1, 6.1), (12.0, 9.0), 9),
        ],
    )
    def test_tiling_covers_target(self, target, expected_cell, expected_atoms):
        slab = make_interface_slab(
            self.specie_with_cell([4.0, 3.0, 5.0]),
            *target,
        )

        assert np.allclose(slab.atoms.cell.lengths()[:2], expected_cell)
        assert len(slab.atoms) == expected_atoms
        assert slab.atoms.cell.lengths()[0] >= target[0]
        assert slab.atoms.cell.lengths()[1] >= target[1]

    def test_nonorthogonal_cell_is_tiled_and_orthogonalized(self, caplog):
        with caplog.at_level("WARNING", logger="mdinterface.build.box"):
            slab = make_interface_slab(
                self.specie_with_cell([[4.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.0, 0.0, 5.0]]),
                8.1,
                6.1,
            )

        assert np.allclose(slab.atoms.cell.lengths(), [12.0, 9.0, 5.0])
        assert np.allclose(slab.atoms.cell.angles(), [90.0, 90.0, 90.0])
        assert len(slab.atoms) == 9
        assert "converted to an orthogonal cell" in caplog.text

    @pytest.mark.parametrize("target", [(0.0, 5.0), (5.0, -1.0)])
    def test_nonpositive_target_raises(self, target):
        with pytest.raises(ValueError, match="target dimensions must be positive"):
            make_interface_slab(self.specie_with_cell([4.0, 3.0, 5.0]), *target)

    def test_nonpositive_cell_projection_raises(self):
        with pytest.raises(ValueError, match="positive X and Y"):
            make_interface_slab(
                self.specie_with_cell([[0.0, 4.0, 0.0], [3.0, 0.0, 0.0], [0.0, 0.0, 5.0]]),
                5.0,
                5.0,
            )
