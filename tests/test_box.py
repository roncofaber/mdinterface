#!/usr/bin/env python3
"""Unit tests for mdinterface.build.box.populate_box, focused on the region-aware constraint paths."""

import pytest

from mdinterface.build.box import populate_box
from mdinterface.build.regions import Sphere
from mdinterface.database import Ion


@pytest.fixture
def na():
    return Ion("Na")


class TestPopulateBoxRegion:

    def test_region_instruction_packs_inside_sphere(self, na):
        sphere = Sphere(center=(10.0, 10.0, 10.0), radius=4.0)
        universe = populate_box(
            [20.0, 20.0, 20.0],
            [(na.to_universe(), 5, "region", sphere)],
        )
        assert universe is not None
        assert len(universe.atoms) == 5
        for pos in universe.atoms.positions:
            dist = ((pos[0] - 10.0) ** 2 + (pos[1] - 10.0) ** 2 + (pos[2] - 10.0) ** 2) ** 0.5
            assert dist <= 4.1

    def test_box_with_extra_outside_constraint(self, na):
        # Pack into the full box but exclude a sphere in the middle - no atom
        # should land inside the excluded sphere.
        sphere = Sphere(center=(10.0, 10.0, 10.0), radius=5.0)
        universe = populate_box(
            [20.0, 20.0, 20.0],
            [(na.to_universe(), 8, "box", None, [sphere.packmol_line("outside")])],
        )
        assert universe is not None
        assert len(universe.atoms) == 8
        for pos in universe.atoms.positions:
            dist = ((pos[0] - 10.0) ** 2 + (pos[1] - 10.0) ** 2 + (pos[2] - 10.0) ** 2) ** 0.5
            assert dist >= 4.9

    def test_box_without_extra_constraints_unchanged(self, na):
        # 3-tuple form (no bounds, no extra constraints) still works exactly as before.
        universe = populate_box([20.0, 20.0, 20.0], [(na.to_universe(), 3, "box")])
        assert universe is not None
        assert len(universe.atoms) == 3

    def test_unknown_instruction_type_raises(self, na):
        with pytest.raises(ValueError, match="Wrong instructions"):
            populate_box([20.0, 20.0, 20.0], [(na.to_universe(), 1, "diagonal")])
