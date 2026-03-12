#!/usr/bin/env python3
"""Unit tests for mdinterface.build.solvent.populate_solutes."""

import warnings
import pytest
import numpy as np

from mdinterface.build.solvent import populate_solutes
from mdinterface.database import Ion


@pytest.fixture
def na():
    return Ion("Na")


@pytest.fixture
def cl():
    return Ion("Cl")


@pytest.fixture
def vol():
    return [20.0, 20.0, 30.0]


# ---------------------------------------------------------------------------
# PACKMOL paths (default, "packmol", "left", "right")
# ---------------------------------------------------------------------------

class TestPopulateSolutesPackmol:

    def test_default_returns_box_instruction(self, na, vol):
        instr = populate_solutes([na], 3, vol)
        assert len(instr) == 1
        _, count, mode, bounds = instr[0]
        assert mode == "box"
        assert count == 3

    def test_packmol_explicit(self, na, vol):
        instr = populate_solutes([na], 2, vol, solute_pos="packmol")
        _, count, mode, bounds = instr[0]
        assert mode == "box"
        assert count == 2

    def test_default_bounds_full_box(self, na, vol):
        instr = populate_solutes([na], 1, vol)
        _, _, _, bounds = instr[0]
        # zmax should be near vol[2] - 1
        assert bounds[5] == pytest.approx(vol[2] - 1.0)
        assert bounds[2] == pytest.approx(1.0)  # zmin

    def test_left_bounds_half_z(self, na, vol):
        instr = populate_solutes([na], 1, vol, solute_pos="left")
        _, _, _, bounds = instr[0]
        assert bounds[5] == pytest.approx(vol[2] / 2)
        assert bounds[2] == pytest.approx(1.0)

    def test_right_bounds_half_z(self, na, vol):
        instr = populate_solutes([na], 1, vol, solute_pos="right")
        _, _, _, bounds = instr[0]
        assert bounds[2] == pytest.approx(vol[2] / 2)
        assert bounds[5] == pytest.approx(vol[2] - 1.0)

    def test_multi_species_list_nsolute(self, na, cl, vol):
        instr = populate_solutes([na, cl], [3, 5], vol)
        assert len(instr) == 2
        assert instr[0][1] == 3
        assert instr[1][1] == 5

    def test_multi_species_int_nsolute(self, na, cl, vol):
        instr = populate_solutes([na, cl], 4, vol)
        assert instr[0][1] == 4
        assert instr[1][1] == 4

    def test_volume_not_mutated(self, na, vol):
        original = list(vol)
        populate_solutes([na], 2, vol, solute_pos="left")
        assert vol == original


# ---------------------------------------------------------------------------
# Center path (fixed coord)
# ---------------------------------------------------------------------------

class TestPopulateSolutesCenter:

    def test_center_returns_fixed_instructions(self, na, vol):
        instr = populate_solutes([na], 3, vol, solute_pos="center")
        assert len(instr) == 3
        for _, coord, mode in instr:
            assert mode == "fixed"

    def test_center_coord_is_midpoint(self, na, vol):
        instr = populate_solutes([na], 1, vol, solute_pos="center")
        _, coord, _ = instr[0]
        assert coord == pytest.approx([v / 2 for v in vol])


# ---------------------------------------------------------------------------
# Deprecated alias
# ---------------------------------------------------------------------------

class TestPopulateSolutesDeprecation:

    def test_box_alias_warns(self, na, vol):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            populate_solutes([na], 1, vol, solute_pos="box")
        assert any(issubclass(x.category, DeprecationWarning) for x in w)

    def test_box_alias_still_works(self, na, vol):
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            instr = populate_solutes([na], 2, vol, solute_pos="box")
        _, count, mode, _ = instr[0]
        assert mode == "box"
        assert count == 2


# ---------------------------------------------------------------------------
# Unknown placement
# ---------------------------------------------------------------------------

class TestPopulateSolutesUnknown:

    def test_unknown_raises(self, na, vol):
        with pytest.raises(ValueError, match="Unknown solute_pos"):
            populate_solutes([na], 1, vol, solute_pos="diagonal")
