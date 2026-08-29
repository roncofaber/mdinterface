#!/usr/bin/env python3
"""Unit tests for mdinterface.build.solvent.populate_solutes."""

import warnings
import pytest
import numpy as np

from mdinterface.build.solvent import populate_solutes
from mdinterface.database import Ion, Water


@pytest.fixture
def na():
    return Ion("Na")


@pytest.fixture
def cl():
    return Ion("Cl")


@pytest.fixture
def vol():
    return [20.0, 20.0, 30.0]


@pytest.fixture
def water():
    return Water()


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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            instr = populate_solutes([na], 1, vol, solute_pos="left")
        _, _, mode, region = instr[0]
        assert mode == "region"
        xmin, ymin, zmin, xmax, ymax, zmax = region.bounding_box()
        assert zmax == pytest.approx(vol[2] / 2)
        assert zmin == pytest.approx(1.0)

    def test_right_bounds_half_z(self, na, vol):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            instr = populate_solutes([na], 1, vol, solute_pos="right")
        _, _, mode, region = instr[0]
        assert mode == "region"
        xmin, ymin, zmin, xmax, ymax, zmax = region.bounding_box()
        assert zmin == pytest.approx(vol[2] / 2)
        assert zmax == pytest.approx(vol[2] - 1.0)

    def test_left_warns_deprecated(self, na, vol):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            populate_solutes([na], 1, vol, solute_pos="left")
        assert any(issubclass(x.category, DeprecationWarning) for x in w)

    def test_right_warns_deprecated(self, na, vol):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            populate_solutes([na], 1, vol, solute_pos="right")
        assert any(issubclass(x.category, DeprecationWarning) for x in w)

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


# ---------------------------------------------------------------------------
# Region-based placement
# ---------------------------------------------------------------------------

from mdinterface.build.regions import Sphere, Box


class TestPopulateSolutesRegion:

    def test_region_instance_returns_region_instructions(self, na, vol):
        sphere = Sphere(center=(10.0, 10.0, 15.0), radius=5.0)
        instr = populate_solutes([na], 2, vol, solute_pos=sphere)
        assert len(instr) == 1
        _, count, mode, region = instr[0]
        assert mode == "region"
        assert region is sphere
        assert count == 2


class TestMakeSolventBoxRegions:

    @pytest.mark.integration
    def test_bulk_and_region_species_both_present(self, water, na):
        from mdinterface.build.solvent import make_solvent_box

        pocket = Sphere(center=(10.0, 10.0, 15.0), radius=4.0)
        universe = make_solvent_box(
            species=[water.to_universe(), na.to_universe()],
            solvent=water,
            solute=None,
            volume=[20.0, 20.0, 30.0],
            density=1.0,
            nsolute=None,
            concentration=None,
            conmodel=None,
            solute_pos=None,
            nsolvent=None,
            regions=[pocket.fill(na, density=1.5)],
        )
        assert universe is not None
        resnames = {res.resname for res in universe.residues}
        assert water.resname in resnames
        assert na.resname in resnames

    @pytest.mark.integration
    def test_region_content_confined_to_region(self, water, na):
        from mdinterface.build.solvent import make_solvent_box

        pocket = Sphere(center=(10.0, 10.0, 15.0), radius=4.0)
        universe = make_solvent_box(
            species=[water.to_universe(), na.to_universe()],
            solvent=water,
            solute=None,
            volume=[20.0, 20.0, 30.0],
            density=1.0,
            nsolute=None,
            concentration=None,
            conmodel=None,
            solute_pos=None,
            nsolvent=None,
            regions=[pocket.fill(na, nsolvent=3)],
        )
        na_atoms = universe.select_atoms(f"resname {na.resname}")
        assert len(na_atoms) == 3
        for pos in na_atoms.positions:
            dist = sum((pos[i] - pocket.center[i]) ** 2 for i in range(3)) ** 0.5
            assert dist <= pocket.radius + 0.2

    def test_region_conmodel_raises(self, water, na):
        from mdinterface.build.solvent import make_solvent_box

        pocket = Sphere(center=(10.0, 10.0, 15.0), radius=4.0)
        filled = pocket.fill(na, nsolvent=1)
        filled.conmodel = {0: ([0, 1], [1, 1])}
        with pytest.raises(ValueError, match="conmodel is not yet supported"):
            make_solvent_box(
                species=[water.to_universe(), na.to_universe()],
                solvent=water, solute=None, volume=[20.0, 20.0, 30.0], density=1.0,
                nsolute=None, concentration=None, conmodel=None, solute_pos=None,
                nsolvent=None, regions=[filled],
            )

    def test_region_out_of_bounds_raises(self, water, na):
        from mdinterface.build.solvent import make_solvent_box

        # radius 10 centered near the edge - extends outside the 20x20x30 box
        pocket = Sphere(center=(2.0, 2.0, 2.0), radius=10.0)
        with pytest.raises(ValueError, match="extends outside"):
            make_solvent_box(
                species=[water.to_universe(), na.to_universe()],
                solvent=water, solute=None, volume=[20.0, 20.0, 30.0], density=1.0,
                nsolute=None, concentration=None, conmodel=None, solute_pos=None,
                nsolvent=None, regions=[pocket.fill(na, nsolvent=1)],
            )

    @pytest.mark.integration
    def test_overlapping_regions_warn(self, water, na):
        from mdinterface.build.solvent import make_solvent_box

        r1 = Sphere(center=(8.0, 10.0, 15.0), radius=4.0)
        r2 = Sphere(center=(10.0, 10.0, 15.0), radius=4.0)
        with pytest.warns(UserWarning, match="overlapping"):
            make_solvent_box(
                species=[water.to_universe(), na.to_universe()],
                solvent=water, solute=None, volume=[20.0, 20.0, 30.0], density=1.0,
                nsolute=None, concentration=None, conmodel=None, solute_pos=None,
                nsolvent=None, regions=[r1.fill(na, nsolvent=1), r2.fill(na, nsolvent=1)],
            )

    @pytest.mark.integration
    def test_nested_region_species_present(self, water, na, cl):
        from mdinterface.build.solvent import make_solvent_box

        outer = Sphere(center=(10.0, 10.0, 15.0), radius=6.0)
        inner = Sphere(center=(10.0, 10.0, 15.0), radius=2.0)
        outer_fr = outer.fill(na, nsolvent=3, regions=[inner.fill(cl, nsolvent=2)])
        universe = make_solvent_box(
            species=[water.to_universe(), na.to_universe(), cl.to_universe()],
            solvent=water, solute=None, volume=[20.0, 20.0, 30.0], density=1.0,
            nsolute=None, concentration=None, conmodel=None, solute_pos=None,
            nsolvent=None, regions=[outer_fr],
        )
        assert universe is not None
        resnames = {res.resname for res in universe.residues}
        assert water.resname in resnames
        assert na.resname in resnames
        assert cl.resname in resnames

    @pytest.mark.integration
    def test_nested_region_confined_to_inner_bounds(self, water, na, cl):
        from mdinterface.build.solvent import make_solvent_box

        outer = Sphere(center=(10.0, 10.0, 15.0), radius=6.0)
        inner = Sphere(center=(10.0, 10.0, 15.0), radius=2.0)
        outer_fr = outer.fill(na, nsolvent=3, regions=[inner.fill(cl, nsolvent=2)])
        universe = make_solvent_box(
            species=[water.to_universe(), na.to_universe(), cl.to_universe()],
            solvent=water, solute=None, volume=[20.0, 20.0, 30.0], density=1.0,
            nsolute=None, concentration=None, conmodel=None, solute_pos=None,
            nsolvent=None, regions=[outer_fr],
        )
        cl_atoms = universe.select_atoms(f"resname {cl.resname}")
        assert len(cl_atoms) == 2
        for pos in cl_atoms.positions:
            dist = sum((pos[i] - inner.center[i]) ** 2 for i in range(3)) ** 0.5
            assert dist <= inner.radius + 0.2

    def test_nested_region_out_of_bounds_raises(self, water, na, cl):
        from mdinterface.build.solvent import make_solvent_box

        outer = Sphere(center=(10.0, 10.0, 15.0), radius=4.0)
        # inner radius 5 centered at the outer's own center - extends past the outer's own bbox
        inner = Sphere(center=(10.0, 10.0, 15.0), radius=5.0)
        outer_fr = outer.fill(na, nsolvent=1, regions=[inner.fill(cl, nsolvent=1)])
        with pytest.raises(ValueError, match="extends outside"):
            make_solvent_box(
                species=[water.to_universe(), na.to_universe(), cl.to_universe()],
                solvent=water, solute=None, volume=[20.0, 20.0, 30.0], density=1.0,
                nsolute=None, concentration=None, conmodel=None, solute_pos=None,
                nsolvent=None, regions=[outer_fr],
            )


class TestResolveRandomRegions:

    def test_places_within_bounds(self):
        from mdinterface.build.solvent import _resolve_random_regions

        region = Sphere(center="random", radius=3.0)
        fr = region.fill(density=1.0)
        rng = np.random.default_rng(0)
        resolved = _resolve_random_regions([fr], (0.0, 0.0, 0.0, 20.0, 20.0, 30.0), rng)

        assert len(resolved) == 1
        assert resolved[0].region.center != "random"
        assert fr.region.center == "random"  # original untouched
        xmin, ymin, zmin, xmax, ymax, zmax = resolved[0].region.bounding_box()
        assert xmin >= 0.0 and ymin >= 0.0 and zmin >= 0.0
        assert xmax <= 20.0 and ymax <= 20.0 and zmax <= 30.0

    def test_avoids_sibling_overlap(self):
        from mdinterface.build.solvent import _resolve_random_regions
        from mdinterface.build.solvent import _bboxes_overlap

        fixed = Sphere(center=(5.0, 5.0, 5.0), radius=3.0)
        random_region = Sphere(center="random", radius=3.0)
        frs = [fixed.fill(density=1.0), random_region.fill(density=1.0)]
        rng = np.random.default_rng(1)
        resolved = _resolve_random_regions(frs, (0.0, 0.0, 0.0, 20.0, 20.0, 30.0), rng)

        fixed_bbox = resolved[0].region.bounding_box()
        random_bbox = resolved[1].region.bounding_box()
        assert not _bboxes_overlap(fixed_bbox, random_bbox)

    def test_raises_when_no_space(self):
        from mdinterface.build.solvent import _resolve_random_regions

        region = Sphere(center="random", radius=15.0)
        fr = region.fill(density=1.0)
        rng = np.random.default_rng(0)
        with pytest.raises(ValueError, match="too large"):
            _resolve_random_regions([fr], (0.0, 0.0, 0.0, 20.0, 20.0, 30.0), rng)

    def test_seed_reproducible(self):
        from mdinterface.build.solvent import _resolve_random_regions

        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(42)
        fr1 = Sphere(center="random", radius=3.0).fill(density=1.0)
        fr2 = Sphere(center="random", radius=3.0).fill(density=1.0)

        resolved1 = _resolve_random_regions([fr1], (0.0, 0.0, 0.0, 20.0, 20.0, 30.0), rng1)
        resolved2 = _resolve_random_regions([fr2], (0.0, 0.0, 0.0, 20.0, 20.0, 30.0), rng2)

        assert resolved1[0].region.center == resolved2[0].region.center

    def test_no_regions_returns_empty(self):
        from mdinterface.build.solvent import _resolve_random_regions

        rng = np.random.default_rng(0)
        assert _resolve_random_regions([], (0.0, 0.0, 0.0, 20.0, 20.0, 30.0), rng) == []


class TestMakeSolventBoxRandomRegion:

    @pytest.mark.integration
    def test_random_region_places_content_in_bounds(self, water, na):
        from mdinterface.build.solvent import make_solvent_box

        pocket = Sphere(center="random", radius=3.0)
        universe = make_solvent_box(
            species=[water.to_universe(), na.to_universe()],
            solvent=water, solute=None, volume=[20.0, 20.0, 30.0], density=1.0,
            nsolute=None, concentration=None, conmodel=None, solute_pos=None,
            nsolvent=None, regions=[pocket.fill(na, nsolvent=3)], seed=7,
        )
        na_atoms = universe.select_atoms(f"resname {na.resname}")
        assert len(na_atoms) == 3
        assert pocket.center == "random"  # original region untouched
        for pos in na_atoms.positions:
            assert 0.0 <= pos[0] <= 20.0
            assert 0.0 <= pos[1] <= 20.0
            assert 0.0 <= pos[2] <= 30.0
