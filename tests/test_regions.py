#!/usr/bin/env python3
"""Unit tests for mdinterface.build.regions."""

import math
import pytest

from mdinterface.build.regions import Region, Sphere, Box, Cylinder, FilledRegion


class TestSphere:

    def test_volume(self):
        s = Sphere(center=(0, 0, 0), radius=2.0)
        assert s.volume() == pytest.approx((4.0 / 3.0) * math.pi * 2.0 ** 3)

    def test_packmol_line_inside(self):
        s = Sphere(center=(1.0, 2.0, 3.0), radius=4.0)
        assert s.packmol_line("inside") == "inside sphere 1.0 2.0 3.0 4.0"

    def test_packmol_line_outside(self):
        s = Sphere(center=(1.0, 2.0, 3.0), radius=4.0)
        assert s.packmol_line("outside") == "outside sphere 1.0 2.0 3.0 4.0"

    def test_packmol_line_invalid_mode_raises(self):
        s = Sphere(center=(0, 0, 0), radius=1.0)
        with pytest.raises(ValueError, match="mode must be"):
            s.packmol_line("sideways")

    def test_bounding_box(self):
        s = Sphere(center=(5.0, 5.0, 5.0), radius=2.0)
        assert s.bounding_box() == pytest.approx((3.0, 3.0, 3.0, 7.0, 7.0, 7.0))

    def test_nonpositive_radius_raises(self):
        with pytest.raises(ValueError, match="radius"):
            Sphere(center=(0, 0, 0), radius=0)


class TestBox:

    def test_volume(self):
        b = Box(center=(0, 0, 0), size=(2.0, 3.0, 4.0))
        assert b.volume() == pytest.approx(24.0)

    def test_bounding_box(self):
        b = Box(center=(10.0, 10.0, 10.0), size=(2.0, 4.0, 6.0))
        assert b.bounding_box() == pytest.approx((9.0, 8.0, 7.0, 11.0, 12.0, 13.0))

    def test_packmol_line_inside(self):
        b = Box(center=(10.0, 10.0, 10.0), size=(2.0, 2.0, 2.0))
        assert b.packmol_line("inside") == "inside box 9.0 9.0 9.0 11.0 11.0 11.0"

    def test_from_bounds(self):
        b = Box.from_bounds(0.0, 0.0, 0.0, 10.0, 20.0, 30.0)
        assert b.center == pytest.approx((5.0, 10.0, 15.0))
        assert b.size == pytest.approx((10.0, 20.0, 30.0))
        assert b.volume() == pytest.approx(6000.0)

    def test_nonpositive_size_raises(self):
        with pytest.raises(ValueError, match="size"):
            Box(center=(0, 0, 0), size=(1.0, 0.0, 1.0))


class TestCylinder:

    def test_volume(self):
        c = Cylinder(center=(0, 0, 0), radius=2.0, height=10.0)
        assert c.volume() == pytest.approx(math.pi * 2.0 ** 2 * 10.0)

    def test_packmol_line_z_axis(self):
        c = Cylinder(center=(5.0, 5.0, 10.0), radius=3.0, height=4.0, axis="z")
        assert c.packmol_line("inside") == "inside cylinder 5.0 5.0 8.0 0.0 0.0 1.0 3.0 4.0"

    def test_packmol_line_x_axis(self):
        c = Cylinder(center=(5.0, 5.0, 5.0), radius=1.0, height=2.0, axis="x")
        assert c.packmol_line("inside") == "inside cylinder 4.0 5.0 5.0 1.0 0.0 0.0 1.0 2.0"

    def test_bounding_box_z_axis(self):
        c = Cylinder(center=(0.0, 0.0, 0.0), radius=2.0, height=10.0, axis="z")
        assert c.bounding_box() == pytest.approx((-2.0, -2.0, -5.0, 2.0, 2.0, 5.0))

    def test_invalid_axis_raises(self):
        with pytest.raises(ValueError, match="axis"):
            Cylinder(center=(0, 0, 0), radius=1.0, height=1.0, axis="w")

    def test_nonpositive_radius_raises(self):
        with pytest.raises(ValueError, match="radius"):
            Cylinder(center=(0, 0, 0), radius=0, height=1.0)

    def test_nonpositive_height_raises(self):
        with pytest.raises(ValueError, match="height"):
            Cylinder(center=(0, 0, 0), radius=1.0, height=0)


class TestFill:

    def test_fill_returns_filled_region(self):
        s = Sphere(center=(0, 0, 0), radius=1.0)
        fr = s.fill(solvent="CO2", density=0.8)
        assert isinstance(fr, FilledRegion)
        assert fr.region is s
        assert fr.solvent == "CO2"
        assert fr.density == 0.8

    def test_fill_defaults(self):
        b = Box(center=(0, 0, 0), size=(1, 1, 1))
        fr = b.fill()
        assert fr.solvent is None
        assert fr.solute is None
        assert fr.nsolute is None
        assert fr.density is None
        assert fr.nsolvent is None
        assert fr.concentration is None
        assert fr.conmodel is None
        assert fr.ratio is None
