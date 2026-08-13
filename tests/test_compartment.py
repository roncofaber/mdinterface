#!/usr/bin/env python3
"""Unit tests for mdinterface.build.compartment."""

import pytest

from mdinterface.build.compartment import (
    Compartment, SlabCompartment, SolventCompartment, VacuumCompartment,
)
from mdinterface.database import Water


@pytest.fixture
def water():
    return Water()


class TestVacuumCompartment:

    def test_is_a_compartment(self):
        assert isinstance(VacuumCompartment(zdim=5.0), Compartment)

    def test_default_zdim_is_zero(self):
        assert VacuumCompartment().zdim == 0.0

    def test_build_returns_none(self):
        vac = VacuumCompartment(zdim=5.0)
        assert vac.build(20.0, 20.0, []) is None


class TestSlabCompartment:

    def test_is_a_compartment(self):
        slab = SlabCompartment(species=Water(), nlayers=2)
        assert isinstance(slab, Compartment)

    def test_default_label_is_slab(self):
        slab = SlabCompartment(species=Water(), nlayers=1)
        assert slab.label == "slab"

    def test_explicit_label_stored(self):
        slab = SlabCompartment(species=Water(), nlayers=1, label="prebuilt")
        assert slab.label == "prebuilt"

    def test_build_returns_none_before_slab_is_tiled(self, water):
        """_slab is only populated by SimCell._fit_slabs (Task 3); before that, build() is a no-op."""
        slab = SlabCompartment(species=water, nlayers=1)
        assert slab.build(20.0, 20.0, []) is None

    def test_native_xy_defaults_to_none(self, water):
        slab = SlabCompartment(species=water, nlayers=1)
        assert slab._native_xy is None
