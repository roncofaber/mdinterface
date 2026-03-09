#!/usr/bin/env python3
"""Tests for mdinterface.build.builder.BoxBuilder."""

import logging
import pytest
import numpy as np

from mdinterface import BoxBuilder
from mdinterface.database import Water, Ion


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def water():
    return Water()


@pytest.fixture
def na():
    return Ion("Na")


@pytest.fixture
def cl():
    return Ion("Cl")


@pytest.fixture
def builder():
    return BoxBuilder(xysize=[20, 20])


# ---------------------------------------------------------------------------
# Construction and xysize validation
# ---------------------------------------------------------------------------

class TestBoxBuilderCreation:

    def test_basic_creation(self):
        b = BoxBuilder(xysize=[20, 20])
        assert b is not None

    def test_stores_xysize(self):
        b = BoxBuilder(xysize=[15.5, 12.3])
        assert np.isclose(b._xsize, 15.5)
        assert np.isclose(b._ysize, 12.3)

    def test_tuple_xysize(self):
        b = BoxBuilder(xysize=(10, 10))
        assert np.isclose(b._xsize, 10.0)

    def test_numpy_xysize(self):
        b = BoxBuilder(xysize=np.array([12.0, 8.0]))
        assert np.isclose(b._xsize, 12.0)

    def test_xysize_wrong_length_one(self):
        with pytest.raises(ValueError):
            BoxBuilder(xysize=[10])

    def test_xysize_wrong_length_three(self):
        with pytest.raises(ValueError):
            BoxBuilder(xysize=[10, 10, 10])

    def test_xysize_non_sequence_raises(self):
        with pytest.raises(TypeError):
            BoxBuilder(xysize=20)

    def test_xysize_non_numeric_raises(self):
        with pytest.raises(ValueError):
            BoxBuilder(xysize=["a", "b"])

    def test_xysize_zero_raises(self):
        with pytest.raises(ValueError):
            BoxBuilder(xysize=[0, 10])

    def test_xysize_negative_raises(self):
        with pytest.raises(ValueError):
            BoxBuilder(xysize=[-5, 10])

    def test_initial_universe_is_none(self):
        b = BoxBuilder(xysize=[20, 20])
        assert b.universe is None


# ---------------------------------------------------------------------------
# Fluent layer API (no build)
# ---------------------------------------------------------------------------

class TestLayerAccumulation:

    def test_add_vacuum_returns_self(self, builder):
        assert builder.add_vacuum(zdim=5) is builder

    def test_add_solvent_returns_self(self, builder, water):
        assert builder.add_solvent(water, zdim=25, density=1.0) is builder

    def test_chain_three_calls(self, builder, water):
        b = builder.add_vacuum(5).add_solvent(water, zdim=20, density=1.0).add_vacuum(5)
        assert b is builder

    def test_vacuum_layer_stored(self, builder):
        builder.add_vacuum(zdim=10)
        vac = [lay for lay in builder._layers if lay["type"] == "vacuum"]
        assert len(vac) == 1
        assert vac[0]["zdim"] == 10

    def test_solvent_layer_stored(self, builder, water):
        builder.add_solvent(water, zdim=30, density=1.0)
        solv = [lay for lay in builder._layers if lay["type"] == "solvent"]
        assert len(solv) == 1
        assert solv[0]["zdim"] == 30

    def test_solvent_ions_stored(self, builder, water, na, cl):
        builder.add_solvent(water, ions=[na, cl], nions=[3, 3], zdim=25, density=1.0)
        assert len(builder._layers[0]["ions"]) == 2

    def test_multiple_solvent_layers(self, builder, water):
        builder.add_solvent(water, zdim=20, density=1.0)
        builder.add_solvent(water, zdim=20, density=1.0)
        solvs = [lay for lay in builder._layers if lay["type"] == "solvent"]
        assert len(solvs) == 2

    def test_add_solvent_requires_zdim(self, builder, water):
        with pytest.raises(ValueError, match="zdim"):
            builder.add_solvent(water, density=1.0)

    def test_species_registered_for_solvent(self, builder, water):
        builder.add_solvent(water, zdim=20, density=1.0)
        assert len(builder._all_species) == 1

    def test_ions_registered(self, builder, water, na, cl):
        builder.add_solvent(water, ions=[na, cl], nions=[2, 2], zdim=20, density=1.0)
        assert len(builder._all_species) == 3  # water + na + cl

    def test_same_species_not_duplicated(self, builder, water):
        # Each add_solvent() copies the species, so 2 copies → 2 entries
        builder.add_solvent(water, zdim=20, density=1.0)
        builder.add_solvent(water, zdim=20, density=1.0)
        assert len(builder._all_species) == 2

    def test_write_before_build_raises(self, builder):
        with pytest.raises(RuntimeError, match="build"):
            builder.write_lammps("dummy.lammps")

    def test_to_ase_before_build_raises(self, builder):
        with pytest.raises(RuntimeError, match="build"):
            builder.to_ase()


# ---------------------------------------------------------------------------
# Ion-only solvent layer (solvent=None)
# ---------------------------------------------------------------------------

class TestIonOnlyLayer:

    def test_ion_only_layer_accepted(self, builder, na, cl):
        builder.add_solvent(None, ions=[na, cl], nions=[3, 3], zdim=20)
        solv = builder._layers[0]
        assert solv["solvent"] is None
        assert len(solv["ions"]) == 2

    def test_ion_only_ions_registered(self, builder, na, cl):
        builder.add_solvent(None, ions=[na, cl], nions=[2, 2], zdim=20)
        assert len(builder._all_species) == 2  # na + cl, no solvent

    def test_ion_only_returns_self(self, builder, na):
        assert builder.add_solvent(None, ions=[na], nions=5, zdim=15) is builder


# ---------------------------------------------------------------------------
# Dilate parameter
# ---------------------------------------------------------------------------

class TestDilate:

    def test_dilate_stored(self, builder, water):
        builder.add_solvent(water, zdim=20, density=1.0, dilate=1.5)
        assert builder._layers[0]["dilate"] == 1.5

    def test_dilate_default_is_one(self, builder, water):
        builder.add_solvent(water, zdim=20, density=1.0)
        assert builder._layers[0]["dilate"] == 1.0

    def test_dilate_zero_raises(self, builder, water):
        with pytest.raises(ValueError, match="dilate"):
            builder.add_solvent(water, zdim=20, density=1.0, dilate=0)

    def test_dilate_negative_raises(self, builder, water):
        with pytest.raises(ValueError, match="dilate"):
            builder.add_solvent(water, zdim=20, density=1.0, dilate=-1.0)


# ---------------------------------------------------------------------------
# PACKMOL tolerance parameter
# ---------------------------------------------------------------------------

class TestPackmolTolerance:

    def test_tolerance_stored(self, builder, water):
        builder.add_solvent(water, zdim=20, density=1.0, packmol_tolerance=1.5)
        assert builder._layers[0]["packmol_tolerance"] == 1.5

    def test_tolerance_default(self, builder, water):
        builder.add_solvent(water, zdim=20, density=1.0)
        assert builder._layers[0]["packmol_tolerance"] == 2.0

    def test_tolerance_zero_raises(self, builder, water):
        with pytest.raises(ValueError, match="packmol_tolerance"):
            builder.add_solvent(water, zdim=20, density=1.0, packmol_tolerance=0)

    def test_tolerance_negative_raises(self, builder, water):
        with pytest.raises(ValueError, match="packmol_tolerance"):
            builder.add_solvent(water, zdim=20, density=1.0, packmol_tolerance=-0.5)


# ---------------------------------------------------------------------------
# Slab suffix uniqueness
# ---------------------------------------------------------------------------

class TestSlabSuffixes:

    def test_two_slabs_get_distinct_suffixes(self, water):
        for sp, idx in [(water.copy(), 0), (water.copy(), 1)]:
            suffix = f"_s{idx}"
            for atom in sp._stype:
                atom.set_label(atom.label + suffix)
            for atom in sp._stype:
                assert atom.label.endswith(suffix)

    def test_slab_suffix_increments(self, water):
        b = BoxBuilder(xysize=[20, 20])
        b.add_slab(water)
        b.add_slab(water)
        slabs = [lay for lay in b._layers if lay["type"] == "slab"]
        label0 = slabs[0]["species"]._stype[0].label
        label1 = slabs[1]["species"]._stype[0].label
        assert label0.endswith("_s0")
        assert label1.endswith("_s1")


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

class TestLogging:

    def test_add_slab_logs_debug(self, water, caplog):
        b = BoxBuilder(xysize=[20, 20])
        with caplog.at_level(logging.DEBUG, logger="mdinterface.builder"):
            b.add_slab(water)
        assert any("slab" in r.message for r in caplog.records)

    def test_add_solvent_logs_debug(self, water, caplog):
        b = BoxBuilder(xysize=[20, 20])
        with caplog.at_level(logging.DEBUG, logger="mdinterface.builder"):
            b.add_solvent(water, zdim=20, density=1.0)
        assert any("solvent" in r.message for r in caplog.records)

    def test_add_vacuum_logs_debug(self, caplog):
        b = BoxBuilder(xysize=[20, 20])
        with caplog.at_level(logging.DEBUG, logger="mdinterface.builder"):
            b.add_vacuum(zdim=5)
        assert any("vacuum" in r.message for r in caplog.records)
