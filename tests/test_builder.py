#!/usr/bin/env python3
"""Tests for mdinterface.build.builder.SimCell."""

import logging
import pytest
import numpy as np

from mdinterface import SimCell
from mdinterface.utils.logger import set_verbosity
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
    return SimCell(xysize=[20, 20])


# ---------------------------------------------------------------------------
# Construction and xysize validation
# ---------------------------------------------------------------------------

class TestSimCellCreation:

    def test_basic_creation(self):
        b = SimCell(xysize=[20, 20])
        assert b is not None

    def test_stores_xysize(self):
        b = SimCell(xysize=[15.5, 12.3])
        assert np.isclose(b._xsize, 15.5)
        assert np.isclose(b._ysize, 12.3)

    def test_tuple_xysize(self):
        b = SimCell(xysize=(10, 10))
        assert np.isclose(b._xsize, 10.0)

    def test_numpy_xysize(self):
        b = SimCell(xysize=np.array([12.0, 8.0]))
        assert np.isclose(b._xsize, 12.0)

    def test_xysize_wrong_length_one(self):
        with pytest.raises(ValueError):
            SimCell(xysize=[10])

    def test_xysize_wrong_length_three(self):
        with pytest.raises(ValueError):
            SimCell(xysize=[10, 10, 10])

    def test_xysize_non_sequence_raises(self):
        with pytest.raises(TypeError):
            SimCell(xysize=20)

    def test_xysize_non_numeric_raises(self):
        with pytest.raises(ValueError):
            SimCell(xysize=["a", "b"])

    def test_xysize_zero_raises(self):
        with pytest.raises(ValueError):
            SimCell(xysize=[0, 10])

    def test_xysize_negative_raises(self):
        with pytest.raises(ValueError):
            SimCell(xysize=[-5, 10])

    def test_initial_universe_is_none(self):
        b = SimCell(xysize=[20, 20])
        assert b.universe is None


# ---------------------------------------------------------------------------
# Fluent layer API (no build)
# ---------------------------------------------------------------------------

class TestLayerAccumulation:

    def test_add_vacuum_returns_none(self, builder):
        assert builder.add_vacuum(zdim=5) is None

    def test_add_solvent_returns_none(self, builder, water):
        assert builder.add_solvent(water, zdim=25, density=1.0) is None

    def test_chain_three_calls(self, builder, water):
        builder.add_vacuum(5)
        builder.add_solvent(water, zdim=20, density=1.0)
        builder.add_vacuum(5)
        assert len(builder._layers) == 3

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

    def test_solvent_solute_stored(self, builder, water, na, cl):
        builder.add_solvent(water, solute=[na, cl], nsolute=[3, 3], zdim=25, density=1.0)
        assert len(builder._layers[0]["solute"]) == 2

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

    def test_solute_registered(self, builder, water, na, cl):
        builder.add_solvent(water, solute=[na, cl], nsolute=[2, 2], zdim=20, density=1.0)
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
# Solute-only solvent layer (solvent=None)
# ---------------------------------------------------------------------------

class TestSoluteOnlyLayer:

    def test_solute_only_layer_accepted(self, builder, na, cl):
        builder.add_solvent(None, solute=[na, cl], nsolute=[3, 3], zdim=20)
        solv = builder._layers[0]
        assert solv["solvent"] == []
        assert len(solv["solute"]) == 2

    def test_solute_only_registered(self, builder, na, cl):
        builder.add_solvent(None, solute=[na, cl], nsolute=[2, 2], zdim=20)
        assert len(builder._all_species) == 2  # na + cl, no solvent

    def test_solute_only_returns_none(self, builder, na):
        assert builder.add_solvent(None, solute=[na], nsolute=5, zdim=15) is None


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
# Multi-solvent support
# ---------------------------------------------------------------------------

class TestMultiSolvent:

    def test_solvent_list_stored(self, builder, water, na):
        """A list of solvents is stored as a list, not wrapped in another list."""
        builder.add_solvent([water, na], zdim=20, nsolvent=[10, 5])
        solv = builder._layers[0]
        assert isinstance(solv["solvent"], list)
        assert len(solv["solvent"]) == 2

    def test_single_solvent_stored_as_list(self, builder, water):
        """A single Specie is normalised to a 1-element list."""
        builder.add_solvent(water, zdim=20, density=1.0)
        solv = builder._layers[0]
        assert isinstance(solv["solvent"], list)
        assert len(solv["solvent"]) == 1

    def test_ratio_stored(self, builder, water, na):
        builder.add_solvent([water, na], zdim=20, ratio=[3, 1], density=1.0)
        assert builder._layers[0]["ratio"] == [3, 1]

    def test_ratio_none_by_default(self, builder, water):
        builder.add_solvent(water, zdim=20, density=1.0)
        assert builder._layers[0]["ratio"] is None

    def test_nsolvent_list_stored(self, builder, water, na):
        builder.add_solvent([water, na], zdim=20, nsolvent=[50, 10])
        assert builder._layers[0]["nsolvent"] == [50, 10]

    def test_multi_solvent_registers_all(self, builder, water, na, cl):
        builder.add_solvent([water, na, cl], zdim=20, nsolvent=[100, 5, 5])
        assert len(builder._all_species) == 3


# ---------------------------------------------------------------------------
# match_cell resolution
# ---------------------------------------------------------------------------

class TestMatchCell:

    def test_false_gives_no_ref_no_match(self, builder):
        cell_ref, do_match = builder._resolve_match_cell(False)
        assert cell_ref is None
        assert do_match is False

    def test_true_gives_no_ref_with_match(self, builder):
        cell_ref, do_match = builder._resolve_match_cell(True)
        assert cell_ref is None
        assert do_match is True

    def test_species_gives_ref_and_match(self, builder, water):
        cell_ref, do_match = builder._resolve_match_cell(water)
        assert cell_ref is water
        assert do_match is True

    def test_species_ref_is_exact_object(self, builder, na, cl):
        cell_ref, _ = builder._resolve_match_cell(na)
        assert cell_ref is na
        assert cell_ref is not cl


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
        b = SimCell(xysize=[20, 20])
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
        b = SimCell(xysize=[20, 20])
        with caplog.at_level(logging.DEBUG, logger="mdinterface"):
            b.add_slab(water)
        assert any("slab" in r.message for r in caplog.records)

    def test_add_solvent_logs_debug(self, water, caplog):
        b = SimCell(xysize=[20, 20])
        with caplog.at_level(logging.DEBUG, logger="mdinterface"):
            b.add_solvent(water, zdim=20, density=1.0)
        assert any("solvent" in r.message for r in caplog.records)

    def test_add_vacuum_logs_debug(self, caplog):
        b = SimCell(xysize=[20, 20])
        with caplog.at_level(logging.DEBUG, logger="mdinterface"):
            b.add_vacuum(zdim=5)
        assert any("vacuum" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# verbose parameter and _configure_logger
# ---------------------------------------------------------------------------

class TestVerbose:

    def _clean_handlers(self):
        """Remove StreamHandlers added by previous tests from the parent logger."""
        lg = logging.getLogger("mdinterface")
        lg.handlers = [h for h in lg.handlers
                       if not isinstance(h, logging.StreamHandler)]

    def test_verbose_true_sets_info(self):
        self._clean_handlers()
        set_verbosity(True)
        assert logging.getLogger("mdinterface").level == logging.INFO

    def test_verbose_false_sets_warning(self):
        self._clean_handlers()
        set_verbosity(False)
        assert logging.getLogger("mdinterface").level == logging.WARNING

    def test_verbose_string_debug(self):
        self._clean_handlers()
        set_verbosity("DEBUG")
        assert logging.getLogger("mdinterface").level == logging.DEBUG

    def test_verbose_int_level(self):
        self._clean_handlers()
        set_verbosity(logging.WARNING)
        assert logging.getLogger("mdinterface").level == logging.WARNING

    def test_verbose_in_constructor(self):
        self._clean_handlers()
        b = SimCell(xysize=[20, 20], verbose="DEBUG")
        assert logging.getLogger("mdinterface").level == logging.DEBUG
        assert b is not None

    def test_verbose_none_does_not_add_handler(self):
        self._clean_handlers()
        n_before = len(logging.getLogger("mdinterface").handlers)
        SimCell(xysize=[20, 20], verbose=None)
        n_after = len(logging.getLogger("mdinterface").handlers)
        assert n_after == n_before

    def test_configure_logger_adds_handler(self):
        self._clean_handlers()
        set_verbosity(True)
        handlers = logging.getLogger("mdinterface").handlers
        assert any(isinstance(h, logging.StreamHandler) for h in handlers)

    def test_configure_logger_not_duplicated(self):
        self._clean_handlers()
        set_verbosity(True)
        set_verbosity(True)
        handlers = [h for h in logging.getLogger("mdinterface").handlers
                    if isinstance(h, logging.StreamHandler)]
        assert len(handlers) == 1
