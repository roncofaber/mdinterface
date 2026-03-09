#!/usr/bin/env python3
"""Tests for mdinterface.database species (Water, ions, noble gases)."""

import pytest
import numpy as np

from mdinterface.database.molecules import Water, Oxygen, Hydrogen
from mdinterface.database.nobles import (
    NobleGas, Argon, Neon, Krypton, Xenon,
    lookup_noble_gas_parameters,
)
from mdinterface.database.ions import Ion, Perchlorate, Hydronium, Hydroxide, lookup_parameters


class TestWater:
    def test_default_model_creates(self):
        w = Water()
        assert w is not None

    def test_has_three_atoms(self):
        w = Water()
        assert len(w.atoms) == 3

    def test_ewald_charges(self):
        w = Water("ewald")
        charges = w.charges
        assert charges[0] == pytest.approx(-0.83)
        assert charges[1] == pytest.approx(0.415)
        assert charges[2] == pytest.approx(0.415)

    def test_charge_sum_is_zero(self):
        for model in ("ewald", "charmm", "spce"):
            w = Water(model)
            assert sum(w.charges) == pytest.approx(0.0, abs=1e-6)

    def test_charmm_model(self):
        w = Water("charmm")
        assert len(w.atoms) == 3

    def test_spce_model(self):
        w = Water("spce")
        assert len(w.atoms) == 3

    def test_resname_is_four_chars(self):
        w = Water()
        assert len(w.resname) <= 4

    def test_bonds_present(self):
        w = Water()
        bonds = w.bonds
        assert bonds is not None
        assert bonds != [[], []]

    def test_angles_present(self):
        w = Water()
        angles = w.angles
        assert angles is not None
        assert angles != [[], []]

    def test_copy_is_independent(self):
        w1 = Water()
        w2 = w1.copy()
        w2.cutoff = 999.0
        assert w1.cutoff != 999.0


class TestOxygen:
    def test_creates(self):
        o = Oxygen()
        assert len(o.atoms) == 2

    def test_charge_sum_is_zero(self):
        o = Oxygen()
        assert sum(o.charges) == pytest.approx(0.0, abs=1e-6)


class TestHydrogen:
    def test_creates_std(self):
        h = Hydrogen()
        assert len(h.atoms) == 2

    def test_creates_alt(self):
        h = Hydrogen(Hset="alt")
        assert len(h.atoms) == 2


class TestNobleGasLookup:
    def test_argon_vrabec(self):
        charge, lj = lookup_noble_gas_parameters("Ar", "vrabec")
        assert charge == pytest.approx(0.0)
        assert lj[0] == pytest.approx(0.2321)  # epsilon
        assert lj[1] == pytest.approx(3.3952)  # sigma

    def test_invalid_element_raises(self):
        with pytest.raises(ValueError):
            lookup_noble_gas_parameters("Unobtanium", "vrabec")

    def test_invalid_ffield_raises(self):
        with pytest.raises(ValueError):
            lookup_noble_gas_parameters("Ar", "nonexistent_ff")


class TestNobleGas:
    @pytest.mark.parametrize("cls,symbol", [
        (Argon, "Ar"),
        (Neon, "Ne"),
        (Krypton, "Kr"),
        (Xenon, "Xe"),
    ])
    def test_creation(self, cls, symbol):
        ng = cls()
        assert len(ng.atoms) == 1
        assert ng.atoms.get_chemical_symbols()[0] == symbol

    def test_argon_charge_zero(self):
        ar = Argon()
        assert ar.charges[0] == pytest.approx(0.0)

    def test_noble_gas_generic(self):
        ng = NobleGas("Kr")
        assert ng.atoms.get_chemical_symbols()[0] == "Kr"

    def test_invalid_noble_gas_raises(self):
        with pytest.raises(ValueError):
            NobleGas("Unobtanium")


class TestIons:
    def test_sodium_creates(self):
        na = Ion("Na")
        assert len(na.atoms) == 1
        assert na.atoms.get_chemical_symbols()[0] == "Na"

    def test_sodium_charge(self):
        na = Ion("Na")
        assert na.charges[0] == pytest.approx(1.0)

    def test_chloride_creates(self):
        cl = Ion("Cl")
        assert len(cl.atoms) == 1
        assert cl.atoms.get_chemical_symbols()[0] == "Cl"

    def test_chloride_charge(self):
        cl = Ion("Cl")
        assert cl.charges[0] == pytest.approx(-1.0)

    def test_potassium_creates(self):
        k = Ion("K")
        assert len(k.atoms) == 1
        assert k.atoms.get_chemical_symbols()[0] == "K"

    def test_charge_scaling(self):
        na = Ion("Na", chg_scaling=0.8)
        assert na.charges[0] == pytest.approx(0.8)

    def test_invalid_element_raises(self):
        with pytest.raises(ValueError):
            Ion("Unobtanium")

    def test_lookup_parameters(self):
        charge, lj = lookup_parameters("Na", "jorgensen")
        assert charge == pytest.approx(1.0)
        assert len(lj) == 2

    def test_perchlorate_creates(self):
        p = Perchlorate()
        assert len(p.atoms) == 5  # ClO4-

    def test_hydronium_creates(self):
        h3o = Hydronium()
        assert len(h3o.atoms) == 4  # OH3+

    def test_hydroxide_creates(self):
        oh = Hydroxide()
        assert len(oh.atoms) == 2  # OH-
