#!/usr/bin/env python3
"""Tests for mdinterface.core.topology classes."""

import pytest
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral, Improper


class TestAtom:
    def test_creation_minimal(self):
        a = Atom("C")
        assert a.symbol == "C"
        assert a.label == "C"

    def test_creation_with_lj(self):
        a = Atom("O", eps=0.102, sig=3.188)
        assert a.eps == pytest.approx(0.102)
        assert a.sig == pytest.approx(3.188)

    def test_custom_label(self):
        a = Atom("C", label="CT")
        assert a.symbol == "C"
        assert a.label == "CT"

    def test_copy_is_independent(self):
        a = Atom("N", eps=0.1, sig=3.0)
        b = a.copy()
        b.eps = 0.9
        assert a.eps == pytest.approx(0.1)

    def test_equality(self):
        a = Atom("O", eps=0.102, sig=3.188)
        b = Atom("O", eps=0.102, sig=3.188)
        assert a == b

    def test_inequality(self):
        a = Atom("O", eps=0.102, sig=3.188)
        b = Atom("O", eps=0.200, sig=3.188)
        assert a != b


class TestBond:
    def test_creation(self):
        b = Bond("O", "H", kr=450, r0=0.9572)
        assert b.symbols == ("O", "H")
        assert b.kr == pytest.approx(450)
        assert b.r0 == pytest.approx(0.9572)

    def test_equality_reversed(self):
        b1 = Bond("O", "H", kr=450, r0=0.9572)
        b2 = Bond("H", "O", kr=450, r0=0.9572)
        assert b1 == b2

    def test_copy_is_independent(self):
        b = Bond("O", "H", kr=450, r0=0.9572)
        c = b.copy()
        c.kr = 999
        assert b.kr == pytest.approx(450)

    def test_values_property(self):
        b = Bond("C", "C", kr=300, r0=1.54)
        assert b.values == (300, 1.54)


class TestAngle:
    def test_creation(self):
        a = Angle("H", "O", "H", kr=55, theta0=104.52)
        assert a.symbols == ("H", "O", "H")
        assert a.kr == pytest.approx(55)
        assert a.theta0 == pytest.approx(104.52)

    def test_equality_reversed(self):
        a1 = Angle("H", "O", "H", kr=55, theta0=104.52)
        a2 = Angle("H", "O", "H", kr=55, theta0=104.52)
        assert a1 == a2


class TestDihedral:
    def test_creation(self):
        d = Dihedral("H", "C", "C", "H", A1=0.0, A2=0.0, A3=0.3)
        assert d.symbols == ("H", "C", "C", "H")
        assert d.values[2] == pytest.approx(0.3)

    def test_equality(self):
        d1 = Dihedral("H", "C", "C", "H", A1=1.0, A2=2.0)
        d2 = Dihedral("H", "C", "C", "H", A1=1.0, A2=2.0)
        assert d1 == d2


class TestImproper:
    def test_creation(self):
        i = Improper("C", "O", "N", "H", K=1.0, d=1, n=2)
        assert i.values == (1.0, 1, 2)
        assert i.symbols == ("C", "O", "N", "H")

    def test_invalid_d_raises(self):
        with pytest.raises(AssertionError):
            Improper("C", "O", "N", "H", K=1.0, d=0, n=2)

    def test_invalid_n_raises(self):
        with pytest.raises(AssertionError):
            Improper("C", "O", "N", "H", K=1.0, d=1, n=9)
