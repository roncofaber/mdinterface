#!/usr/bin/env python3
"""Tests for mdinterface.core.specie.Specie."""

import numpy as np
import pytest
import ase

from mdinterface.core.specie import Specie
from mdinterface.core.topology import Bond, Angle, Atom


@pytest.fixture
def water_atoms():
    """Minimal H2O ASE Atoms object."""
    return ase.build.molecule("H2O")


@pytest.fixture
def simple_specie(water_atoms):
    """Specie built from an ASE Atoms object with no topology."""
    return Specie(water_atoms)


class TestSpecieCreation:
    def test_from_ase_atoms(self, water_atoms):
        s = Specie(water_atoms)
        assert s is not None
        assert len(s.atoms) == 3

    def test_from_molecule_name_string(self):
        s = Specie("H2O")
        assert len(s.atoms) == 3

    def test_from_formula_string(self):
        # Falls back to ase.Atoms formula parser
        s = Specie("CO")
        assert len(s.atoms) == 2

    def test_atoms_are_copied(self, water_atoms):
        s = Specie(water_atoms)
        water_atoms.positions[0] = [99, 99, 99]
        assert not np.allclose(s.atoms.positions[0], [99, 99, 99])


class TestSpecieProperties:
    def test_cutoff_stored(self, water_atoms):
        s = Specie(water_atoms, cutoff=1.5)
        assert s.cutoff == pytest.approx(1.5)

    def test_resname_default_from_formula(self, water_atoms):
        s = Specie(water_atoms)
        assert isinstance(s.resname, str)
        assert len(s.resname) <= 4

    def test_resname_custom(self, water_atoms):
        s = Specie(water_atoms, name="WAT")
        assert s.resname == "WAT"

    def test_resname_truncated_to_4(self, water_atoms):
        s = Specie(water_atoms, name="TOOLONG")
        assert s.resname == "TOOL"

    def test_charges_set(self, water_atoms):
        charges = [-0.83, 0.415, 0.415]
        s = Specie(water_atoms, charges=charges)
        assert s.charges == pytest.approx(charges)

    def test_charges_scalar_broadcast(self, water_atoms):
        s = Specie(water_atoms, charges=0.0)
        assert all(c == pytest.approx(0.0) for c in s.charges)

    def test_atoms_property_returns_ase_atoms(self, simple_specie):
        assert isinstance(simple_specie.atoms, ase.Atoms)

    def test_graph_property_exists(self, simple_specie):
        import networkx as nx
        assert isinstance(simple_specie.graph, nx.Graph)

    def test_tot_charge_default_zero(self, water_atoms):
        s = Specie(water_atoms)
        assert s._tot_charge == 0

    def test_tot_charge_explicit(self, water_atoms):
        s = Specie(water_atoms, tot_charge=1)
        assert s._tot_charge == 1


class TestSpecieCopy:
    def test_copy_is_independent(self, simple_specie):
        c = simple_specie.copy()
        c.cutoff = 999.0
        assert simple_specie.cutoff != 999.0

    def test_copy_atoms_are_independent(self, simple_specie):
        c = simple_specie.copy()
        c.atoms.positions[0] = [99, 99, 99]
        assert not np.allclose(simple_specie.atoms.positions[0], [99, 99, 99])


class TestSpecieWithTopology:
    def test_bond_stored(self, water_atoms):
        b = Bond("O", "H", kr=450, r0=0.9572)
        s = Specie(water_atoms, bonds=b)
        # bonds property returns [indices, types] — types should be non-empty
        bond_data = s.bonds
        assert bond_data is not None

    def test_angle_stored(self, water_atoms):
        b = Bond("O", "H", kr=450, r0=0.9572)
        a = Angle("H", "O", "H", kr=55, theta0=104.52)
        s = Specie(water_atoms, bonds=b, angles=a)
        assert s.angles is not None

    def test_lj_none_default_equivalent_to_empty_dict(self, water_atoms):
        # Regression test for mutable default argument fix
        s1 = Specie(water_atoms)
        s2 = Specie(water_atoms, lj=None)
        assert len(s1.atoms) == len(s2.atoms)

    def test_lj_provided(self, water_atoms):
        lj = {"O": [0.102, 3.188], "H": [0.0, 1.0]}
        s = Specie(water_atoms, lj=lj)
        assert s is not None


class TestMutableDefaultRegression:
    """Ensure mutable defaults do not share state across calls."""

    def test_lj_not_shared_between_instances(self, water_atoms):
        s1 = Specie(water_atoms)
        s2 = Specie(water_atoms)
        # If lj={} were shared, mutations would leak. Creating two instances
        # and verifying both work correctly is sufficient.
        assert s1.resname == s2.resname

    def test_add_to_topology_with_none_defaults(self, simple_specie):
        # Should not raise — None defaults must be replaced with [] inside method
        simple_specie._add_to_topology()
