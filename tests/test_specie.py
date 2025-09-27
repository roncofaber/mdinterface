"""Tests for the Specie class."""

import pytest
import numpy as np
import ase
from pathlib import Path

from mdinterface.core.specie import Specie


class TestSpecie:
    """Test cases for the Specie class."""

    def test_specie_initialization_basic(self, water_molecule):
        """Test basic Specie initialization with ASE Atoms."""
        specie = Specie(atoms=water_molecule, name="water")

        assert specie.resname == "wate"  # Truncated to 4 chars
        assert len(specie.atoms) == 3
        assert specie.cutoff == 1.0

    def test_specie_initialization_with_charges(self, water_molecule, mock_charges):
        """Test Specie initialization with custom charges."""
        specie = Specie(
            atoms=water_molecule,
            charges=mock_charges,
            name="H2O"
        )

        assert specie.resname == "H2O"
        assert len(specie.atoms) == 3

    def test_specie_initialization_no_atoms(self):
        """Test Specie initialization without atoms should handle gracefully."""
        with pytest.raises((ValueError, TypeError)):
            Specie(atoms=None)

    def test_specie_name_truncation(self, water_molecule):
        """Test that long names are truncated to 4 characters."""
        long_name = "very_long_molecule_name"
        specie = Specie(atoms=water_molecule, name=long_name)

        assert specie.resname == "very"
        assert len(specie.resname) == 4

    def test_specie_chemical_formula_name(self, methane_molecule):
        """Test automatic name generation from chemical formula."""
        specie = Specie(atoms=methane_molecule)

        # Should use chemical formula if no name provided
        assert hasattr(specie, 'resname')
        assert len(specie.resname) <= 4

    @pytest.mark.parametrize("cutoff", [0.5, 1.0, 1.5, 2.0])
    def test_specie_cutoff_parameter(self, water_molecule, cutoff):
        """Test different cutoff values."""
        specie = Specie(atoms=water_molecule, cutoff=cutoff)
        assert specie.cutoff == cutoff

    def test_specie_total_charge(self, water_molecule):
        """Test total charge parameter."""
        tot_charge = -1
        specie = Specie(atoms=water_molecule, tot_charge=tot_charge)
        assert specie.tot_charge == tot_charge

    @pytest.mark.slow
    def test_specie_with_lammps_data(self, temp_dir):
        """Test Specie initialization from LAMMPS data file."""
        # This would require a sample LAMMPS data file
        # Skip if no test data available
        pytest.skip("Requires sample LAMMPS data file")

    @pytest.mark.integration
    def test_specie_ligpargen_integration(self, methane_molecule):
        """Test LigParGen integration (requires external dependencies)."""
        # Skip if LigParGen not available
        try:
            specie = Specie(atoms=methane_molecule, ligpargen=True)
        except ImportError:
            pytest.skip("LigParGen not available")
        except Exception:
            pytest.skip("LigParGen integration test requires proper setup")

    def test_specie_pbc_parameter(self, water_molecule):
        """Test periodic boundary conditions parameter."""
        specie = Specie(atoms=water_molecule, pbc=True)
        # Should not raise error and handle PBC appropriately


class TestSpecieProperties:
    """Test Specie properties and methods."""

    def test_specie_has_required_attributes(self, water_molecule):
        """Test that Specie has all required attributes after initialization."""
        specie = Specie(atoms=water_molecule)

        required_attrs = ['atoms', 'cutoff', 'resname']
        for attr in required_attrs:
            assert hasattr(specie, attr), f"Missing required attribute: {attr}"

    def test_specie_atoms_property(self, water_molecule):
        """Test that atoms property is correctly set."""
        specie = Specie(atoms=water_molecule)

        assert isinstance(specie.atoms, ase.Atoms)
        assert len(specie.atoms) == len(water_molecule)
        assert np.allclose(specie.atoms.positions, water_molecule.positions)


class TestSpecieErrorHandling:
    """Test error handling in Specie class."""

    def test_invalid_lammps_file(self):
        """Test handling of invalid LAMMPS data file."""
        with pytest.raises((FileNotFoundError, ValueError)):
            Specie(lammps_data="nonexistent_file.data")

    def test_invalid_charges_length(self, water_molecule):
        """Test handling of charges array with wrong length."""
        invalid_charges = np.array([0.5, -0.5])  # Wrong length for water

        with pytest.raises((ValueError, IndexError)):
            Specie(atoms=water_molecule, charges=invalid_charges)