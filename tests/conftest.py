"""PyTest configuration and fixtures for mdinterface tests."""

import pytest
import tempfile
import shutil
from pathlib import Path
from typing import Generator

import ase
import ase.build
import numpy as np


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test files."""
    temp_path = Path(tempfile.mkdtemp())
    try:
        yield temp_path
    finally:
        shutil.rmtree(temp_path)


@pytest.fixture
def water_molecule() -> ase.Atoms:
    """Create a simple water molecule for testing."""
    return ase.build.molecule('H2O')


@pytest.fixture
def methane_molecule() -> ase.Atoms:
    """Create a methane molecule for testing."""
    return ase.build.molecule('CH4')


@pytest.fixture
def simple_atomic_positions() -> np.ndarray:
    """Simple atomic positions for testing."""
    return np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])


@pytest.fixture
def simple_atomic_numbers() -> np.ndarray:
    """Simple atomic numbers for testing."""
    return np.array([1, 6, 8])  # H, C, O


@pytest.fixture
def mock_charges() -> np.ndarray:
    """Mock atomic charges for testing."""
    return np.array([-0.834, 0.417, 0.417])