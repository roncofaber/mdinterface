"""
mdinterface: Build Interface Systems for Molecular Dynamics Simulations

`mdinterface` is a Python package designed to build systems for Molecular Dynamics (MD) simulations.
Initially developed to construct electrolyte/electrode interfaces, it is also well-suited for generating MD boxes of liquids, electrolyte systems, and polymer networks.

"""

__version__ = "1.3.0"
__date__ = "17 Jul. 2025"
__author__ = "Fabrice Roncoroni"
__all__ = ["SimulationBox", "Specie", "Polymer", "setup_logging"]

# Load configuration file
from .config import load_config
from .core.polymer import Polymer
from .core.specie import Specie

# Setup logging first
from .logging_config import setup_logging

# Import main classes
from .simulationbox import SimulationBox

load_config()
