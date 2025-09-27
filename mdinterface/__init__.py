"""
mdinterface: Build Interface Systems for Molecular Dynamics Simulations.

A Python package designed to build systems for Molecular Dynamics (MD) simulations.
Initially developed to construct electrolyte/electrode interfaces, it is also well-suited
for generating MD boxes of liquids, electrolyte systems, and polymer networks.

Author: Fabrice Roncoroni
Created: 2024-04-19
"""

__version__ = "1.3.0"
__date__ = "2024-04-19"
__author__ = "Fabrice Roncoroni"
__all__ = ["SimulationBox", "Specie", "Polymer"]

# Load configuration file
from .config import load_config
from .core.polymer import Polymer
from .core.specie import Specie

# Import main classes
from .simulationbox import SimulationBox

load_config()
