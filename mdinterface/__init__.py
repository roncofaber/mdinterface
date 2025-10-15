"""
mdinterface: Build Interface Systems for Molecular Dynamics Simulations

`mdinterface` is a Python package designed to build systems for Molecular Dynamics (MD) simulations.
Initially developed to construct electrolyte/electrode interfaces, it is also well-suited for generating MD boxes of liquids, electrolyte systems, and polymer networks.

"""

__version__ = '1.4.0'
__date__ = '14 Oct. 2025'
__author__ = 'Fabrice Roncoroni'
__all__ = ['SimulationBox', "Specie", "Polymer"]

from .simulationbox import SimulationBox
from .core.specie import Specie
from .core.polymer import Polymer

# load configuration file
from .config import load_config
load_config()
