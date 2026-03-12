"""
mdinterface: Build Interface Systems for Molecular Dynamics Simulations

`mdinterface` is a Python package designed to build systems for Molecular Dynamics (MD) simulations.
Initially developed to construct electrolyte/electrode interfaces, it is also well-suited for generating MD boxes of liquids, electrolyte systems, and polymer networks.

"""

__version__ = '1.5.0'
__date__ = '09 Mar. 2026'
__author__ = 'Fabrice Roncoroni'
__all__ = ['SimulationBox', 'SimCell', 'BoxBuilder', "Specie", "Polymer"]

from .simulationbox import SimulationBox
from .build.builder import SimCell, BoxBuilder
from .core.specie import Specie
from .core.polymer import Polymer
from .utils.logger import set_verbosity

# load configuration file
from .config import load_config
load_config()
