"""
mdinterface: Build Interface Systems for Molecular Dynamics Simulations

`mdinterface` is a Python package designed to build systems for Molecular Dynamics (MD) simulations.
Initially developed to construct electrolyte/electrode interfaces, it is also well-suited for generating MD boxes of liquids,
electrolyte systems, and polymer networks.

"""

__version__ = '1.1.0'
__date__ = '14 Jan. 2025'
__author__ = 'Fabrice Roncoroni'
__all__ = ['SimulationBox']

from .simulationbox import SimulationBox