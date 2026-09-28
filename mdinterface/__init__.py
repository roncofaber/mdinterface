"""
mdinterface: Build Interface Systems for Molecular Dynamics Simulations

`mdinterface` is a Python package designed to build systems for Molecular Dynamics (MD) simulations.
Initially developed to construct electrolyte/electrode interfaces, it is also well-suited for generating MD boxes of liquids, electrolyte systems, and polymer networks.

"""

__version__ = '2.0.0'
__date__ = '28 Sep. 2026'
__author__ = 'Fabrice Roncoroni'
__all__ = ['SimCell', "Specie", "Polymer", "PackmolError", "LigParGenError"]

from .core.specie import Specie
from .core.polymer import Polymer
from .build.box import PackmolError
from .externals import LigParGenError
from .build.builder import SimCell
from .utils.logger import set_verbosity


def load_config():
    """Load user configuration without performing configuration I/O during import."""
    from .config import load_config as _load_config

    return _load_config()
