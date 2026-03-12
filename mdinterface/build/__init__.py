"""
Build subsystem: polymer assembly, PACKMOL packing, slab tiling, and SimCell.

Public API is re-exported here for convenience; import directly from submodules
for finer-grained access.
"""

from .polymerize import (build_polymer, start_chain, attach_to_chain,
                         calculate_overlap_and_gradient,
                         generate_random_normalized_vector,
                         generate_spherical_points,
                         optimize_monomer_rotation_gradient)
from .box import (populate_box, make_interface_slab, add_component)
from .solvent import (make_solvent_box, populate_solutes)
from .continuum2sim import discretize_concentration
from .builder import SimCell, BoxBuilder
