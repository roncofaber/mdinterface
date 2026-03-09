#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BoxBuilder example: multi-layer cell with heterogeneous electrolyte regions.

This setup is impractical to express with the old SimulationBox API (which
only supports one interface, one enderface, and one miderface slot). With
BoxBuilder you can stack as many slabs and solvent regions as needed.

Layout (bottom to top):
    Au(111) slab        — bottom electrode
    NaCl solution       — left-half electrolyte
    Pt(111) slab        — middle electrode
    KF solution         — right-half electrolyte
    Pt(111) slab        — top electrode
    vacuum gap

Author: roncofaber
"""

from mdinterface import BoxBuilder
from mdinterface.database import Water, Ion, Metal111

#%% Define species

water = Water(model="ewald")

na = Ion("Na", ffield="Cheatham")
cl = Ion("Cl", ffield="Cheatham")
k  = Ion("K",  ffield="Cheatham")
f  = Ion("F",  ffield="Dang")

gold     = Metal111("Au")
platinum = Metal111("Pt")

#%% Build the five-layer stack

system = (
    BoxBuilder(xysize=[20, 20])
        # bottom electrode
        .add_slab(gold, nlayers=3)
        # first electrolyte region
        .add_solvent(water, ions=[na, cl], nions=[4, 4], zdim=20, density=1.0)
        # middle electrode
        .add_slab(platinum, nlayers=2)
        # second electrolyte region (different salt)
        .add_solvent(water, ions=[k, f],  nions=[4, 4], zdim=20, density=1.0)
        # top electrode
        .add_slab(platinum, nlayers=2)
        # vacuum gap to avoid interaction across the periodic boundary
        .add_vacuum(zdim=10)
        .build(padding=0.5)
        .write_lammps("data_multilayer.lammps", atom_style="full", write_coeff=True)
)
