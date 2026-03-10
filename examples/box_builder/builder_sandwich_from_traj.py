#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BoxBuilder example: Pt(111) / polymer slab / Pt(111) sandwich loaded from a
LAMMPS dump trajectory.

The polymer slab positions come from a prior MD run rather than being built
from scratch.  The topology object (Specie / Polymer) is defined here — in a
real workflow it would be loaded from a pickle saved during chain building.

A minimal single-frame dump is provided in data/water_frame.dump so the
example is self-contained.  Replace DUMP_FILE with the path to your own
production dump to reproduce the real workflow.

Workflow
--------
1. Define / load the polymer topology (Specie or Polymer with FF params).
2. Point DUMP_FILE at a LAMMPS dump from a prior MD run.
3. Read the target frame with read_lammps_nth_frame — streams the file and
   never loads the full trajectory into memory.
4. Inject those positions into the topology object with update_positions().
5. Stack everything with BoxBuilder as usual.

Author: roncofaber
"""

from pathlib import Path

from mdinterface import BoxBuilder
from mdinterface.database import Water, Metal111
from mdinterface.io import read_lammps_nth_frame

#%% Paths

# Minimal 1-frame dump shipped with this example.
# In a real workflow: DUMP_FILE = "/path/to/production.dump"
DUMP_FILE = Path(__file__).parent / "data" / "water_frame.dump"

#%% Define species
# Water is used here as a stand-in for a Nafion polymer.
# In practice you would do:  polymer = pickle.load(open("nafi_20.pkl", "rb"))

water    = Water(model="ewald")
platinum = Metal111("Pt")

#%% Read one frame from the dump
# frame=-1  -> last frame (default)
# frame=0   -> first frame
# The file is streamed; only the requested frame is parsed.

last_frame = read_lammps_nth_frame(DUMP_FILE, frame=-1)

#%% Inject trajectory positions into the topology object
# prune_z=True trims any vacuum that accumulated during MD along z.
# The cell is also updated, so BoxBuilder picks up the correct XY dimensions.

water.update_positions(atoms=last_frame, prune_z=True)

#%% Set up simulation box
# xysize is a tiling hint for the electrode slabs.
# match_cell=water overrides the final XY to the polymer's exact cell,
# so the starting guess here only needs to be in the right ballpark.

simbox = BoxBuilder(
    xysize=[20, 20],
    verbose=True,
)

#%% Add layers (bottom -> top)

simbox.add_slab(platinum, nlayers=2)

simbox.add_slab(water, nlayers=1)       # polymer: cell already matches, no tiling

simbox.add_slab(platinum, nlayers=2)

#%% Build
# match_cell=water locks XY to the polymer's cell.
# Pt slabs are stretched to conform; the polymer itself is left unscaled.

simbox.build(
    padding=0.5,
    match_cell=platinum,
)

#%% Output

simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True)

atoms    = simbox.to_ase()
universe = simbox.universe
