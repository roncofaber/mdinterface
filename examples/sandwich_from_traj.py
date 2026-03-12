#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SimCell example: Pt(111) / polymer membrane / Pt(111) sandwich loaded from a
prior MD run.

The polymer positions come from the last frame of a production LAMMPS dump.
The topology (bonds, angles, dihedrals, impropers) is read from the
corresponding LAMMPS data file that was used to run that simulation.

This workflow is typical when building electrode/membrane/electrode sandwich
structures for a new simulation starting from an already-equilibrated membrane.

Workflow
--------
1. Read the polymer topology from the LAMMPS data file used in the prior run.
2. Read the last frame of the production dump with read_lammps_nth_frame.
3. Centre the polymer in XY so it sits in the middle of the new cell.
4. Inject those positions into the topology Specie with update_positions().
5. Stack everything with SimCell: Pt / polymer / Pt.

Author: roncofaber
"""

from pathlib import Path

import numpy as np

from mdinterface import SimCell, Specie
from mdinterface.database import Metal111
from mdinterface.io.read import read_lammps_data_file
from mdinterface.io import read_lammps_nth_frame

#%% Paths -- replace with your own files

# LAMMPS data file from the prior simulation (provides topology)
DATA_FILE = Path("data/membrane.data.lammps")

# Production dump from the prior simulation (provides equilibrated positions)
DUMP_FILE = Path("data/membrane.production.dump")

#%% Load polymer topology from the LAMMPS data file

monomer, atoms, bonds, angles, dihedrals, impropers = \
    read_lammps_data_file(DATA_FILE, pbc=True)

polymer = Specie(
    atoms     = monomer,
    name      = "MEMB",
    bonds     = bonds,
    angles    = angles,
    dihedrals = dihedrals,
    impropers = impropers,
    fix_missing = False,
)

#%% Read last frame of the production dump

# read_lammps_nth_frame seeks from the end of the file -- fast even for
# multi-GB trajectories.
last_frame = read_lammps_nth_frame(DUMP_FILE, frame=-1)

#%% Centre the polymer in XY

# Shift the heavy-atom centroid to the cell centre so it sits symmetrically
# between the two electrodes after stacking.
C_atoms = np.array(last_frame.get_chemical_symbols()) == "C"
centroid = last_frame.get_positions()[C_atoms].mean(axis=0)
cell_center = np.diag(last_frame.cell) / 2
last_frame.translate(cell_center - centroid)

#%% Inject equilibrated positions into the topology object

# prune_z=True removes any vacuum gap that accumulated during NPT MD.
# The cell is updated automatically, so SimCell picks up the correct XY.
polymer.update_positions(atoms=last_frame, prune_z=True)

#%% Set up simulation box

platinum = Metal111("Pt")


xsize = polymer.atoms.cell[0][0]
ysize = polymer.atoms.cell[1][1]

simbox = SimCell(xysize=[xsize, ysize], verbose=True)

#%% Add layers: bottom electrode / membrane / top electrode

simbox.add_slab(platinum, nlayers=2)

simbox.add_prebuilt(polymer)   # positions and cell already set -- no tiling

simbox.add_slab(platinum, nlayers=2)

#%% Build
# match_cell=polymer locks the final XY to the membrane cell.
# Pt slabs are stretched to match; the membrane itself is left unscaled.

simbox.build(
    padding    = 1.5,
    center     = False,
    match_cell = polymer,
)

#%% Output

simbox.write_lammps("data.lammps", atom_style="full", write_coeff=True)

ase_atoms = simbox.to_ase()
universe  = simbox.universe
