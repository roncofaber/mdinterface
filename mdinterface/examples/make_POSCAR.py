#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 10:39:29 2026

@author: roncofaber
"""

from mdinterface import SimulationBox
from mdinterface.database import Water, Metal111, Ion

import ase
import ase.io
import ase.visualize

#%% set up simulation box

# solvent https://docs.lammps.org/Howto_tip3p.html (Ewald model, check Price paper)
mol = Water()

#ions FF 
# Na  = Ion("Na", ffield="Merz")
# Li  = Ion("Li", ffield="Merz")

# interface = Metal111("Au")

#%%

simobj = SimulationBox(
    solvent   = mol,
    # solute    = [Na, F],
    # interface = interface,
    # enderface = interface.copy()
    )

#%%

lat  = 15

# is cube, everything is size lat
xysize = [lat, lat]
zdim   = lat

layering = [
    # {"type": "interface", "nlayers" : 2 },
    {"type": "solvent", "rho": 1.0, "zdim": zdim},#, "concentration": 0.5},# "ion_pos": "left"},
    # {"type": "enderface", "nlayers" : 2 },
    ]

system = simobj.make_simulation_box(xysize,
                                    layering,
                                    to_ase = True,
                                    write_data = False,
                                    # ion_pos="fixed",
                                    # padding=0.5,
                                    # center_electrode = True,
                                    # vacuum=20,
                                    # match_cell=True,
                                    # filename=f"data.lammps_{zdim}",
                                    # atom_style="full",
                                    # write_coeff=True
                                    )

#%% write poscar

ase.io.write("POSCAR", system, sort=True, direct=True, vasp6=True)

