#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 09:35:01 2024

@author: roncofaber
"""


from simulationbox import SimulationBox
from specie import Specie
from database import Water

import ase
import ase.visualize


#%% set up simulation box

# solvent
mol = Water()


#particle
particle = ase.io.read("Au19_Inital_30Vac_Center.cif")
particle.set_pbc(None)
particle.set_cell(None)
gold     = Specie(particle, )


# interface slab
slab = ase.io.read("TiO2_5Layers_Lowest3LayersFIxed_Optimized.cif")
slab.cell[2,2] = 18.5 # cut it here
interface = Specie(slab)


#%%

# frame = ase.io.read("/home/roncofaber/HPC_SCRATCH/17_Salmeron/01_LAMMPS/10_metadynamics/06_Cs2SO4_H2O_bulk_30x30x30/positions.xyz")

simobj = SimulationBox(
    solvent   = mol,
    solute    = [gold],
    # interface = interface
    )


#%%

solvent_vol = [15, 15, 15] # size of solvent \AA
solvent_rho = 0.99713          # density g/cm3
nions       = [1]            # number of ions
nlayers     = 1           # number of interface layers

system = simobj.make_simulation_box(solvent_vol, solvent_rho,
                                    nions = nions,
                                    # conmodel = cm,
                                    layers=nlayers, 
                                    to_ase = True,
                                    # padding=-1.5,
                                    write_data=True,
                                    layered=True,
                                    ion_pos='center', #can be center (center of box), left closer to iterface or None
                                    # hijack=lmptraj[-1]
                                    )

ase.io.write("POSCAR", system)

# simobj.write_lammps_file(system)#, coeff_file="coeff.lammps")
