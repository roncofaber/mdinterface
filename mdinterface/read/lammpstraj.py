#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 11:51:21 2024

@author: roncofaber
"""

# ASE stuff
import ase
import ase.io
import ase.visualize

# package stuff
from .read import read_data_file, traj2chunks, dump2ase
from .trajectory import Trajectory

# parallel computation
import multiprocessing
from functools import partial
from multiprocessing.pool import Pool

#%%

class LammpsTraj(Trajectory):
    
    @staticmethod
    def _read_trajectory(filename, datafile, index, parallel=False, every=1):
        
        if datafile is not None:
            lmpdata = read_data_file(datafile, type2sym=True)
            
            data = {
                "labels"       : lmpdata[0],
                "symbols"      : lmpdata[1],
                "mol_idx"      : lmpdata[2],
                "mol_typ"      : lmpdata[3],
                "mol_names"    : lmpdata[4],
                "ato_typ"      : lmpdata[5],
                "connectivity" : lmpdata[6],
                "charges"      : lmpdata[7],
                }
            
        else:
            data = {
                "labels"  : None,
                "charges" : None,
                }
        
        if parallel:
            # Determine the number of processes to use (often based on CPU cores)
            num_processes = multiprocessing.cpu_count() - 2
            
            # Make chunks iterable
            chunks = traj2chunks(filename, every=every)
            
            pfunc = partial(dump2ase, specorder=data["labels"])
            # Create a Pool of processes
            
            traj = []
            with Pool(processes=num_processes) as pool:
                # Use pool.map() to apply the function to each item in parallel
                for frame in pool.imap(pfunc, chunks):
                    traj.append(frame)
            
        else:
            traj = ase.io.read(filename, format="lammps-dump-text",
                                index=index, specorder=data["labels"])
            
        if 'initial_charges' not in traj[0].arrays and data is not None:
            [ii.set_initial_charges(data["charges"]) for ii in traj]
    
        return traj, data