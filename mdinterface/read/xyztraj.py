#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 30 14:36:46 2025

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

class XYZTraj(Trajectory):
    
    @staticmethod
    def _read_trajectory(filename, datafile, index, parallel=False, every=1):
        
        traj = ase.io.read(filename)
        data = None
    
        return traj, data