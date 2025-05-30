#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:21:04 2024

@author: roncofaber
"""

import ase
import ase.io
import ase.visualize
from ase import units

import numpy as np
import glob
import os

from ase.calculators.singlepoint import SinglePointCalculator
#%%

def find_duplicate_frames(traj):
    tot_frames = len(traj)
    nstep = np.inf
    to_remove = []
    for cc, frame in enumerate(reversed(traj)):
        new_step = frame.info["i"]
        if new_step >= nstep:
            to_remove.append(tot_frames-cc-1)
        else:
            nstep = new_step
    return to_remove

def del_list_inplace(l, id_to_del):
    for i in sorted(id_to_del, reverse=True):
        del(l[i])

def get_md_data(filename, symbols, natoms):

    with open(filename, "r") as fin:
        
        results = []
        for line in fin:
            cline = line.split()
    
            if cline[0] in symbols:
                results.append(cline[1:])
    
    return np.array(results, dtype=float).reshape((-1, natoms, 3))
        
class CP2KTraj(object):
    
    def __init__(self, filename):
        
        self._filename = os.path.basename(filename)
        self._pathname = os.path.dirname(filename) + "/"
        self._prefix   = self._filename.split("-")[0]
        self._md_str   = self._filename.split("-")[1]
        
        self._traj = self.read_trajectory()
        
        
        return
    
    def read_trajectory(self):
        
        traj = ase.io.read(self._pathname + self._filename , index=":")
        
        natoms = len(traj[0])
        symbols = np.unique(traj[0].get_chemical_symbols())
        nframes = len(traj)
        
        read = {"cell": False, "vel": False, "frc": False, "stress": False}
        md_files = glob.glob(self._pathname + "/{}-{}*".format(self._prefix, self._md_str))
        
        cell, stress, vels, frcs = nframes*[None], nframes*[None],\
            nframes*[None], nframes*[None]
        for file in md_files:
            if "cell" in file and not read["cell"]:
                cell = np.atleast_2d(np.loadtxt(file))
                cell = cell[:,2:11].reshape((-1,3,3))
                read["cell"] = True
            elif "stress" in file and not read["stress"]:
                stress = np.atleast_2d(np.loadtxt(file))
                stress = units.bar*stress[:,2:].reshape((-1,3,3))
                read["stress"] = True
            elif "vel" in file and not read["vel"]:
                conv = units.Bohr/(units.AUT/units.fs) #convert to Ang/fs
                vels = conv*get_md_data(file, symbols, natoms)
                read["vel"] = True
            elif "frc" in file and not read["frc"]:
                conv = units.Ha/units.Bohr # convert to eV/Ang
                frcs = conv*get_md_data(file, symbols, natoms)
                read["frc"] = True
        
        # re-add last cell for structure opt
        if len(cell) < len(traj):
            cell = np.vstack([cell, cell[-1][np.newaxis, ...]])
        if len(stress) < len(traj):
            stress = np.vstack([stress, stress[-1][np.newaxis, ...]])
        
        
        # update traj
        for frame, tcell, tvel, tfrc, tstress in zip(traj, cell, vels, frcs, stress):
            
            frame.set_cell(tcell)
            # frame.set_pbc(True)
            
            energy = units.Ha*frame.info["E"]
            
            calc = SinglePointCalculator(frame,
                                         energy = energy,
                                         forces = tfrc,
                                         stress = tstress,
                                         )
            frame.calc = calc
            
            if tvel is not None:
                frame.set_velocities(tvel)
                
            
        # remove duplicate frames
        to_remove = find_duplicate_frames(traj)
        del_list_inplace(traj, to_remove)
        
        return traj

    @property
    def traj(self):
        return self._traj
    
    def get_energy(self, relative=True):
        energy = np.array([frame.get_potential_energy() for frame in self.traj])
        if relative:
            energy -= energy[-1]
        return energy
    
    def view(self):
        ase.visualize.view(self.traj)
        return
