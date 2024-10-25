#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 17:17:30 2024

@author: roncofaber
"""

import numpy as np
import ase
import ase.build

#%%

def build_polymer(system, substitute, nrep, target_distance=1.600):
    """
    Build a polymer by replicating a monomer structure and replacing 'X' atoms with a substitute.

    Parameters:
    - system (ase.Atoms): The monomer structure as an ASE Atoms object.
    - substitute (str): Chemical symbol to replace 'X' with after replication.
    - nrep (int): Number of times to replicate the monomer.
    - target_distance (float, optional): Desired bond distance at the joining points (in angstroms). 
      If None, the original distance between 'X' atoms is used.

    Returns:
    - sout (ase.Atoms): The replicated and modified polymer structure.
    """
    
    system = system.copy()
    
    # Check if nrep is at least 1
    if nrep < 1:
        raise ValueError("nrep must be at least 1")

    # Find the indices of the 'X' atoms
    x_idxs = np.where(np.array(system.get_chemical_symbols()) == "X")[0]
    if len(x_idxs) != 2:
        raise ValueError("The monomer should contain exactly two 'X' atoms.")

    # Find atoms connected to X
    ini_idx = ase.build.connected_indices(system, x_idxs[0])[1]
    end_idx = ase.build.connected_indices(system, x_idxs[1])[1]

    bnd_vec = system.get_distance(end_idx, x_idxs[1], vector=True)
    bnd_vec = target_distance * bnd_vec / np.linalg.norm(bnd_vec)

    ini_pos = system.get_positions()[ini_idx]
    end_pos = system.get_positions()[end_idx] + bnd_vec

    # Get the vector between the two 'X' atoms
    X_vec = end_pos - ini_pos
    
    # remember connecting points
    is_connected = np.array(len(system)*[False])
    is_connected[ini_idx] = True
    is_connected[end_idx] = True
    system.new_array("is_connected", is_connected)
    
    # Start with a copy of the system
    sout = system.copy()
    sout.new_array("mon_id", np.array(len(sout)*[0]))
    
    # Replicate the monomer
    if nrep > 1:
        del sout[x_idxs[1]]  # Remove the second 'X' atom in the first monomer

    for ii in range(nrep - 1):
        sadd = system.copy()

        if ii < nrep - 2:
            del sadd[x_idxs]  # Remove both 'X' atoms in the intermediate monomers
        else:
            del sadd[x_idxs[0]]  # Remove the first 'X' atom in the last monomer

        sadd.translate((ii + 1) * X_vec)
        sadd.set_array("mon_id", np.array(len(sadd)*[ii+1]))
        sout += sadd

    # Replace 'X' with the substitute
    symbols = np.array(sout.get_chemical_symbols())
    symbols[symbols == "X"] = substitute
    sout.set_chemical_symbols(symbols)

    return ase.build.sort(sout)
