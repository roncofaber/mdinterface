#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 19:43:33 2025

@author: roncofaber
"""

import os
import random
import shutil
import subprocess
import numpy as np
import ase
import ase.io
from mdinterface.io.read import read_lammps_data_file

def select_connection_elements(specie):

    atoms = specie.atoms
    
    # Get elements where there is a connection
    centers = np.argwhere(atoms.arrays["is_connected"]).flatten()

    poi = []
    for center in centers:
        idxs = set(specie.find_relevant_distances(1, centers=center).flatten())
        if not any([ii in poi for ii in idxs]):
            poi.append(center)
        
    return poi

def cleanup(path, check_file):
    """
    Cleans up files or directories based on the specified path and check_file.

    Parameters:
    path (str): The path to the file or directory to be cleaned up.
    check_file (str): The file name to check for existence in the directory.
    """
    if os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        log_file_path = os.path.join(path, f"{check_file}.log")
        if os.path.exists(log_file_path):
            shutil.rmtree(path)
        else:
            print(f"Log file {check_file} not found in {path}. Directory not removed.")
    else:
        print(f"{path} is neither a file nor a directory")

def run_ligpargen(xyzfile):
    """
    Runs the ligpargen command for the given xyz file.

    Parameters:
    xyzfile (str): The name of the xyz file to be processed.

    Returns:
    ase.Atoms: The atoms object from the generated lammps file.
    """
    os.environ['BOSSdir'] = '/home/roncofaber/Software/boss'

    random_number = random.randint(10000000, 99999999)
    folder_name = f"test_{random_number}"
    ligpargen_command = f"ligpargen -i {xyzfile}.xyz -p {folder_name} -debug -o 0 -c 0 -cgen CM1A"

    os.makedirs(folder_name, exist_ok=True)

    try:
        subprocess.run(ligpargen_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the command: {e}")
        raise

    snippet, _, _, _, _, _ = read_lammps_data_file(f"{folder_name}/{xyzfile}.lammps.lmp")
    
    cleanup(folder_name, xyzfile)
    os.remove(f"{xyzfile}.xyz")
    
    return snippet

def make_snippet(specie, centers, Nmax, ending="F"):
    
    idxs = list(set(np.concatenate(specie.find_relevant_distances(Nmax, centers=centers))))
    edxs = list(set(np.concatenate(specie.find_relevant_distances(Nmax+1, centers=centers, Nmin=Nmax))) - set([centers]))
    
    snippet_idxs = np.array(idxs + edxs)

    chain = specie.atoms[idxs].copy()
    term = specie.atoms[edxs].copy()
    term.set_chemical_symbols(len(term) * [ending])
    snippet = ase.Atoms(chain + term)
    
    return snippet, snippet_idxs

def update_charges(specie, center, Nmax, charges, snippet_cache, ending="F"):
    """
    Updates the charges for the specified species.

    Parameters:
    specie (ase.Atoms): The species object containing the atoms.
    center (int): The center atom index.
    Nmax (int): The maximum number of neighbors to consider.
    charges (np.ndarray): The array of charges to be updated.
    ending (str): The chemical symbol for the terminal atoms.
    """
    
    ldxs = list(set(np.concatenate(specie.find_relevant_distances(4, centers=center))))
    
    snippet, snippet_idxs = make_snippet(specie, center, Nmax, ending=ending)

    # Check if snippet already exists in cache
    snippet_hash = ''.join(snippet.get_chemical_symbols())

    if snippet_hash in snippet_cache:
        cached_charges = snippet_cache[snippet_hash]
        mapping = [np.argwhere(snippet_idxs == ll)[0][0] for ll in ldxs]
        charges[ldxs] = cached_charges[mapping]
        return

    random_number = random.randint(10000000, 99999999)
    xyzfile = f"{random_number}"
    ase.io.write(f"{xyzfile}.xyz", snippet)

    nchain = run_ligpargen(xyzfile)

    mapping = [np.argwhere(snippet_idxs == ll)[0][0] for ll in ldxs]
    new_charges = nchain.get_initial_charges()[mapping]
    charges[ldxs] = new_charges

    # Store the new snippet and its charges in cache
    snippet_cache[snippet_hash] = nchain.get_initial_charges()
    
    return

def refine_charges(specie, Nmax=12, offset=False):
    """
    Refines the charges for the specified species.

    Parameters:
    specie (ase.Atoms): The species object containing the atoms.
    Nmax (int): The maximum number of neighbors to consider.
    offset (bool): Whether to apply an offset to the charges.

    Returns:
    np.ndarray: The refined charges.
    """
    charges = specie.atoms.get_initial_charges()
    centers = select_connection_elements(specie)
    print(centers)
    # Initialize a cache for storing snippets and their charges
    snippet_cache = {}

    for center in centers:
        update_charges(specie, center, Nmax, charges, snippet_cache)

    if offset:
        charges -= (charges.sum() / len(charges))

    return charges
