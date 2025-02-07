#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 19:43:33 2025

@author: roncofaber
"""

# repo
from mdinterface.io.read import read_lammps_data_file

# not repo
import os
import ase
import ase.io
import random
import shutil
import subprocess

#%%

# helper folder to cleanup mess
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

# main ligpargen driver
def run_ligpargen(system, charge=None):
    """
    Runs the ligpargen command for the given xyz file.

    Parameters:
    system (ase.Atoms): The atoms system to be processed.

    Returns:
    ase.Atoms: The atoms object from the generated lammps file.
    """
    
    if "BOSSdir" not in os.environ:
        mdint = os.environ["MDINT_CONFIG_DIR"]
        print("BOSS was NOT found. Please either:")
        print("os.environ['BOSSdir'] = '/path/to/your/boss'")
        print("add 'BOSSdir = /path/to/your/boss' in [settings] in the")
        print(f"{mdint}/config.ini file.")
        
    # write file and run ligpargen
    random_number = random.randint(10000000, 99999999)
    folder_name = f"test_{random_number}"
    os.makedirs(folder_name, exist_ok=True)
    
    ase.io.write(f"{random_number}.xyz", system)
    
    # define ligpargen command
    if charge is None:
        ligpargen_command = f"ligpargen -i {random_number}.xyz -p {folder_name} -debug -o 0 -cgen CM1A"
    else:
        ligpargen_command = f"ligpargen -i {random_number}.xyz -p {folder_name} -debug -o 0 -c {charge} -cgen CM1A"

    try:
        subprocess.run(ligpargen_command, shell=True, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the command: {e}")
        # cleanup stuff
        cleanup(folder_name, f"{random_number}")
        os.remove(f"{random_number}.xyz")
        raise

    system, atoms, bonds, angles, dihedrals, impropers = read_lammps_data_file(
        f"{folder_name}/{random_number}.lammps.lmp")
    
    # cleanup stuff
    cleanup(folder_name, f"{random_number}")
    os.remove(f"{random_number}.xyz")
    
    return system, atoms, bonds, angles, dihedrals, impropers
