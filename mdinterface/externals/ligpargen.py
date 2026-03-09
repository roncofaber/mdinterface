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
import tempfile
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

def run_ligpargen(system, charge=None, is_snippet=False):
    """
    Runs the ligpargen command for the given xyz file.

    Parameters:
    system (ase.Atoms): The atoms system to be processed.

    Returns:
    tuple: Containing system, atoms, bonds, angles, dihedrals, impropers.
    """
    
    if "BOSSdir" not in os.environ:
        mdint = os.environ["MDINT_CONFIG_DIR"]
        print("BOSS was NOT found. Please either:")
        print("os.environ['BOSSdir'] = '/path/to/your/boss'")
        print("add 'BOSSdir = /path/to/your/boss' in [settings] in the")
        print(f"{mdint}/config.ini file.")
        
    # Write the XYZ file and prepare to run ligpargen
    folder_name = tempfile.mkdtemp(prefix="ligpargen_")
    mol_name = os.path.basename(folder_name)
    xyz_file = f"{mol_name}.xyz"

    ase.io.write(xyz_file, system)

    # Define ligpargen command as a list to avoid shell injection
    ligpargen_command = ['ligpargen', '-i', xyz_file, '-p', folder_name,
                         '-debug', '-o', '0', '-cgen', 'CM1A']
    if charge is not None:
        ligpargen_command.extend(['-c', str(charge)])

    error_log_path = os.path.join(folder_name, "error_log.txt")  # Path for the error log

    try:
        # Run the command and capture both stdout and stderr
        result = subprocess.run(ligpargen_command, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Optionally log stdout for debugging
        output = result.stdout.decode()
        print("Command executed successfully, output captured:")
        print(output)
        
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the command: {e}")
        print(e.stderr.decode())
        
        # Write stderr output to the error log
        with open(error_log_path, 'w') as error_log_file:
            # To capture both stdout and stderr in the log file for further analysis
            error_log_file.write("STDOUT:\n" + e.stdout.decode() + "\n")
            error_log_file.write("STDERR:\n" + e.stderr.decode() + "\n")

        # Additional cleanup code can be executed here
        cleanup(folder_name, mol_name)
        os.remove(xyz_file)
        
        raise

    # Proceed to read lammps data file
    system, atoms, bonds, angles, dihedrals, impropers = read_lammps_data_file(
        f"{folder_name}/{mol_name}.lammps.lmp", is_snippet=is_snippet)

    # Additional cleanup code can be executed here
    cleanup(folder_name, mol_name)
    os.remove(xyz_file)
    
    return system, atoms, bonds, angles, dihedrals, impropers
