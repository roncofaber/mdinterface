#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 19:43:33 2025

@author: roncofaber
"""

# repo
from mdinterface.io.read import read_lammps_data_file

# not repo
import logging
import os
import ase
import ase.io
import tempfile
import shutil
import subprocess

logger = logging.getLogger(__name__)

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
            logger.warning("Log file %s not found in %s; directory not removed.", check_file, path)
    else:
        logger.warning("%s is neither a file nor a directory", path)

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
        logger.warning(
            "BOSSdir not set. Please either:\n"
            "  os.environ['BOSSdir'] = '/path/to/your/boss'\n"
            "  or add 'BOSSdir = /path/to/boss' in [settings] in %s/config.ini",
            mdint,
        )
        
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
        
        logger.debug("ligpargen completed successfully")
        logger.debug("ligpargen stdout:\n%s", result.stdout.decode())

    except subprocess.CalledProcessError as e:
        logger.error("ligpargen failed: %s", e)
        logger.debug("ligpargen stderr:\n%s", e.stderr.decode())
        
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
