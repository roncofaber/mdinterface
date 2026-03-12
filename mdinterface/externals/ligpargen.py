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

    # all ligpargen files go in a temp dir; kept on failure for inspection
    tmpdir   = tempfile.mkdtemp(prefix="ligpargen_")
    mol_name = os.path.basename(tmpdir)
    xyz_file = os.path.join(tmpdir, f"{mol_name}.xyz")

    ase.io.write(xyz_file, system)

    # use relative filenames and cwd=tmpdir -- ligpargen does not accept
    # absolute paths for -i
    ligpargen_command = ["ligpargen", "-i", f"{mol_name}.xyz", "-p", tmpdir,
                         "-debug", "-o", "0", "-cgen", "CM1A"]
    if charge is not None:
        ligpargen_command.extend(["-c", str(charge)])

    try:
        result = subprocess.run(ligpargen_command, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                cwd=tmpdir)
        logger.debug("ligpargen completed successfully")
        logger.debug("ligpargen stdout:\n%s", result.stdout.decode())

    except subprocess.CalledProcessError as e:
        # write stdout/stderr into the temp dir for inspection, then keep it
        error_log = os.path.join(tmpdir, "error_log.txt")
        with open(error_log, "w") as fh:
            fh.write("STDOUT:\n" + e.stdout.decode() + "\n")
            fh.write("STDERR:\n" + e.stderr.decode() + "\n")
        logger.error("ligpargen failed; temp files kept at: %s", tmpdir)
        logger.debug("ligpargen stderr:\n%s", e.stderr.decode())
        raise

    # read result
    system, atoms, bonds, angles, dihedrals, impropers = read_lammps_data_file(
        os.path.join(tmpdir, f"{mol_name}.lammps.lmp"), is_snippet=is_snippet)

    # success -- clean up
    shutil.rmtree(tmpdir, ignore_errors=True)

    return system, atoms, bonds, angles, dihedrals, impropers
