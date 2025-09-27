#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenBabel integration for charge calculation.

Interface to OpenBabel library for calculating atomic charges using
various charge models including EEM, MMFF94, Gasteiger, and others.

Author: Fabrice Roncoroni
Created: 2025-02-03
"""

import os
import random

import ase
import ase.io

# not repo
import numpy as np

# %%


# charge_models = [ "eem", "mmff94", "gasteiger", "qeq", "qtpie",
#                   "eem2015ha", "eem2015hm", "eem2015hn",
#                   "eem2015ba", "eem2015bm", "eem2015bn" ]
def run_OBChargeModel(atoms, charge_type="eem"):

    try:
        import openbabel as ob
    except ImportError:
        print("openbabel NOT found. Install it.")
        raise

    # Create an OBConversion object
    obConversion = ob.OBConversion()
    obConversion.SetInFormat("xyz")

    # Create an OBMol object
    mol = ob.OBMol()

    # Write and read the XYZ file
    # Generate a random 8-digit integer
    random_number = random.randint(10000000, 99999999)

    # Convert to obabel format
    filename = f"tmp_{random_number}.xyz"
    ase.io.write(filename, atoms)
    obConversion.ReadFile(mol, filename)
    os.remove(filename)

    ob_charge_model = ob.OBChargeModel.FindType(charge_type)
    ob_charge_model.ComputeCharges(mol)
    charges = ob_charge_model.GetPartialCharges()

    return np.array(charges)
