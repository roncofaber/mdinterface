#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:34:20 2025

@author: roncofaber
"""

# not repo
import numpy as np
import os
import ase
import ase.io
import random

#%%

# charge_models = [ "eem", "mmff94", "gasteiger", "qeq", "qtpie", 
#                   "eem2015ha", "eem2015hm", "eem2015hn", 
#                   "eem2015ba", "eem2015bm", "eem2015bn" ]
def run_OBChargeModel(atoms, charge_type="eem"):
    
    try:
        import openbabel as ob
    except:
        print("openbabel NOT found. Install it.")
    
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
