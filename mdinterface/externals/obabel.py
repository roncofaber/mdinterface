#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:34:20 2025

@author: roncofaber
"""

# not repo
import logging
import numpy as np
import os
import tempfile
import ase
import ase.io

logger = logging.getLogger(__name__)

#%%

# charge_models = [ "eem", "mmff94", "gasteiger", "qeq", "qtpie", 
#                   "eem2015ha", "eem2015hm", "eem2015hn", 
#                   "eem2015ba", "eem2015bm", "eem2015bn" ]
def run_OBChargeModel(atoms, charge_type="eem"):
    
    try:
        import openbabel as ob
    except ImportError:
        raise ImportError("openbabel NOT found. Install it.")
    
    # Create an OBConversion object
    obConversion = ob.OBConversion()
    obConversion.SetInFormat("xyz")

    # Create an OBMol object
    mol = ob.OBMol()

    # Write and read the XYZ file using a secure temp file
    with tempfile.NamedTemporaryFile(suffix=".xyz", delete=False) as tmp:
        filename = tmp.name
    try:
        ase.io.write(filename, atoms)
        obConversion.ReadFile(mol, filename)
    finally:
        os.remove(filename)

    logger.info("OBabel charges: %s model,  %d atoms", charge_type, mol.NumAtoms())
    ob_charge_model = ob.OBChargeModel.FindType(charge_type)
    ob_charge_model.ComputeCharges(mol)
    charges = np.array(ob_charge_model.GetPartialCharges())
    logger.debug("  >> charges: sum=%.4f, min=%.4f, max=%.4f", charges.sum(), charges.min(), charges.max())

    return charges
