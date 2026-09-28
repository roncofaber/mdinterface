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
def run_OBChargeModel(atoms, charge_type="eem", charge=None):
    """Estimate partial charges with Open Babel while preserving stored chemistry.

    Parameters
    ----------
    atoms : ase.Atoms
        Molecular coordinates and optional RDKit graph.
    charge_type : str, default "eem"
        Open Babel charge model identifier.
    charge : int, optional
        Total molecular charge. Defaults to the stored formal-charge sum.

    Returns
    -------
    numpy.ndarray
        Partial charges in elementary-charge units, in input atom order.

    Raises
    ------
    ValueError
        If the charge conflicts with stored chemistry or the model fails.
    """
    
    try:
        from openbabel import openbabel as ob
    except ImportError:
        raise ImportError("openbabel NOT found. Install it.")
    
    from rdkit import Chem
    from mdinterface.core.chemistry import stored_molecule

    # Create an OBConversion object
    obConversion = ob.OBConversion()
    obConversion.SetInFormat("xyz")

    # Create an OBMol object
    mol = ob.OBMol()

    chemical_mol = stored_molecule(atoms)
    if charge is None:
        charge = Chem.GetFormalCharge(chemical_mol) if chemical_mol is not None else int(np.sum(atoms.arrays.get("nominal_charge", [0])))
    if chemical_mol is not None and charge != Chem.GetFormalCharge(chemical_mol):
        raise ValueError("charge conflicts with the stored molecular charge.")
    if chemical_mol is not None:
        obConversion.SetInFormat("mol")
        if not obConversion.ReadString(mol, Chem.MolToMolBlock(chemical_mol)):
            raise ValueError("Open Babel could not read the molecular graph.")
    else:
        with tempfile.NamedTemporaryFile(suffix=".xyz", delete=False) as tmp:
            filename = tmp.name
        try:
            ase.io.write(filename, atoms)
            if not obConversion.ReadFile(mol, filename):
                raise ValueError("Open Babel could not read the molecular coordinates.")
        finally:
            os.remove(filename)

    mol.SetTotalCharge(int(charge))

    logger.info("OBabel charges: %s model,  %d atoms", charge_type, mol.NumAtoms())
    ob_charge_model = ob.OBChargeModel.FindType(charge_type)
    if ob_charge_model is None or not ob_charge_model.ComputeCharges(mol):
        raise ValueError(f"Open Babel charge model {charge_type!r} failed.")
    charges = np.array(ob_charge_model.GetPartialCharges())
    logger.debug("  >> charges: sum=%.4f, min=%.4f, max=%.4f", charges.sum(), charges.min(), charges.max())

    return charges
