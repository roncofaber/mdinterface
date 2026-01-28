#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 11:23:00 2025

@author: roncofaber
"""

# ase stuff
import ase

#%%

def get_default_config():
    
    default_pyscfargs = {
        "basis"     : '6-31gs',#'def2-svpd',
        "xc"        : "b3lyp",
        "calc_type" : "UKS",
        "is_gpu"    : True,
        "spin"      : None,
        "charge"    : None,
        }
    
    
    return default_pyscfargs

def make_ase_calc_pyscf(specie, **pyscfargs):
    
    try:
        # pymbxas stuff
        from pymbxas.build.structure import ase_to_mole
        from pymbxas.build.input_pyscf import make_pyscf_calculator

        # gpu4pyscf stuff
        from gpu4pyscf.tools.ase_interface import PySCF
    except:
        raise ValueError("No PySCF and PyMBXAS stuff")
    
    # load default values
    pyscf_arguments = get_default_config()
    
    # update values
    pyscf_arguments.update(pyscfargs)
    
    if pyscf_arguments["charge"] is None:
        pyscf_arguments["charge"] = specie._tot_charge

    # make pyscf mol
    mol = ase_to_mole(specie.atoms, **pyscf_arguments)

    mol.cart = True    # PySCF uses spherical basis by default

    # generate pyscf calculator
    mf = make_pyscf_calculator(mol, **pyscfargs)
    
    # return ASE version
    return PySCF(method=mf)

def make_ase_calc_uma(specie, task_name="omol", model_name="uma-s-1p1"):
    
    try:
        from fairchem.core import pretrained_mlip, FAIRChemCalculator

        predictor = pretrained_mlip.get_predict_unit(model_name=model_name, device="cuda")
        calc = FAIRChemCalculator(predictor, task_name="omol")
    except:
        raise ValueError("No UMA and OMol stuff")
    
    return calc