#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 18:12:07 2025

@author: roncofaber
"""

import numpy as np

# internal stuff
from mdinterface.utils.graphs import find_equivalent_atoms

#%%

def calculate_RESP_charges(specie, basis='def2-svpd', xc="b3lyp", calc_type="RKS",
                           gpu=True, optimize=False, maxit=250):
    
    try:
        from gpu4pyscf.pop import esp
        from pymbxas.build.structure import ase_to_mole, mole_to_ase
        from pymbxas.build.input_pyscf import make_pyscf_calculator
        from pymbxas.md.solvers import Geometry_optimizer
    except:
        raise ImportError("You need to install gpu4pyscf AND pymbxas to run this.")
    
    # make pyscf mol
    mol = ase_to_mole(specie.atoms, basis=basis)

    # generate calculator
    mf = make_pyscf_calculator(mol, xc=xc, calc_type=calc_type, gpu=gpu)
    
    atoms = specie._atoms
    if optimize:
        # optimize molecule
        gopt = Geometry_optimizer(mf)
        gopt.optimize(100)
        
        mol = gopt.mol_eq
        
        atoms = mole_to_ase(mol)
    
    mol.cart = True    # PySCF uses spherical basis by default
    
    # run calculation and calculate density matrix
    mf = make_pyscf_calculator(mol, xc=xc, calc_type=calc_type, gpu=gpu)
    
    mf.kernel()
    dm = mf.make_rdm1()
    
    # now, start RESP with ESP first!

    # ESP charge (do I need this?)
    # q0 = esp.esp_solve(mol, dm)

    # RESP charge // first stage fitting
    q1 = esp.resp_solve(mol, dm, maxit=maxit)

    # Add constraint: fix those charges in the second stage
    sum_constraints   = []
    equal_constraints = []

    _, idxs = find_equivalent_atoms(specie.graph)
    for idx in set(idxs):
        tmp_idx = np.concatenate(np.argwhere(idxs == idx))
        
        if len(tmp_idx) == 1:
            sum_constraints.append([q1[tmp_idx[0]], tmp_idx[0]])
        else:
            equal_constraints.append(tmp_idx.tolist())
            
    # RESP charge // second stage fitting
    q2 = esp.resp_solve(mol, dm, resp_a=5e-4, resp_b=0.1, tol=1e-7,
                        sum_constraints=sum_constraints,
                        equal_constraints=equal_constraints, maxit=maxit)

    return q2, atoms
