#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 18:12:07 2025

@author: roncofaber
"""

import logging
import numpy as np

# internal stuff
from mdinterface.utils.graphs import find_equivalent_atoms

logger = logging.getLogger(__name__)

#%%

def calculate_RESP_charges(specie, basis='def2-svpd', xc="b3lyp", calc_type="RKS",
                           gpu=True, optimize=False, maxit=250):
    
    try:
        from gpu4pyscf.pop import esp
        from pymbxas.build.structure import ase_to_mole, mole_to_ase
        from pymbxas.build.input_pyscf import make_pyscf_calculator
        from pymbxas.md.solvers import Geometry_optimizer
    except ImportError:
        raise ImportError("You need to install gpu4pyscf AND pymbxas to run this.")
    
    logger.info("RESP charges: %d atoms,  basis=%s,  xc=%s,  optimize=%s",
                len(specie.atoms), basis, xc, optimize)

    # make pyscf mol
    mol = ase_to_mole(specie.atoms, basis=basis)

    # generate calculator
    mf = make_pyscf_calculator(mol, xc=xc, calc_type=calc_type, gpu=gpu)

    atoms = specie._atoms
    if optimize:
        logger.info("  ├> geometry optimization (max 100 steps)")
        gopt = Geometry_optimizer(mf)
        gopt.optimize(100)
        mol = gopt.mol_eq
        atoms = mole_to_ase(mol)
        logger.info("  ├> geometry optimized")

    mol.cart = True    # PySCF uses spherical basis by default

    # run calculation and calculate density matrix
    mf = make_pyscf_calculator(mol, xc=xc, calc_type=calc_type, gpu=gpu)
    logger.debug("  ├> SCF kernel")
    mf.kernel()
    dm = mf.make_rdm1()

    # RESP charge // first stage fitting
    logger.debug("  ├> RESP stage 1")
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
    logger.debug("  ├> RESP stage 2 (%d equal constraints)", len(equal_constraints))
    q2 = esp.resp_solve(mol, dm, resp_a=5e-4, resp_b=0.1, tol=1e-7,
                        sum_constraints=sum_constraints,
                        equal_constraints=equal_constraints, maxit=maxit)

    logger.info("  └─> done: sum=%.4f,  min=%.4f,  max=%.4f", q2.sum(), q2.min(), q2.max())
    return q2, atoms
