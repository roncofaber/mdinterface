#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structure relaxation functions using ASE and UMA

Created on XXX

@author: roncofaber
"""

import ase
import numpy as np
from ase.optimize import BFGS, LBFGS, FIRE
import warnings

#%%

def relax_structure(atoms, optimizer='FIRE', fmax=0.05, steps=200,
                       trajectory=None, logfile=None, **kwargs):
    """
    Perform structure relaxation using ASE optimizers.

    Parameters
    ----------
    atoms : ase.Atoms
        The atomic structure to relax
    optimizer : str, default 'BFGS'
        Optimizer to use ('BFGS', 'LBFGS', 'FIRE')
    fmax : float, default 0.05
        Maximum force threshold for convergence (eV/Å)
    steps : int, default 200
        Maximum number of optimization steps
    trajectory : str, optional
        Path to save optimization trajectory
    logfile : str, optional
        Path to save optimization log
    **kwargs
        Additional arguments passed to the optimizer

    Returns
    -------
    relaxed_atoms : ase.Atoms
        The relaxed atomic structure
    converged : bool
        Whether the optimization converged
    """

    # Make a copy to avoid modifying the original
    relaxed_atoms = atoms.copy()

    # Select optimizer
    optimizer_dict = {
        'BFGS': BFGS,
        'LBFGS': LBFGS,
        'FIRE': FIRE
    }

    if optimizer not in optimizer_dict:
        raise ValueError(f"Unknown optimizer: {optimizer}. Available: {list(optimizer_dict.keys())}")
    
    try:
        from fairchem.core import pretrained_mlip, FAIRChemCalculator

        predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
        calc = FAIRChemCalculator(predictor, task_name="omol")
    except:
        raise ValueError("No UMA and OMol stuff")
    
    relaxed_atoms.calc = calc
    
    # Initialize optimizer
    opt_class = optimizer_dict[optimizer]
    
    # Create
    dyn = opt_class(relaxed_atoms, trajectory=trajectory, logfile=logfile, **kwargs)
    
    dyn.run(fmax, steps)

    return relaxed_atoms

