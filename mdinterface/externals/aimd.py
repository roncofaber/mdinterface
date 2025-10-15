#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 11:59:31 2025

@author: roncofaber
"""

from ase import units
from ase.io import Trajectory
from ase.md import MDLogger
from ase.md.langevin import Langevin
from ase.md.bussi import Bussi
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

# Optional fairchem import
try:
    from fairchem.core import pretrained_mlip, FAIRChemCalculator
    FAIRCHEM_AVAILABLE = True
except ImportError:
    FAIRCHEM_AVAILABLE = False

#%%

def run_aimd(atoms, timestep=0.5, temperature_K=300, friction=0.1, steps=1000,
             trajectory=None, logfile=None, **kwargs):
    """
    Perform Ab Initio Molecular Dynamics (AIMD) using ASE and FAIRChem.

    Parameters
    ----------
    atoms : ase.Atoms
        The atomic structure to run the AIMD on.
    timestep : float, default 0.1
        Time step for the simulation in fs.
    temperature_K : float, default 300
        Target temperature for the Langevin dynamics in Kelvin.
    friction : float, default 0.001
        Frictional damping coefficient in 1/fs.
    steps : int, default 1000
        Number of time steps to run the AIMD.
    trajectory : str, optional
        Path to save the MD trajectory.
    logfile : str, optional
        Path to save the MD log.
    **kwargs
        Additional arguments passed to the Langevin integrator.

    Returns
    -------
    ase.Atoms
        The atomic structure after AIMD simulation.

    Raises
    ------
    ImportError
        If fairchem is not installed.
    ValueError
        If there's an error loading FAIRChem models.
    """

    if not FAIRCHEM_AVAILABLE:
        raise ImportError(
            "FAIRChem is required for AIMD simulations but is not installed. "
            "Install it with: pip install fairchem-core"
        )

    # Make a copy to avoid modifying the original
    aimd_atoms = atoms.copy()

    # Ensure the computation setup with a fairchem calculator
    try:
        predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
        calc = FAIRChemCalculator(predictor, task_name="omol")
    except Exception as e:
        raise ValueError(f"Error loading UMA or OMol models: {e}")

    aimd_atoms.calc = calc
    
    # start velocities
    MaxwellBoltzmannDistribution(aimd_atoms, temperature_K=temperature_K)
    
    # Initialize Bussi Dynamics
    dyn = Langevin(
            aimd_atoms,
            timestep=timestep * units.fs,  # convert to effective unit
            temperature_K=temperature_K,
            friction=friction / units.fs,  # convert to effective unit
            **kwargs
            )

    # define traj and logger
    if trajectory:
        npt_traj = Trajectory(trajectory, mode="w", atoms=aimd_atoms)
        dyn.attach(npt_traj.write, interval=1)
    if logfile:
        npt_log  = MDLogger(dyn, aimd_atoms, logfile, stress=False, peratom=False, mode="w")
        dyn.attach(npt_log, interval=1)

    # Run the dynamics
    dyn.run(steps)

    return aimd_atoms  # Returning the updated structure after simulation