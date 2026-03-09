#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Noble gas parameters for classical MD simulations.

@author: roncofaber
"""

import numpy as np

from mdinterface.core.specie import Specie
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral, Improper

import ase

#%%

# dict with standard noble gas parameters
# ffield format: epsilon [kcal/mol], sigma [Ang]
# Vrabec: https://pubs.acs.org/doi/10.1021/jp012542o

noble_gas_parameters = {
    "Ne": {
        "charge": 0.0,
        "ffield": {
            "vrabec": [0.0674, 2.8010],
        }
    },
    "Ar": {
        "charge": 0.0,
        "ffield": {
            "vrabec": [0.2321, 3.3952],
        }
    },
    "Kr": {
        "charge": 0.0,
        "ffield": {
            "vrabec": [0.3232, 3.6274],
        }
    },
    "Xe": {
        "charge": 0.0,
        "ffield": {
            "vrabec": [0.4522, 3.9011],
        }
    },
}

def lookup_noble_gas_parameters(element, ffield):
    """
    Look up LJ parameters for a noble gas.
    
    Args:
        element (str): Chemical symbol of noble gas
        ffield (str): Force field name
        
    Returns:
        tuple: (charge, lj_parameters)
    """
    try:
        charge = noble_gas_parameters[element]["charge"]
        lj = noble_gas_parameters[element]["ffield"][ffield.lower()]
        return charge, lj
    except KeyError as e:
        raise ValueError(f"Invalid element or force field: {e}")


class NobleGas(Specie):
    """
    Class representing a noble gas atom with LJ parameters.
    
    Args:
        element (str): The chemical symbol of the noble gas (He, Ne, Ar, Kr, Xe, Rn).
        ffield (str): The force field to use for parameters (default is "opls").
        **kwargs: Additional keyword arguments to pass to the Specie superclass.
    """
    
    def __init__(self, element, ffield="vrabec", **kwargs):
        
        charge, lj_params = lookup_noble_gas_parameters(element, ffield)
        
        # Create single atom
        atom = ase.Atoms(element, positions=[(0.0, 0.0, 0.0)])
        
        # Set up LJ parameters
        lj = {element: lj_params}
        
        super().__init__(atoms=atom, charges=charge, lj=lj, **kwargs)
        
        return


# Convenience classes for each noble gas

class Neon(NobleGas):
    def __init__(self, ffield="vrabec", **kwargs):
        super().__init__("Ne", ffield=ffield, **kwargs)


class Argon(NobleGas):
    def __init__(self, ffield="vrabec", **kwargs):
        super().__init__("Ar", ffield=ffield, **kwargs)


class Krypton(NobleGas):
    def __init__(self, ffield="vrabec", **kwargs):
        super().__init__("Kr", ffield=ffield, **kwargs)


class Xenon(NobleGas):
    def __init__(self, ffield="vrabec", **kwargs):
        super().__init__("Xe", ffield=ffield, **kwargs)
