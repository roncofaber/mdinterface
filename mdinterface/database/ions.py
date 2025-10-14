#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 11:47:01 2025

@author: roncofaber
"""

import numpy as np

from mdinterface.core.specie import Specie
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral, Improper

import ase
#%%

# dict with standard ions parameters. Normally for TIP3P water!
# ffield format: epsilon [kcal/mol], sigma [Ang]
# Aqvist   : https://pubs.acs.org/doi/10.1021/jp8001614 (check Simulation Conditions paragraph)
# Jorgensen: https://pubs.acs.org/doi/10.1021/ct600252r
# Cheatham : https://pubs.acs.org/doi/10.1021/jp8001614
# Sengupta : https://pubs.acs.org/doi/10.1021/acs.jcim.0c01390 (12-6 HFE param)
# Dang     : https://doi.org/10.1063/1.462555, https://doi.org/10.1063/1.466363
# OPLS-AA (Cl) : https://pubs.acs.org/doi/10.1021/ct900009a
    
ions_parameters = {
    "F": {
        "charge": -1.0,
        "ffield": {
            "jorgensen": [0.71000, 3.0500],
            "cheatham" : [0.00336, 4.1035],
            "sengupta" : [0.24140, 3.2678],
            "dang"     : [0.20000, 3.1680] # was [0.18000, 3.1180]
        }
    },
    "Cl": {
        "charge": -1.0,
        "ffield": {
            "aqvist"   : [   None,   None],
            "jorgensen": [0.71000, 4.0200],
            "cheatham" : [0.03559, 4.4777],
            "sengupta" : [0.63803, 4.0981],
            "dang"     : [0.10000, 4.4500],
            "opls-aa"  : [0.14800, 3.7700]
        }
    },
    "Br": {
        "charge": -1.0,
        "ffield": {
            "jorgensen": [0.71000, 4.2800],
            "cheatham" : [0.05866, 4.6469],
            "sengupta" : [0.75027, 4.4206],
        }
    },
    "I": {
        "charge": -1.0,
        "ffield": {
            "jorgensen": [0.71000, 4.8100],
            "cheatham" : [0.05368, 5.0959],
            "sengupta" : [0.86006, 4.8856],
        }
    },
    "Li": {
        "charge": +1.0,
        "ffield": {
            "aqvist"   : [0.01830, 2.02590],
            "jorgensen": [0.00050, 2.87000],
            "cheatham" : [0.02799, 1.82634],
            "sengupta" : [0.00282, 2.24506],
        }
    },
    "Na": {
        "charge": +1.0,
        "ffield": {
            "aqvist"   : [0.00277, 3.3284],
            "jorgensen": [0.00050, 4.0700],
            "cheatham" : [0.08744, 2.4393],
            "sengupta" : [0.02759, 2.5996],
            "dang"     : [0.13000, 2.3500]
        }
    },
    "K": {
        "charge": +1.0,
        "ffield": {
            "aqvist"   : [0.00033, 4.7360],
            "jorgensen": [0.00050, 5.1700],
            "cheatham" : [0.19368, 3.0380],
            "sengupta" : [0.14158, 3.0380],
        }
    },
    "Rb": {
        "charge": +1.0,
        "ffield": {
            "aqvist"   : [0.00017, 5.2670],
            "jorgensen": [0.00050, 5.6000],
            "cheatham" : [0.32782, 3.2304],
            "sengupta" : [0.21476, 3.2108],
        }
    },
    "Cs": {
        "charge": +1.0,
        "ffield": {
            "aqvist"   : [0.00008, 6.0492],
            "jorgensen": [0.00050, 6.2000],
            "cheatham" : [0.40654, 3.5208],
            "sengupta" : [0.36217, 3.5101],
        }
    }
}

def lookup_parameters(element, ffield):
    try:
        if ffield.lower() == "merz":
            ffield = "sengupta"
        charge = ions_parameters[element]["charge"]
        lj = ions_parameters[element]["ffield"][ffield.lower()]
        return charge, lj
    except KeyError as e:
        raise ValueError(f"Invalid element or force field: {e}")

# metal parameters from literature
class Ion(Specie):
    """
    Class representing an ion with specific parameters.
    
    Args:
        element (str): The chemical symbol of the ion.
        ffield (str): The force field to use for parameters (default is "Jorgensen").
        chg_scaling (float): Scaling factor for the charge (default is 0.8).
        **kwargs: Additional keyword arguments to pass to the Specie superclass.
    """
    
    def __init__(self, element, ffield="Jorgensen", chg_scaling=1.0, **kwargs):
        
        charge, lj = lookup_parameters(element, ffield)
        lj = {element: lj}
        super().__init__(element, charges=chg_scaling*charge, lj=lj, **kwargs)
        
        return


#perchlorate https://pubs.acs.org/doi/full/10.1021/jp801280s
class Perchlorate(Specie):
    def __init__(self, **kwargs):
        
        # Bond length
        bl = 1.506
        
        # Coordinates of ClO4^- anion
        coordinates = [
            (0.0, 0.0, 0.0),  # Chlorine atom at the origin
            (bl * np.sqrt(8/9), 0.0, -bl / 3),  # Oxygen atom 1
            (-bl * np.sqrt(2/9), bl * np.sqrt(2/3), -bl / 3),  # Oxygen atom 2
            (-bl * np.sqrt(2/9), -bl * np.sqrt(2/3), -bl / 3),  # Oxygen atom 3
            (0.0, 0.0, bl)   # Oxygen atom 4
        ]
        
        charges = [1.176, -0.544, -0.544, -0.544, -0.544]
        
        pclo = ase.Atoms("ClO4", positions=coordinates, charges=charges)
        
        b1 = Bond("Cl", "O", kr=757.286, r0=1.506)
        b2 = Bond("O",  "O", kr=33.445,  r0=2.459)
        a1 = Angle("O", "Cl", "O", kr=207.9, theta0=109.5)
        lj = {
            "Cl" : [0.1177, 3.5000],
            "O"  : [0.2099, 2.9000]
            }

        super().__init__(atoms=pclo, lj=lj, bonds=[b1, b2], angles=a1, cutoff=1.5, **kwargs)
        return


# hydronium parameters:
    # std:  https://pubs.acs.org/doi/pdf/10.1021/jp036842c
    # netz: https://refubium.fu-berlin.de/bitstream/handle/fub188/15473/1.4942771.pdf
class Hydronium(Specie):
    def __init__(self, **kwargs):
        
        # make ion by cheating and making NH3 first
        hyd = ase.build.molecule("NH3")
        hyd.set_chemical_symbols(["O", "H", "H", "H"])

        # factor 1/2 included in LAMMPS harmonic!
        b1 = Bond("O", "H", kr=1085.9565/2, r0=0.9820)
        a1 = Angle("H", "O", "H", kr=79.0263/2, theta0=113.4)
        
        # changed to -0.3819 -> -0.3818 to give +1 
        charges = [-0.3818, 0.4606, 0.4606, 0.4606]
        
        # divide by (2**(1/6)) to go from R0 to sig
        lj = {"O": [0.1848, 3.1655], "H": [0.010, 0.8018]}

        super().__init__(hyd, charges=charges, bonds=b1, angles=a1, lj=lj, **kwargs)
        return

# hydroxide parameters:
    # netz: https://refubium.fu-berlin.de/bitstream/handle/fub188/15473/1.4942771.pdf
class Hydroxide(Specie):
    def __init__(self, **kwargs):
        
        # make ion 
        hoh = ase.build.molecule("OH")

        # bond is fixed in paper #FIXME used hydronium bond with r0=1
        b1 = Bond("O", "H", kr=1085.9565/2, r0=1.000)
        
        charges = [-1.000, 0.000]
        
        # converted J to cal
        lj = {"O": [0.01195, 3.8100], "H": [0.000, 0.000]}

        super().__init__(hoh, charges=charges, bonds=b1, lj=lj, **kwargs)
        return