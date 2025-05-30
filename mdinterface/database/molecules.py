#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 11:40:58 2025

@author: roncofaber
"""

from mdinterface.core.specie import Specie
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral, Improper

#%%
# solvent https://docs.lammps.org/Howto_tip3p.html (Ewald model)
class Water(Specie):
    def __init__(self, model="ewald", **kwargs):
        
        if model.lower() == "ewald":
            b1 = Bond("O", "H", kr=450, r0=0.9572)
            a1 = Angle("H", "O", "H", kr=55, theta0=104.52)
            charges = [-0.83, 0.415, 0.415]
            lj = {"O": [0.102, 3.188], "H": [0.0, 1.0]}
        
        if model.lower() == "charmm":
            b1 = Bond("O", "H", kr=450, r0=0.9572)
            a1 = Angle("H", "O", "H", kr=55, theta0=104.52)
            charges = [-0.834, 0.417, 0.417]
            lj = {"O": [0.1521, 3.1507], "H": [0.0460, 0.4]}

        super().__init__("H2O", charges=charges, bonds=b1, angles=a1, lj=lj, **kwargs)
        return

# 22 Apr. 2025 correction: all bond terms have been divider by 2:
#oxygen https://pubs.acs.org/doi/10.1021/acs.jctc.0c01132 /!\: divide sig by 2**(1/6)
class Oxygen(Specie):
    def __init__(self, **kwargs):
        
        b1 = Bond("O", "O", kr=1640.4/2, r0=1.2074)
        lj = {"O" : [0.1047, 2.9373]}

        super().__init__("O2", charges = 0.0, lj=lj, bonds=b1, **kwargs)
        return

#hydrogen https://pubs.acs.org/doi/10.1021/acs.jctc.0c01132 /!\: divide sig by 2**(1/6)
class Hydrogen(Specie):
    def __init__(self, Hset="std", **kwargs):
        
        b1 = Bond("H", "H", kr=700/2, r0=0.7414)
        
        if Hset.lower() == "std":   # standard 12-6 set
            lj = {"H" : [0.0153, 2.5996]}
        elif Hset.lower() == "alt": # alternative 12-6 set
            lj = {"H" : [0.0145, 2.8001]}

        super().__init__("H2", charges = 0.0, lj=lj, bonds=b1, **kwargs)
        return

#nitrogen https://pubs.acs.org/doi/10.1021/acs.jctc.0c01132 /!\: divide sig by 2**(1/6)
class Nitrogen(Specie):
    def __init__(self, **kwargs):
        
        b1 = Bond("N", "N", kr=3190/2, r0=1.0977)
        lj = {"N" : [0.0797, 3.2197]}

        super().__init__("N2", charges = 0.0, lj=lj, bonds=b1, **kwargs)
        return
