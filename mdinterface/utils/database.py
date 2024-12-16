#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 17:13:50 2024

@author: roncofaber
"""

import numpy as np

# internal modules
from mdinterface.core.specie import Specie
from mdinterface.core.topology import Bond, Angle, Atom, Dihedral, Improper

import ase.build

#%%

# graphene https://onlinelibrary.wiley.com/doi/10.1002/adma.201705791
class Graphene(Specie):
    def __init__(self, **kwargs):
        system = ase.build.graphene()
        system.cell[-1][-1] = 3.35
        g_b = Bond("C", "C", kr=469, r0=1.4)
        g_d = Dihedral("C", "C", "C", "C", A1=7.25, A2=0, A3=-7.25, A4=0, A5=0)
        g_i = Improper(a1="C", K=1.1, d=-1, n=2)
        g_a = Angle("C", "C", "C", kr=63, theta0=120)
        lj = {"C": [0.07, 3.54996412]}

        super().__init__(system, bonds=g_b, dihedrals=g_d, impropers=g_i, angles=g_a, lj=lj, **kwargs)
        
        return

# solvent https://docs.lammps.org/Howto_tip3p.html (Ewald model)
class Water(Specie):
    def __init__(self, model="ewald", **kwargs):
        
        if model.lower() == "ewald":
            b1 = Bond("O", "H", kr=450, r0=0.9572)
            a1 = Angle("H", "O", "H", kr=55, theta0=104.52)
            charges = [-0.83, 0.415, 0.415]
            lj = {"O": [0.102, 3.188], "H": [0.0, 1.0]}

        super().__init__("H2O", charges=charges, bonds=b1, angles=a1, lj=lj, **kwargs)
        return


#oxygen https://pubs.acs.org/doi/10.1021/acs.jctc.0c01132 /!\: divide sig by 2**(1/6)
class Oxygen(Specie):
    def __init__(self, **kwargs):
        
        b1 = Bond("O", "O", kr=1640.4, r0=1.2074)
        lj = {"O" : [0.1047, 2.9373]}

        super().__init__("O2", charges = 0.0, lj=lj, bonds=b1, **kwargs)
        return

#hydrogen https://pubs.acs.org/doi/10.1021/acs.jctc.0c01132 /!\: divide sig by 2**(1/6)
class Hydrogen(Specie):
    def __init__(self, Hset="std", **kwargs):
        
        b1 = Bond("H", "H", kr=700, r0=0.7414)
        
        if Hset.lower() == "std":   # standard 12-6 set
            lj = {"H" : [0.0153, 2.5996]}
        elif Hset.lower() == "alt": # alternative 12-6 set
            lj = {"H" : [0.0145, 2.8001]}

        super().__init__("H2", charges = 0.0, lj=lj, bonds=b1, **kwargs)
        return

#nitrogen https://pubs.acs.org/doi/10.1021/acs.jctc.0c01132 /!\: divide sig by 2**(1/6)
class Nitrogen(Specie):
    def __init__(self, **kwargs):
        
        b1 = Bond("N", "N", kr=3190, r0=1.0977)
        lj = {"N" : [0.0797, 3.2197]}

        super().__init__("N2", charges = 0.0, lj=lj, bonds=b1, **kwargs)
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
