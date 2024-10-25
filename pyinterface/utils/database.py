#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 17:13:50 2024

@author: roncofaber
"""

# internal modules
from pyinterface.core.specie import Specie
from pyinterface.core.topology import Bond, Angle, Atom, Dihedral, Improper

import ase.build

#%%

# graphene https://onlinelibrary.wiley.com/doi/10.1002/adma.201705791
class Graphene(Specie):
    def __init__(self):
        system = ase.build.graphene()
        system.cell[-1][-1] = 3.35
        g_b = Bond("C", "C", kr=469, r0=1.4)
        g_d = Dihedral("C", "C", "C", "C", A1=7.25, A2=0, A3=-7.25, A4=0, A5=0)
        g_i = Improper("C", K=1.1, d=-1, n=2)
        g_a = Angle("C", "C", "C", kr=63, theta0=120)
        lj = {"C": [0.07, 3.54996412]}

        super().__init__(system, bonds=g_b, dihedrals=g_d, impropers=g_i, angles=g_a, lj=lj)
        
        return

# solvent https://docs.lammps.org/Howto_tip3p.html (Ewald model)
class Water(Specie):
    def __init__(self, model="ewald"):
        
        if model.lower() == "ewald":
            b1 = Bond("O", "H", kr=450, r0=0.9572)
            a1 = Angle("H", "O", "H", kr=55, theta0=104.52)
            charges = [-0.83, 0.415, 0.415]
            lj = {"O": [0.102, 3.188], "H": [0.0, 0.0]}

        super().__init__("H2O", charges=charges, bonds=b1, angles=a1, lj=lj)
        return


#oxygen https://pubs.acs.org/doi/10.1021/acs.jctc.0c01132 /!\: divide sig by 2**(1/6)
class Oxygen(Specie):
    def __init__(self):
        
        b1 = Bond("O", "O", kr=1640.4, r0=1.2074)
        lj = {"O" : [0.1047, 2.9373]}

        super().__init__("O2", charges = 0.0, lj=lj, bonds=b1)
        return

#hydrogen https://pubs.acs.org/doi/10.1021/acs.jctc.0c01132 /!\: divide sig by 2**(1/6)
class Hydrogen(Specie):
    def __init__(self, Hset="std"):
        
        b1 = Bond("H", "H", kr=700, r0=0.7414)
        
        if Hset.lower() == "std":   # standard 12-6 set
            lj = {"H" : [0.0153, 2.5996]}
        elif Hset.lower() == "alt": # alternative 12-6 set
            lj = {"H" : [0.0145, 2.8001]}

        super().__init__("H2", charges = 0.0, lj=lj, bonds=b1)
        return

#nitrogen https://pubs.acs.org/doi/10.1021/acs.jctc.0c01132 /!\: divide sig by 2**(1/6)
class Nitrogen(Specie):
    def __init__(self):
        
        b1 = Bond("N", "N", kr=3190, r0=1.0977)
        lj = {"N" : [0.0797, 3.2197]}

        super().__init__("N2", charges = 0.0, lj=lj, bonds=b1)
        return
