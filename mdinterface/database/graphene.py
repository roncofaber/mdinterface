#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Graphene species with OPLS-AA force-field parameters.

Bond, angle, dihedral, and improper parameters from
Jorgensen & Severance, J. Am. Chem. Soc. 1990;
LJ parameters (epsilon=0.07 kcal/mol, sigma=3.55 Å) from
Koenig et al., Adv. Mater. 2018.
"""

import numpy as np

# internal modules
from mdinterface.core.specie import Specie
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral, Improper

import ase
#%%

# graphene https://onlinelibrary.wiley.com/doi/10.1002/adma.201705791
class Graphene(Specie):
    """
    Single graphene sheet with OPLS-AA force-field parameters.

    The unit cell is built with :func:`ase.build.graphene` and the
    interlayer spacing is set to 3.35 Å.  Bond, angle, dihedral, and
    improper parameters follow the OPLS-AA graphene parametrisation.

    Parameters
    ----------
    **kwargs
        Forwarded to :class:`~mdinterface.core.specie.Specie`.

    Examples
    --------
    ::

        from mdinterface.database import Graphene
        grap = Graphene()
    """
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




