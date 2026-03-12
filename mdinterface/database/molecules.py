#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pre-parameterised small-molecule species: water models and diatomic gases.

Water models available: SPC/E (``"ewald"``), CHARMM TIP3P (``"charmm"``),
and SPC/E (``"spce"``).  Diatomic species (O2, H2, N2) use parameters from
Lim et al., J. Chem. Theory Comput. 2021, 17, 821.
"""

from mdinterface.core.specie import Specie
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral, Improper

#%%
# solvent https://docs.lammps.org/Howto_tip3p.html (Ewald model)
class Water(Specie):
    """
    Water molecule with pre-parameterised force-field parameters.

    Parameters
    ----------
    model : str, default ``"ewald"``
        Water model to use:

        - ``"ewald"``  -- SPC/E with Ewald-compatible charges
          (q_O=-0.83, q_H=0.415)
        - ``"charmm"`` -- CHARMM TIP3P
          (q_O=-0.834, q_H=0.417)
        - ``"spce"``   -- SPC/E
          (q_O=-0.8476, q_H=0.4238)
    **kwargs
        Forwarded to :class:`~mdinterface.core.specie.Specie`.

    Examples
    --------
    ::

        from mdinterface.database import Water
        water = Water(model="ewald")
    """
    def __init__(self, model="ewald", **kwargs):
        
        if model.lower() == "ewald":
            b1 = Bond("O", "H", kr=450, r0=0.9572)
            a1 = Angle("H", "O", "H", kr=55, theta0=104.52)
            charges = [-0.83, 0.415, 0.415]
            lj = {"O": [0.102, 3.188], "H": [0.0, 1.0]}
        
        elif model.lower() == "charmm":
            b1 = Bond("O", "H", kr=450, r0=0.9572)
            a1 = Angle("H", "O", "H", kr=55, theta0=104.52)
            charges = [-0.834, 0.417, 0.417]
            lj = {"O": [0.1521, 3.1507], "H": [0.0460, 0.4]}
            
        elif model.lower() == "spce":
            b1 = Bond("O", "H", kr=1, r0=1.0)
            a1 = Angle("H", "O", "H", kr=1, theta0=109.47)
            charges = [-0.8476, 0.4238, 0.4238]
            lj = {"O": [0.1553, 3.1660], "H": [0.0, 0.0]}

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
