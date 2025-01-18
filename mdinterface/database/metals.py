#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 11:39:16 2025

@author: roncofaber
"""

import numpy as np

from mdinterface.core.specie import Specie

from ase.build import fcc111
#%%

# those are in kcal/mol and r0, so it should be divided by 2**1/6
metal_params = {
    "Ag" : [4.56, 2.955],
    "Al" : [4.02, 2.925],
    "Au" : [5.29, 2.951],
    "Cu" : [4.72, 2.616],
    "Ni" : [5.65, 2.552],
    "Pb" : [2.93, 3.565],
    "Pd" : [6.15, 2.819],
    "Pt" : [7.80, 2.845],
    "Ac" : [6.51, 3.843],
    "Ca" : [3.36, 4.025],
    "Ce" : [6.38, 3.734],
    "Es" : [2.88, 4.133],
    "Ir" : [9.20, 2.785],
    "Fe" : [6.00, 2.590],
    "Rh" : [7.84, 2.757],
    "Sr" : [3.40, 4.379],
    "Th" : [8.47, 3.683],
    "Yb" : [2.71, 3.942]
}
        
#metal parameters from https://pubs.acs.org/doi/full/10.1021/jp801931d
class Metal111(Specie):
    
    def __init__(self, metal, lj=None, **kwargs):
        
        slab = fcc111(metal, size=(1,2,3), orthogonal=True, periodic=True)
        
        if lj is None:
            lj = metal_params[metal]
            lj_copy = {metal: [lj[0], round(lj[1] / np.power(2, 1/6), 5)]}
        else:
            lj_copy = lj
        
        super().__init__(atoms=slab, lj=lj_copy, **kwargs)
        
        return
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.resname})"
