#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Here you can find how to create a new specie

Created on Fri Mar 28 14:16:16 2025

@author: roncofaber
"""

from mdinterface import Specie

#%% Make a specie, and use LigParGen to estimate FF parameters

my_specie = Specie("CH3ONO", ligpargen=True) # molecule is in ASE database

#%% Re-calculate the charge on atoms using ab-initio pyscf calculations

q0 = my_specie.estimate_charges("resp", # charge scheme: "obabel", "ligpargen" or "resp"
                                calc_type = "RKS",
                                optimize  = True, # also run geo optimization of molecule (slower)
                                basis     = "6311++gss", 
                                xc        = "pbe",
                                assign    = True # assign charges to specie
                                )
