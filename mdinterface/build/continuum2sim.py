#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 10:34:26 2025

@author: roncofaber
"""

from ase import units
import numpy as np
from scipy import integrate

#%%
# return a discrete concentration profile given a "continuum" conc. profile
# for a specie. Return positions of the bins where ion should be.
def discretize_concentration(specie, conc_profile, z_coords, volume):
    
    xsize, ysize, zsize = volume

    conv = units.mol/1e27
    area = xsize*ysize

    # integrate conc. prof. over area
    cumsum = integrate.cumulative_trapezoid(conc_profile, z_coords, initial=0, axis=0)
    cumsum_conv = cumsum*conv*area

    # generate bins (1 atom in each one) for K
    Nions = 0
    bins = []
    for ii, con in enumerate(cumsum_conv):
        if con > np.maximum(cumsum_conv.max()/100, Nions):
            Nions += 1
            bins.append(z_coords[ii])        
    
    # now find middle point of bins
    z_pos = []
    for ii, pos in enumerate(bins[:-1]):
        z_pos.append((bins[ii] + bins[ii+1])/2)
        
    # add one last ion
    if cumsum_conv.max() % 1 > 0.5:
        z_pos.append((z_coords.max()+bins[-1])/2)
        Nions += 1
    
    return z_pos
