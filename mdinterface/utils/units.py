#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 09:35:05 2024

@author: roncofaber
"""

from ase import units
from scipy import constants as const

#%% Physical constants, units and stuff

# assuming a volume element dV of 1 A^3
molar_to_number = units.mol*10e-28

# Boltzmann's constant J/K
kBoltzmann = units.kB/units.J

# elementary charge in C
elementary_charge = 1/units.C

# Temperature K
Temperature = 298.0

# kT in J
kT = kBoltzmann * Temperature

# eV in kT at 298K
eV_per_kT = kT/elementary_charge

# vacuum permittivity in F/m = C / (V m)
epsilon_0_e = const.epsilon_0

#vacuum permittivity in  q/(V Ang)
eps00 = (units.C/units.m )*epsilon_0_e

# Pascal in kT/A**3
Pascal = 1/(kT*units.m**3)
atmosphere = 101325*Pascal

au_per_ang3_to_g_per_cm3 = 1/(units.mol*1e-24)
