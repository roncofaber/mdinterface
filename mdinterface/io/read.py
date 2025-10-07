#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 18:33:06 2024

@author: roncofaber
"""

import numpy as np

import ase
import ase.visualize
import ase.io.lammpsdata
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral, Improper
from mdinterface.utils.auxiliary import mass2symbol

#%%

def read_lammps_data_file(filename, pbc=False, ato_start_idx=0):
    system = ase.io.lammpsdata.read_lammps_data(filename)
    
    if not pbc:
        system.set_pbc(False)
        system.set_cell(None)

    symbols = system.get_chemical_symbols()
    unique_symbols = set(symbols)
    
    pair_coeff     = []
    bond_coeff     = []
    angle_coeff    = []
    dihedral_coeff = []
    improper_coeff = []

    atoms     = []
    atosym    = []
    bonds     = []
    angles    = []
    dihedrals = []
    impropers = []

    section_map = {
        "masses"     :     "masses",
        'pair coeffs':     'pair_coeff',
        'bond coeffs':     'bond_coeff',
        'angle coeffs':    'angle_coeff',
        'dihedral coeffs': 'dihedral_coeff',
        'improper coeffs': 'improper_coeff',
        'atoms':            'atoms',
        'bonds':            'bonds',
        'angles':           'angles',
        'dihedrals':        'dihedrals',
        'impropers':        'impropers'
    }

    with open(filename, "r") as file:
        to_read = None

        for line in file:
            line = line.strip()

            if not line:
                continue

            section_key = line.lower()

            if section_key in section_map:
                to_read = section_map[section_key]
                continue

            line = line.split('#')[0].split()

            if to_read == "masses":
                mass = float(line[1])
                atosym.append(mass2symbol(mass, unique_symbols))
                
            elif to_read == "pair_coeff":
                eps = float(line[1])
                sig = float(line[2])
                pair_coeff.append([eps, sig])

            elif to_read == "bond_coeff":
                kr = float(line[1])
                r0 = float(line[2])
                bond_coeff.append([kr, r0])

            elif to_read == "angle_coeff":
                kr = float(line[1])
                theta0 = float(line[2])
                angle_coeff.append([kr, theta0])

            elif to_read == "dihedral_coeff":
                values = [float(ii) for ii in line[1:]]
                dihedral_coeff.append(values)

            elif to_read == "improper_coeff":
                values = [float(line[1]), int(line[2]), int(line[3])]
                improper_coeff.append(values)

            elif to_read == "atoms":
                atoidx = ato_start_idx + (int(line[0]) - 1)
                atotyp = int(line[2]) - 1
                eps, sig = pair_coeff[atotyp]
                symbol   = atosym[atotyp]
                atom = Atom(symbol, f"{symbol}_{str(atoidx).zfill(3)}", eps, sig)
                atoms.append(atom)

            elif to_read == "bonds":
                idx1 = int(line[2]) - 1
                idx2 = int(line[3]) - 1
                kr, r0 = bond_coeff[int(line[1]) - 1]
                bond = Bond(atoms[idx1].label, atoms[idx2].label, kr, r0)
                bonds.append(bond)

            elif to_read == "angles":
                idx1 = int(line[2]) - 1
                idx2 = int(line[3]) - 1
                idx3 = int(line[4]) - 1
                kr, theta0 = angle_coeff[int(line[1]) - 1]
                angle = Angle(atoms[idx1].label, atoms[idx2].label, atoms[idx3].label, kr, theta0)
                angles.append(angle)

            elif to_read == "dihedrals":
                idxs = [int(ii) - 1 for ii in line[2:]]
                values = dihedral_coeff[int(line[1]) - 1]
                dihedral = Dihedral(atoms[idxs[0]].label, atoms[idxs[1]].label,
                                    atoms[idxs[2]].label, atoms[idxs[3]].label, *values)
                dihedrals.append(dihedral)

            elif to_read == "impropers":
                idxs = [int(ii) - 1 for ii in line[2:]]
                values = improper_coeff[int(line[1]) - 1]
                improper = Improper(atoms[idxs[0]].label, atoms[idxs[1]].label,
                                    atoms[idxs[2]].label, atoms[idxs[3]].label,
                                    K=values[0], d=values[1], n=values[2])
                impropers.append(improper)

    system.new_array("stype", np.array(atoms))
    return system, atoms, bonds, angles, dihedrals, impropers
