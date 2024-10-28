#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 11:30:04 2024

@author: roncofaber
"""

import numpy as np

#%%

def map_atoms(atoms):
    
    # Create an array to store the type ID of each atom
    atom_type_ids = []
    type_id = 0
    atoms_map = {}
    atoms_list = []
    
    for cc, atom in enumerate(atoms):
        if atom not in atoms_list:
            atom_type_ids.append(atom.label)
            atoms_list.append(atom)
            atoms_map[atom.label] = type_id
            type_id += 1
            
        else:
            idx = atoms_list.index(atom)
            atom_type_ids.append(atom.label)
            atoms_map[atom.label] = idx
    
    atom_type_ids = np.array(atom_type_ids)
    
    return atoms_list, atoms_map, atom_type_ids

def map_bonds(bonds):
    bond_type_ids = []
    type_id = 0
    bond_map = {}
    bonds_list = []
    
    for cc, bond in enumerate(bonds):
        if bond not in bonds_list:
            bond_type_ids.append(type_id)
            bonds_list.append(bond)
            bond_map[bond.symbols] = type_id
            type_id += 1
        else:
            idx = bonds_list.index(bond)
            bond_type_ids.append(idx)
            bond_map[bond.symbols] = idx
    
    bond_type_ids = np.array(bond_type_ids)
    return bonds_list, bond_map, bond_type_ids

def map_angles(angles):
    angle_type_ids = []
    type_id = 0
    angle_map = {}
    angles_list = []
    
    for cc, angle in enumerate(angles):
        if angle not in angles_list:
            angle_type_ids.append(type_id)
            angles_list.append(angle)
            angle_map[angle.symbols] = type_id
            type_id += 1
        else:
            idx = angles_list.index(angle)
            angle_type_ids.append(idx)
            angle_map[angle.symbols] = idx
    
    angle_type_ids = np.array(angle_type_ids)
    return angles_list, angle_map, angle_type_ids

def map_dihedrals(dihedrals):
    dihedral_type_ids = []
    type_id = 0
    dihedral_map = {}
    dihedrals_list = []
    
    for cc, dihedral in enumerate(dihedrals):
        if dihedral not in dihedrals_list:
            dihedral_type_ids.append(type_id)
            dihedrals_list.append(dihedral)
            dihedral_map[dihedral.symbols] = type_id
            type_id += 1
        else:
            idx = dihedrals_list.index(dihedral)
            dihedral_type_ids.append(idx)
            dihedral_map[dihedral.symbols] = idx
    
    dihedral_type_ids = np.array(dihedral_type_ids)
    return dihedrals_list, dihedral_map, dihedral_type_ids

def map_impropers(impropers):
    
    if impropers is None:
        return None, None
    
    improper_type_ids = []
    type_id = 0
    improper_map = {}
    impropers_list = []
    
    for cc, improper in enumerate(impropers):
        if improper not in impropers_list:
            improper_type_ids.append(type_id)
            impropers_list.append(improper)
            improper_map[improper.symbols] = type_id
            type_id += 1
        else:
            idx = impropers_list.index(improper)
            improper_type_ids.append(idx)
            improper_map[improper.symbols] = idx
    
    improper_type_ids = np.array(improper_type_ids)
    return impropers_list, improper_map, improper_type_ids
