#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 11:30:04 2024

@author: roncofaber
"""

import numpy as np
import mdinterface.utils.auxiliary as aux
from mdinterface.core.topology import Bond, Angle, Dihedral, Improper

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

#%%

def find_missing_bonds(nas):
    tmp_bonds, _ = nas.bonds
    all_bonds = aux.find_unique_paths_of_length(nas.graph, 1)

    # Convert tmp_bonds to a set of tuples for efficient membership checking
    tmp_bonds_set = set(tuple(bond) for bond in tmp_bonds)
    tmp_bonds_set.update(tuple(reversed(bond)) for bond in tmp_bonds)

    missing_bonds = []
    for tmp_bond in all_bonds:
        tmp_bond_tuple = tuple(tmp_bond)
        tmp_bond_tuple_rev = tuple(reversed(tmp_bond_tuple))
        if tmp_bond_tuple not in tmp_bonds_set and tmp_bond_tuple_rev not in tmp_bonds_set:
            missing_bonds.append(tuple(nas._sids[ii] for ii in tmp_bond_tuple))

    return missing_bonds

def find_missing_angles(nas):
    tmp_angles, _ = nas.angles
    all_angles = aux.find_unique_paths_of_length(nas.graph, 2)

    # Convert tmp_angles to a set of tuples for efficient membership checking
    tmp_angles_set = set(tuple(angle) for angle in tmp_angles)
    tmp_angles_set.update(tuple(reversed(angle)) for angle in tmp_angles)

    missing = []
    for tmp_ang in all_angles:
        tmp_ang_tuple = tuple(tmp_ang)
        tmp_ang_tuple_rev = tuple(reversed(tmp_ang_tuple))
        if tmp_ang_tuple not in tmp_angles_set and tmp_ang_tuple_rev not in tmp_angles_set:
            missing.append(tuple(nas._sids[ii] for ii in tmp_ang_tuple))
    return missing

def find_missing_dihedrals(nas):
    tmp_dihedrals, _ = nas.dihedrals
    all_dihedrals = aux.find_unique_paths_of_length(nas.graph, 3)

    # Convert tmp_dihedrals to a set of tuples for efficient membership checking
    tmp_dihedrals_set = set(tuple(dihedral) for dihedral in tmp_dihedrals)
    tmp_dihedrals_set.update(tuple(reversed(dihedral)) for dihedral in tmp_dihedrals)

    missing_dihedrals = []
    for tmp_dihedral in all_dihedrals:
        tmp_dihedral_tuple = tuple(tmp_dihedral)
        tmp_dihedral_tuple_rev = tuple(reversed(tmp_dihedral_tuple))
        if tmp_dihedral_tuple not in tmp_dihedrals_set and tmp_dihedral_tuple_rev not in tmp_dihedrals_set:
            missing_dihedrals.append(tuple(nas._sids[ii] for ii in tmp_dihedral_tuple))

    return missing_dihedrals

def find_missing_impropers(nas):
    tmp_impropers, _ = nas.impropers
    all_impropers = aux.find_unique_paths_of_length(nas.graph, 3)  # Assuming path length 3 for impropers

    # Convert tmp_impropers to a set of tuples for efficient membership checking
    tmp_impropers_set = set(tuple(improper) for improper in tmp_impropers)
    tmp_impropers_set.update(tuple(reversed(improper)) for improper in tmp_impropers)

    missing_impropers = []
    for tmp_improper in all_impropers:
        tmp_improper_tuple = tuple(tmp_improper)
        tmp_improper_tuple_rev = tuple(reversed(tmp_improper_tuple))
        if tmp_improper_tuple not in tmp_impropers_set and tmp_improper_tuple_rev not in tmp_impropers_set:
            missing_impropers.append(tuple(nas._sids[ii] for ii in tmp_improper_tuple))

    return missing_impropers

#%%

def generate_missing_interactions(nas, interaction_type):
    mss_interactions = nas.suggest_missing_interactions(interaction_type)

    new_interactions = []
    interaction_type_map = {}

    # Determine the appropriate attributes based on the interaction type
    if interaction_type == "bonds":
        interaction_list = nas._btype
        num_atoms = 2
    elif interaction_type == "angles":
        interaction_list = nas._atype
        num_atoms = 3
    elif interaction_type == "dihedrals":
        interaction_list = nas._dtype
        num_atoms = 4
    elif interaction_type == "impropers":
        interaction_list = nas._itype
        num_atoms = 4
    else:
        raise ValueError("Invalid interaction type")

    # Create a mapping of interaction types for faster lookup
    for itype in interaction_list:
        ctypes = tuple(nas._smap[ii] for ii in itype.symbols)
        interaction_type_map[ctypes] = itype
        interaction_type_map[ctypes[::-1]] = itype  # Add the reversed tuple as well

    for mss_interaction in mss_interactions:
        found_interaction = False
        stypes = tuple(nas._smap[ii] for ii in mss_interaction)
        stypes_rev = stypes[::-1]
        
        # topology attribute already existing in map
        if stypes in interaction_type_map or stypes_rev in interaction_type_map:
            if found_interaction:
                raise ValueError(f"More than one possible {interaction_type[:-1]} found, abort!")
                
            if stypes in interaction_type_map:
                ninteraction = interaction_type_map[stypes].copy()
                
            else:
                ninteraction = interaction_type_map[stypes_rev].copy()
                mss_interaction = mss_interaction[::-1]
            
            a1, a2, *rest = mss_interaction
            
            ninteraction._a1 = a1
            ninteraction._a2 = a2
            
            if num_atoms > 2:
                ninteraction._a3 = rest[0]
            if num_atoms > 3:
                ninteraction._a4 = rest[1]
            new_interactions.append(ninteraction)
            found_interaction = True
        
        # new topology attribute
        else:
            if found_interaction:
                raise ValueError(f"More than one possible {interaction_type[:-1]} found, abort!")
            
            a1, a2, *rest = mss_interaction
            
            if num_atoms == 2:
                ninteraction = Bond(a1, a2)
            elif num_atoms == 3:
                ninteraction = Angle(a1, a2, rest[0])
            elif num_atoms == 4:
                ninteraction = Dihedral(a1, a2, rest[0], rest[1])
                
            new_interactions.append(ninteraction)
            found_interaction = True
            

    return new_interactions
