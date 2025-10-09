#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 16:03:43 2025

@author: roncofaber
"""

import numpy as np
import copy

# other stuff
import ase
from ase.data import atomic_numbers, covalent_radii

#%%

# snip the polymer into a smaller molecule, return as ase.Atoms
def make_snippet(polymer, center, Nmax, ending="H", preserve_rings=True):

    # Get initial atom selection based on distance
    ini_idxs = list(set(np.concatenate(polymer.find_relevant_distances(Nmax, centers=center))))

    # Reconstruct broken rings, if any
    if preserve_rings:
        # Find rings that contain any atoms in our initial selection
        relevant_rings = polymer._get_rings_containing_atoms(ini_idxs)

        # Add all atoms from relevant rings to our selection
        ring_atoms = set()
        for ring in relevant_rings:
            ring_atoms.update(ring)

        # Merge ring atoms with distance-based selection
        ini_idxs = list(set(ini_idxs) | ring_atoms)

    # Add all nodes that are connected by one connection to those nodes (and nothing else)
    # Iterate over the initial nodes
    con_idxs   = set()
    ter_idxs   = set()
    distances  = []
    atom_pairs = []
    for node in ini_idxs:
        # Get the neighbors of the current node
        neighbors = list(polymer.graph.neighbors(node))
        
        # Check each neighbor to see if it connects back to the initial nodes
        for neighbor in neighbors:
            
            # ignore if already in the list
            if neighbor in ini_idxs:
                continue
            
            # If it connects to exactly one node in the initial list, add it to connected_nodes
            if polymer.graph.degree(neighbor) == 1:
                con_idxs.add(neighbor)
            else:
                if neighbor not in ter_idxs:
                    ter_idxs.add(neighbor)
                    d1 = covalent_radii[atomic_numbers[polymer.graph.nodes[node]["element"]]]
                    d2 = covalent_radii[atomic_numbers[ending]]
                    distances.append(d1+d2)
                    atom_pairs.append([node, neighbor])

    # Combine the initial nodes with the newly found connected nodes
    con_idxs = list(set(ini_idxs).union(con_idxs))
    ter_idxs = list(ter_idxs)

    snippet_idxs = np.array(con_idxs + ter_idxs)
    
    # create chain and termination
    chain = polymer.atoms[con_idxs].copy()
    term  = polymer.atoms[ter_idxs].copy()
    term.set_chemical_symbols(len(term) * [ending])

    # make snippet        
    snippet = ase.Atoms(chain + term)
    
    # fix distances for ligpargen
    for cc, (a1, a2) in enumerate(atom_pairs):
        idx1 = int(np.where(snippet_idxs == a1)[0])
        idx2 = int(np.where(snippet_idxs == a2)[0])
        snippet.set_distance(idx1, idx2, distances[cc], fix=0)

    return snippet, snippet_idxs
 
def remap_snippet_topology(original_idxs, sn_atoms, sn_atypes, sn_bonds,
                           sn_angles, sn_dihedrals, sn_impropers, local_idxs):
    
    # Create filtered topology objects
    sn_bonds_out     = []
    sn_angles_out    = []
    sn_dihedrals_out = [] 
    sn_impropers_out = []
    
    # Prepare remapping IDs based on original IDs from snippet indices
    remap_ids = {}
    for cc, atype in enumerate(sn_atypes):
        remap_ids[atype.label] = str(original_idxs[cc])
    
    # Update bonds
    for bond in copy.deepcopy(sn_bonds):
        a1 = remap_ids[bond._a1]
        a2 = remap_ids[bond._a2]
        
        if all([ii in local_idxs for ii in [a1, a2]]):
            bond.update(a1=a1, a2=a2)
            sn_bonds_out.append(bond)
    
    # Update angles
    for angle in copy.deepcopy(sn_angles):
        a1 = remap_ids[angle._a1]
        a2 = remap_ids[angle._a2]
        a3 = remap_ids[angle._a3]
        if all([ii in local_idxs for ii in [a1, a2, a3]]):
            angle.update(a1=a1, a2=a2, a3=a3)
            sn_angles_out.append(angle)
    
    # Update dihedrals
    for dihedral in copy.deepcopy(sn_dihedrals):
        a1 = remap_ids[dihedral._a1]
        a2 = remap_ids[dihedral._a2]
        a3 = remap_ids[dihedral._a3]
        a4 = remap_ids[dihedral._a4]
        if all([ii in local_idxs for ii in [a1, a2, a3, a4]]):
            dihedral.update(a1=a1, a2=a2, a3=a3, a4=a4)
            sn_dihedrals_out.append(dihedral)
    
    # Update impropers
    for improper in copy.deepcopy(sn_impropers):
        a1 = remap_ids[improper._a1]
        a2 = remap_ids[improper._a2]
        a3 = remap_ids[improper._a3]
        a4 = remap_ids[improper._a4]
        if all([ii in local_idxs for ii in [a1, a2, a3, a4]]):
            # K, d, n = improper.values  # Assuming values consist of K, d, n parameters
            improper.update(a1=a1, a2=a2, a3=a3, a4=a4)
            sn_impropers_out.append(improper)
        
    return sn_bonds_out, sn_angles_out, sn_dihedrals_out, sn_impropers_out