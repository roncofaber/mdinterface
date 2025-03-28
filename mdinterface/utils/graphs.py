#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:39:52 2025

@author: roncofaber
"""

# repo
from .auxiliary import atoms_to_indexes, as_list

# not repo
import collections
import numpy as np
from ase import neighborlist

# networking
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

#%%
# bunch of functions to work with NetworkX graphs

def get_nth_neighbors(graph, start_node, n):
    visited = set()
    queue = collections.deque([(start_node, 0)])
    neighbors = []

    while queue:
        current_node, depth = queue.popleft()
        if depth > n:
            break
        if current_node not in visited:
            visited.add(current_node)
            if depth > 0:  # Exclude the start_node itself
                neighbors.append(current_node)
            for neighbor in graph.neighbors(current_node):
                if neighbor not in visited:
                    queue.append((neighbor, depth + 1))

    return neighbors

def molecule_to_graph(molecule, cutoff_scale=1.0):
    
    # Generate cutoff
    cutOff = cutoff_scale*np.array(neighborlist.natural_cutoffs(molecule))
    ignore_atoms = ""  # Assuming ignore_atoms is an empty string
    cutOff[atoms_to_indexes(molecule, ignore_atoms)] = 0

    # Calculate neighbor list
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
    neighborList.update(molecule)

    # Check if it's a polymer
    is_polymer = "is_connected" in molecule.arrays
    
    # Generate graph
    G = nx.Graph()
    G.add_nodes_from(list(range(len(molecule))))
    
    # Iterate through neighbors to add edges
    for atom_index, bonded_atoms in enumerate(neighborList.nl.neighbors):
        for neighbor_index in bonded_atoms:
            if is_polymer:
                # Skip if the atoms belong to different monomers and are not connected
                if molecule.arrays["mon_id"][atom_index] != molecule.arrays["mon_id"][neighbor_index]:
                    if not (molecule.arrays["is_connected"][atom_index] and molecule.arrays["is_connected"][neighbor_index]):
                        continue
            # Add edge between the atoms
            G.add_edge(atom_index, neighbor_index)
    
    return G

def find_atom_types(molecule, max_depth=1):
    
    G = molecule_to_graph(molecule)
    
    # Get chemical symbols
    symbols = molecule.get_chemical_symbols()
    
    # Create a dictionary to store the unique atom types and their IDs
    atom_types = {}
    type_id = 0

    # Create an array to store the type ID of each atom
    atom_type_ids = np.zeros(len(molecule), dtype=int)

    # Iterate over each node in the graph
    for node in G.nodes:
        # Get the element of the current atom
        element = symbols[node]
        
        # Get the elements of the nth nearest neighboring atoms
        nth_neighbors = get_nth_neighbors(G, node, max_depth)
        neighbor_elements = sorted([symbols[neighbor] for neighbor in nth_neighbors])
        
        # Create a unique identifier for the atom type
        atom_type = (element, tuple(neighbor_elements))
        
        # Assign an ID to the atom type if it is not already in the dictionary
        if atom_type not in atom_types:
            atom_types[atom_type] = type_id
            type_id += 1
        
        # Store the type ID in the array
        atom_type_ids[node] = atom_types[atom_type]

    return atom_type_ids, {v: k for k, v in atom_types.items()}

# similar to find atoms types but works for whole graph. Keep an eye on this #TODO
def find_equivalent_atoms(G):

    # Step 1: Find automorphisms
    # Use the GraphMatcher to find all automorphisms of the graph
    GM = GraphMatcher(G, G)
    automorphisms = list(GM.isomorphisms_iter())
    
    # Step 2: Identify topologically equivalent atoms
    # Create a dictionary to store equivalence classes
    equivalence_classes = {}
    
    for i, atom in enumerate(G.nodes()):
        equivalence_classes[atom] = {atom}
    
    for automorphism in automorphisms:
        for atom, mapped_atom in automorphism.items():
            equivalence_classes[atom].add(mapped_atom)
    
    # Convert sets to sorted lists for easier interpretation
    for atom in equivalence_classes:
        equivalence_classes[atom] = sorted(equivalence_classes[atom])

    # Step 3: Assign unique labels to each equivalence class
    class_labels = {}
    label = 0
    for atom, eq_class in equivalence_classes.items():
        eq_class = frozenset(eq_class)
        if eq_class not in class_labels:
            class_labels[eq_class] = label
            label += 1
    
    # Step 4: Create the list with equivalence class labels
    equivalence_list = [0] * len(G.nodes())
    for atom, eq_class in equivalence_classes.items():
        eq_class = frozenset(eq_class)
        equivalence_list[atom] = class_labels[eq_class]
            
    return equivalence_classes, np.array(equivalence_list)

def find_unique_paths_of_length(graph, length):
    def dfs(current_node, current_path):
        if len(current_path) == length + 1:
            # Check if the reverse of the path already exists
            if tuple(current_path[::-1]) not in paths:
                paths.add(tuple(current_path))
            return
        for neighbor in graph.neighbors(current_node):
            if neighbor not in current_path:  # Avoid cycles
                dfs(neighbor, current_path + [neighbor])

    paths = set()
    for node in graph.nodes:
        dfs(node, [node])
    
    # Convert set of tuples back to list of lists
    paths = [list(path) for path in paths]
    
    # Sort paths lexicographically
    paths.sort()
    
    return np.array(paths, dtype=int)

def find_improper_idxs(graph):
    nodes = [node for node, degree in dict(graph.degree()).items() if degree == 3]

    # Find all nodes connected to those nodes
    improper_idxs = []
    for node in nodes:
        improper_idxs.append(sorted([node,* graph.neighbors(node)]))
        
    return improper_idxs

def find_relevant_distances(graph, Nmax, Nmin=0, centers=None, Ninv=0):
    
    # get list of relevant nodes
    if centers is None:
        relevant_nodes = graph.nodes()
    else:
        relevant_nodes = as_list(centers)
    
    # Set to store unique pairs
    unique_pairs = set()

    # Iterate over all nodes in the graph
    for node1 in relevant_nodes:
        # Get the shortest path lengths from node node to all other reachable nodes
        shortest_paths = nx.single_source_shortest_path_length(graph, node1)
        
        longest_path = max([dist for _, dist in shortest_paths.items()])
        
        # Collect pairs where the distance is within N but above Nmin
        for node2, distance in shortest_paths.items():
            if distance <= Nmax and distance > Nmin:
                # Use tuple (min(node, n), max(node, n)) to avoid duplicates
                pair = (min(node1, node2), max(node1, node2))
                unique_pairs.add(pair)
            elif Ninv > longest_path - distance:
                pair = (min(node1, node2), max(node1, node2))
                unique_pairs.add(pair)

    # Convert set to list
    unique_pairs_list = list(unique_pairs)
    unique_pairs_list.sort()
    
    return np.array(unique_pairs_list, dtype=int)
