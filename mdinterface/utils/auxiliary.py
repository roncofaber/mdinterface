#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 14:04:43 2024

@author: roncofaber
"""

import re
import collections
import numpy as np
import networkx as nx

import ase
from ase import neighborlist
#%%

def label_to_element(atostr, atomss):
    """
    Attempts to determine the chemical element symbol corresponding to a given 
    string and atomic mass.

    Args:
        atostr (str): A string potentially representing a chemical element.
        atomss (float): The approximate atomic mass of the element.

    Returns:
        str: The inferred chemical element symbol.

    Raises:
        ValueError: If the function cannot determine a valid element from the input.
    """


    new_label = re.sub(r'[^A-Za-z]', '', atostr).capitalize()    # Clean up input
    
    # load ase info
    atomic_masses = ase.data.atomic_masses
    elements = ase.data.chemical_symbols
    
    # initialize variables
    is_ready = False  # Flag to track if the element is found
    tried_last_resort = False 

    while not is_ready:
        try:
            try_atom = ase.Atom(new_label)  # Attempt to create an Atom object
            existent = True
        except:
            existent = False

        if existent and np.abs(try_atom.mass - atomss) < 1:  # Check mass match
            is_ready = True

        if not is_ready:
            new_label = new_label[:-1]  # Shorten the label for the next attempt

            if not new_label:  # If the label is empty, try a last-resort approach
                
                new_label = elements[np.argmin(np.abs(atomss - atomic_masses))]

                if tried_last_resort:
                    raise ValueError("{} is not a valid element".format(new_label))

                tried_last_resort = True

    return new_label


# return copy of input as list if not one
def as_list(inp):
    if inp is None:
        return []
    elif isinstance(inp, int) or isinstance(inp, np.int64):
        return [inp]
    elif isinstance(inp, collections.abc.Iterable) and not isinstance(inp, str): 
        # Handles lists, tuples, NumPy arrays, etc. (Excludes strings)
        return list(inp)  
    else:
        return [inp] # prone to error?
        # raise TypeError(f"Cannot convert type {type(inp)} to list")
        

def find_smallest_missing(data, start=0):
    """Finds the next smallest integer that is not in the list.
  
    This function efficiently finds the next smallest integer that is not present in the input list.
    It leverages sets for fast membership checks.
    
    Args:
        data: A list of integers.
  
    Returns:
        The next smallest integer that is not in the list.
      """

    data_set = set(data)  # Convert the list to a set for efficient membership checks
    smallest = start          # Start with the smallest possible positive integer
    while smallest in data_set:  # Check if 'smallest' is in the set
        smallest += 1              # If found, increment 'smallest'
    return smallest            # Return the first integer not found in the set 


def remove_inverted_tuples(list_of_tuples):
    seen = set()
    for i in range(len(list_of_tuples) - 1, -1, -1):  # Iterate backwards
        tup = list_of_tuples[i]
        reversed_tup = tup[::-1]
        if reversed_tup in seen:
            del list_of_tuples[i]  # Remove the tuple
        else:
            seen.add(tup)
            
# return list of indexes from mixed input of indexes and string (elements)
def atoms_to_indexes(system, symbols):

    # check if symbols is a list of strings
    if isinstance(symbols, str):
        if symbols == 'all':
            return list(range(len(system.get_chemical_symbols())))

    symbols = as_list(symbols)

    indexes = []
    for symbol in symbols:
        if not isinstance(symbol, str):
            indexes.append(symbol)
        else:
            for cc, atom in enumerate(system.get_chemical_symbols()):
                if atom == symbol:
                    indexes.append(cc)
    return np.unique(indexes).tolist()

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

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


def same_rev_check(list1, list2):
    if list(list1) == list(list2):
        return True
    elif list(list1) == list(reversed(list2)):
        return True
    return False
