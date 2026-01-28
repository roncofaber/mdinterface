#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 11:41:32 2025

@author: roncofaber
"""

def find_rings(graph, max_ring_size=8):
    """
    Find all rings in the molecular graph using NetworkX cycle detection.

    Parameters:
    max_ring_size (int): Maximum ring size to detect (default 8)

    Returns:
    list: List of rings, where each ring is a list of atom indices
    """
    import networkx as nx

    rings = []
    try:
        # Find simple cycles (rings) in the graph
        for cycle in nx.simple_cycles(graph):
            if len(cycle) <= max_ring_size:
                rings.append(sorted(cycle))
    except:
        # Fallback: use minimum cycle basis for undirected graphs
        try:
            cycle_basis = nx.minimum_cycle_basis(graph)
            for cycle in cycle_basis:
                if len(cycle) <= max_ring_size:
                    rings.append(sorted(cycle))
        except:
            # If all else fails, return empty list
            rings = []

    return rings

def get_rings_containing_atoms(graph, atom_indices, max_ring_size=8):
    """
    Find all rings that contain any of the specified atoms.

    Parameters:
    atom_indices (list): List of atom indices to check
    max_ring_size (int): Maximum ring size to detect

    Returns:
    list: List of rings containing any of the specified atoms
    """
    all_rings = find_rings(graph, max_ring_size)
    relevant_rings = []

    atom_set = set(atom_indices)
    for ring in all_rings:
        if atom_set.intersection(set(ring)):
            relevant_rings.append(ring)

    return relevant_rings