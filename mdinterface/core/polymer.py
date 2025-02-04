#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:00:01 2025

@author: roncofaber
"""

# repo stuff
from .specie import Specie
from mdinterface.externals import run_ligpargen

# other stuff
import ase
import random
import numpy as np

#%%

class Polymer(Specie):

    def __init__(self, atoms=None, charges=None, atom_types=None, bonds=None,
                 angles=None, dihedrals=None, impropers=None, lj={}, cutoff=1.0,
                 name=None, lammps_data=None, fix_missing=False, chg_scaling=1.0,
                 pbc=False):
        
        # Initialize the parent class
        super(Polymer, self).__init__(atoms, charges, atom_types, bonds, angles, dihedrals,
                                      impropers, lj, cutoff, name, lammps_data, fix_missing,
                                      chg_scaling, pbc)
        
        # initialize polymer stuff
        self._snippet_cache = {}
        
        return
    
    # return list of elements adjacent to a connection point
    def get_connection_elements(self):

        # Get elements where there is a connection
        centers = np.argwhere(self.atoms.arrays["is_connected"]).flatten()

        poi = []
        for center in centers:
            idxs = set(self.find_relevant_distances(1, centers=center).flatten())
            if not any([ii in poi for ii in idxs]):
                poi.append(center)
        return poi
    
    # snip the polymer into a smaller molecule, return as ase.Atoms
    def make_snippet(self, centers, Nmax, ending="F"):
        
        idxs = list(set(np.concatenate(self.find_relevant_distances(Nmax, centers=centers))))
        edxs = list(set(np.concatenate(self.find_relevant_distances(Nmax+1, centers=centers, Nmin=Nmax))) - set([centers]))
        
        snippet_idxs = np.array(idxs + edxs)

        
        chain = self.atoms[idxs].copy()
        term  = self.atoms[edxs].copy()
        
        term.set_chemical_symbols(len(term) * [ending])
        snippet = ase.Atoms(chain + term)
        
        return snippet, snippet_idxs

    def update_charges(self, center, Nmax, charges, ending="F"):
        """
        Updates the charges for the specified species.
    
        Parameters:
        specie (ase.Atoms): The species object containing the atoms.
        center (int): The center atom index.
        Nmax (int): The maximum number of neighbors to consider.
        charges (np.ndarray): The array of charges to be updated.
        ending (str): The chemical symbol for the terminal atoms.
        """
        
        ldxs = list(set(np.concatenate(self.find_relevant_distances(4, centers=center))))
        
        snippet, snippet_idxs = self.make_snippet(center, Nmax, ending=ending)
    
        # Check if snippet already exists in cache
        snippet_hash = ''.join(snippet.get_chemical_symbols())
        
        if snippet_hash in self._snippet_cache:
            cached_charges = self._snippet_cache[snippet_hash]
            mapping = [np.argwhere(snippet_idxs == ll)[0][0] for ll in ldxs]
            charges[ldxs] = cached_charges[mapping]
            return
    
        # run ligpargen
        new_charges = run_ligpargen(snippet)
        
        # find mapping
        mapping = [np.argwhere(snippet_idxs == ll)[0][0] for ll in ldxs]
        
        # update charges
        charges[ldxs] = new_charges[mapping]
    
        # Store the new snippet and its charges in cache
        self._snippet_cache[snippet_hash] = new_charges
        
        return
    
    # main driver that refines charges across the whole polymer
    def refine_charges(self, Nmax=12, offset=False, ending="F"):
        """
        Refines the charges for the specified species.
    
        Parameters:
        specie (ase.Atoms): The species object containing the atoms.
        Nmax (int): The maximum number of neighbors to consider.
        offset (bool): Whether to apply an offset to the charges.
    
        Returns:
        np.ndarray: The refined charges.
        """
        
        # get charges and connection elements
        charges = self.atoms.get_initial_charges()
        centers = self.get_connection_elements()
        
        # get charges at every point
        for center in centers:
            self.update_charges(center, Nmax, charges, ending=ending)
            
        # bring back to zero
        if offset:
            charges -= (charges.sum() / len(charges))
    
        return charges
