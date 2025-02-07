#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:00:01 2025

@author: roncofaber
"""

# repo stuff
from .specie import Specie
from mdinterface.externals import run_ligpargen
from mdinterface.build.polymerize import build_polymer

# other stuff
import ase
import random
import numpy as np

#%%

class Polymer(Specie):

    def __init__(self, atoms=None, charges=None, atom_types=None, bonds=None,
                 angles=None, dihedrals=None, impropers=None, lj={}, cutoff=1.0,
                 name=None, lammps_data=None, fix_missing=False, chg_scaling=1.0,
                 pbc=False, ligpargen=False, tot_charge=0, nrep=1, start_end_idxs=None,
                 target_distance=1.600, substitute=None, refine_charges=False,
                 offset=False):
        
        # initialize polymer stuff
        self._snippet_cache = {}
        
        # polymerize
        polymer = build_polymer(atoms, substitute, nrep, start_end_idxs=start_end_idxs,
                      target_distance=target_distance)
        
        # Initialize the parent class
        super().__init__(polymer, charges, atom_types, bonds, angles, dihedrals,
                         impropers, lj, cutoff, name, lammps_data, fix_missing,
                         chg_scaling, pbc, ligpargen, tot_charge)
        
        if refine_charges:
            self.refine_charges(offset=offset)
        
        return
    
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
    
    # return list of elements adjacent to a connection point
    def _get_connection_elements(self):

        # Get elements where there is a connection
        centers = np.argwhere(self.atoms.arrays["is_connected"]).flatten()

        poi = []
        for center in centers:
            idxs = set(self.find_relevant_distances(1, centers=center).flatten())
            if not any([ii in poi for ii in idxs]):
                poi.append(center)
        return poi

    def _update_charges(self, center, Nmax, charges, ending="F"):
        """
        Updates the charges for the specified species.
    
        Parameters:
        specie (ase.Atoms): The species object containing the atoms.
        center (int): The center atom index.
        Nmax (int): The maximum number of neighbors to consider.
        charges (np.ndarray): The array of charges to be updated.
        ending (str): The chemical symbol for the terminal atoms.
        """
        
        ldxs = list(set(np.concatenate(self.find_relevant_distances(3, centers=center))))
        
        snippet, snippet_idxs = self.make_snippet(center, Nmax, ending=ending)
    
        # Check if snippet already exists in cache
        snippet_hash = ''.join(snippet.get_chemical_symbols())
        
        if snippet_hash in self._snippet_cache:
            cached_charges = self._snippet_cache[snippet_hash]
            mapping = [np.argwhere(snippet_idxs == ll)[0][0] for ll in ldxs]
            charges[ldxs] = cached_charges[mapping]
            return
        
        # snippet charge
        if "nominal_charge" in snippet.arrays:
            sn_charge = snippet.arrays["nominal_charge"].sum()
        else:
            sn_charge = None
            
        # run ligpargen
        output = run_ligpargen(snippet, charge=sn_charge)
        
        # get charges
        new_charges = output[0].get_initial_charges()
        
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
        charges = self.charges
        centers = self._get_connection_elements()
        
        # get charges at every point
        for center in centers:
            self._update_charges(center, Nmax, charges, ending=ending)
            
        # bring back to zero
        if offset:
            
            if "nominal_charge" in self.atoms.arrays:
                target = self.atoms.arrays["nominal_charge"].sum()
            else:
                target = 0
            
            charges -= ((charges.sum() - target) / len(charges))
            
        self.atoms.set_initial_charges(charges)
        return
