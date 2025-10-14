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
from mdinterface.build.snippets import make_snippet, remap_snippet_topology

# import random
import numpy as np
import copy

#%%

class Polymer(Specie):

    def __init__(self, monomers=None, charges=None, atom_types=None, bonds=None,
                 angles=None, dihedrals=None, impropers=None, lj={}, cutoff=1.0,
                 name=None, lammps_data=None, fix_missing=False, chg_scaling=1.0,
                 pbc=False, ligpargen=False, tot_charge=None, nrep=None,
                 sequence=None, refine_polymer=False, offset=False, ending="H"):
        
        # initialize polymer stuff
        self._snippet_cache = {}
        self._sequence = sequence
        
        # polymerize
        polymer = build_polymer(monomers, sequence=sequence, nrep=nrep)
        
        # Initialize the parent class with polymerized monomers
        super().__init__(atoms=polymer, charges=charges, atom_types=atom_types,
                         bonds=bonds, angles=angles, dihedrals=dihedrals,
                         impropers=impropers, lj=lj, cutoff=cutoff, name=name,
                         lammps_data=lammps_data, fix_missing=fix_missing,
                         chg_scaling=chg_scaling, pbc=pbc, ligpargen=ligpargen,
                         tot_charge=tot_charge)
        
        if refine_polymer:
            self.refine_polymer_topology(Nmax=12, offset=offset, ending=ending)
        
        return
    
    # return list of elements adjacent to a connection point
    def _get_connection_elements(self):

        # Get elements where there is a connection
        centers = np.argwhere(self.atoms.arrays["is_connected"]).flatten()

        pairs = []
        for center in centers:
            if any(center in pp for pp in pairs):
                continue
            distances = self.atoms.get_distances(center, centers)
            idx = centers[np.argsort(distances)[1]]
            pairs.append([center, idx])
            
        return pairs
    

    def _update_connection(self, center, Nmax, charges, ending="H"):

        # make a lil snippet
        snippet, snippet_idxs = make_snippet(self, center, Nmax, ending=ending)
        
        # get local indexes (within dihedral from center)
        ldxs = list(set(np.concatenate(self.find_relevant_distances(4, centers=center))))
        
        # find mapping between indexes
        mapping = [np.argwhere(snippet_idxs == ll)[0][0] for ll in ldxs]
    
        # Check if snippet already exists in cache
        snippet_hash = ''.join(snippet.get_chemical_symbols())
        if snippet_hash in self._snippet_cache:
            
            cached_snippet = copy.deepcopy(self._snippet_cache[snippet_hash])
            
            sn_atoms     = cached_snippet["sn_atoms"]
            sn_atypes    = cached_snippet["sn_atypes"]
            sn_bonds     = cached_snippet["sn_bonds"]
            sn_angles    = cached_snippet["sn_angles"]
            sn_dihedrals = cached_snippet["sn_dihedrals"]
            sn_impropers = cached_snippet["sn_impropers"]
            new_charges  = cached_snippet["new_charges"]
        
        # if not, ligpargen it
        else:

            # snippet charge
            if "nominal_charge" in snippet.arrays:
                sn_charge = snippet.arrays["nominal_charge"].sum()
            else:
                sn_charge = None
                
            # run ligpargen
            sn_atoms, sn_atypes, sn_bonds, sn_angles, sn_dihedrals, sn_impropers =\
                run_ligpargen(snippet, charge=sn_charge)
            
            # get charges
            new_charges = sn_atoms.get_initial_charges()
            
            # Store the new snippet and its charges in cache
            self._snippet_cache[snippet_hash] = {
                "new_charges"  : new_charges,
                "sn_atoms"     : sn_atoms,
                "sn_atypes"    : sn_atypes,
                "sn_bonds"     : sn_bonds,
                "sn_angles"    : sn_angles,
                "sn_dihedrals" : sn_dihedrals,
                "sn_impropers" : sn_impropers,
                "snippet_idxs" : snippet_idxs,
                "local_idxs"   : ldxs
                }
        
        # update topology of section
        original_idxs = self._sids[snippet_idxs]
        local_idxs = self._sids[ldxs]
        new_bonds, new_angles, new_dihedrals, new_impropers = remap_snippet_topology(
            original_idxs, sn_atoms, sn_atypes, sn_bonds, sn_angles, sn_dihedrals,
            sn_impropers, local_idxs)
        
        charges[ldxs] = new_charges[mapping]
        
        self._add_to_topology(bonds=new_bonds, angles=new_angles,
                              dihedrals=new_dihedrals, impropers=new_impropers)
        
        return
    
    # main driver that refines charges across the whole polymer
    def refine_polymer_topology(self, Nmax=12, offset=False, ending="H"):
        """
        Refines the charges for the specified species.

        Parameters:
        specie (ase.Atoms): The species object containing the atoms.
        Nmax (int): The maximum number of neighbors to consider.
        offset (bool): Whether to apply an offset to the charges.

        Returns:
        np.ndarray: The refined charges.
        """

        # Clean topology first to remove any invalid interactions from polymerization
        self._cleanup_topology()

        # get charges and connection elements
        charges = self.charges
        centers = np.array([int(ii[1]) for ii in self._get_connection_elements()])

        # get charges at every point
        for center in centers:
            self._update_connection(center, Nmax, charges, ending=ending)

        # bring back to zero
        if offset:

            if "nominal_charge" in self.atoms.arrays:
                target = self.atoms.arrays["nominal_charge"].sum()
            else:
                target = 0

            charges -= ((charges.sum() - target) / len(charges))

        self.atoms.set_initial_charges(charges)
        return
