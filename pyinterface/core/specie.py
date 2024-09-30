#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 13:58:43 2024

@author: roncofaber
"""

import ase
import ase.build
import ase.visualize
from ase.geometry.analysis import Analysis
from ase.neighborlist import NeighborList, natural_cutoffs
import MDAnalysis as mda

from pyinterface.utils.auxiliary import as_list, find_smallest_missing,\
    remove_inverted_tuples, atoms_to_indexes
    
from pyinterface.core.topology import Atom

import copy
import numpy as np
#%%

class Specie(object):
    
    def __init__(self, atoms, charges = None, bonds=None, angles=None,
                 dihedrals=None, impropers=None, lj={}, cutoff=1.0):
        
        # read atoms
        atoms = self._read_atoms(atoms, charges)
        
        # set up atoms
        self.set_atoms(atoms, cutoff)
        
        # set up internal topology attributes
        self._btype = as_list(bonds)
        self._atype = as_list(angles)
        self._dtype = as_list(dihedrals)
        self._itype = as_list(impropers)
        self._stype = self._setup_atom_types(lj)
        
        # initialize topology info
        self._initialize_topology()
        
        return
    
    # read atoms to return ase.Atoms
    @staticmethod
    def _read_atoms(atoms, charges):
        
        # initialize atoms obj
        if isinstance(atoms, str):
            try:
                atoms = ase.io.read(atoms)
            except:
                atoms = ase.build.molecule(atoms)
        elif isinstance(atoms, ase.Atoms):
            atoms = atoms.copy()
        
        if charges is not None:
            charges = as_list(charges)
            if len(charges) == 1:
                charges = len(atoms)*charges
            atoms.set_initial_charges(charges)

        return atoms
    
    def _setup_atom_types(self, lj):
        
        atoms = self.atoms
    
        formula = atoms.get_chemical_formula()
        
        stype   = []
        for atom_symbol in np.unique(atoms.get_chemical_symbols()):

            if atom_symbol in lj:
                eps, sig = lj[atom_symbol]
            else:
                eps, sig = 0,0
            
            atom = Atom(atom_symbol, label="{}_{}".format(atom_symbol, formula),
                        eps=eps, sig=sig)
            stype.append(atom)
        
        return as_list(stype)
    
    def set_atoms(self, atoms, cutoff=1.0):
        
        self._atoms = atoms
        
        atoms_cutoff = cutoff*np.array(natural_cutoffs(atoms))
        
        nl = NeighborList(atoms_cutoff)#, bothways=True)
        nl.update(atoms)
        self._ana  = Analysis(self.atoms, nl=nl)
        
        return
    
    def _initialize_topology(self):
        
        formula = self.atoms.get_chemical_formula()
        
        for attributes in ["_btype", "_atype", "_dtype", "_itype",
                           "_stype"]:
            attr_type = []
            for attr in self.__getattribute__(attributes):
                attr.set_formula(formula)
                if attr.id is None or attr.id in attr_type:
                    idx = find_smallest_missing(attr_type, start=1)
                    attr.set_id(idx)
                    attr_type.append(idx)
                
        return
    

    # covnert to mda.Universe
    def to_universe(self, charges=True, layered=False):
        
        # empty top object
        top = mda.core.topology.Topology(n_atoms=len(self.atoms))
        
        # empty universe
        uni = mda.Universe(top, self.atoms.get_positions())
        
        # add some stuff
        uni.add_TopologyAttr("masses", self.atoms.get_masses())
        uni.add_TopologyAttr("resnames", [self.atoms.get_chemical_formula()])
        
        # generate type
        types, indexes = self.get_atom_types(return_index=True)
        uni.add_TopologyAttr("types", types)
        uni.add_TopologyAttr("type_index", indexes)
        
        # populate with bonds and angles
        for att in ["bonds", "angles", "dihedrals", "impropers"]:
            attribute, types = self.__getattribute__(att)
            uni._add_topology_objects(att, attribute, types=types)
        
        # add charges
        if charges:
            uni.add_TopologyAttr("charges", self.atoms.get_initial_charges())
            
        # layer it nicely
        if layered:
            layer_idxs, __ = ase.geometry.get_layers(self.atoms, [0,0,1],
                                                      tolerance=0.01)
            groups = []
            for idx in np.unique(layer_idxs):
                idxs = np.where(layer_idxs == idx)[0]
                
                groups.append(uni.atoms[idxs])
                
            uni = mda.Merge(*groups)
            
            # uni.residues[0].resids = layer_idxs
        
        # if has cell info, pass them along
        if self.atoms.get_cell():
            uni.dimensions = self.atoms.cell.cellpar()
                
        return uni
    
    def repeat(self, rep, make_cubic=False):
        
        atoms = self.atoms.repeat(rep)
        
        if make_cubic: #TODO: this can be dangerous if the cell
            
            xsize = [1,0,0]@atoms.cell@[1,0,0]
            ysize = [0,1,0]@atoms.cell@[0,1,0]
            zsize = [0,0,1]@atoms.cell@[0,0,1]
            
            atoms.set_cell([xsize, ysize, zsize, 90, 90, 90])
            atoms.wrap()
        
        self.set_atoms(atoms)
        
        return
    
    @property
    def atoms(self):
        return self._atoms
    
    @property
    def ana(self):
        return self._ana
    
    @property
    def bonds(self):
        
        bond_list = []
        bond_type = []
        for bond in self._btype:
            new_bonds = self.ana.get_bonds(*bond.symbols)[0]
            bond_list.extend(new_bonds)
            bond_type.extend(len(new_bonds)*[bond.id])
            
        return [bond_list, bond_type]
    
    @property
    def angles(self):
        
        ang_list = []
        ang_type = []
        for angle in self._atype:
            new_angles = self.ana.get_angles(*angle.symbols)[0]
            ang_list.extend(new_angles)
            ang_type.extend(len(new_angles)*[angle.id])
                
        return [ang_list, ang_type]
    
    @property
    def dihedrals(self):
        
        dih_list = []
        dih_type = []
        for dihedral in self._dtype:
            new_dihs = self.ana.get_dihedrals(*dihedral.symbols)[0]
            remove_inverted_tuples(new_dihs)
            dih_list.extend(new_dihs)
            dih_type.extend(len(new_dihs)*[dihedral.id])
        
        return [dih_list, dih_type]
    
    @property
    def impropers(self):
        
        imp_list = []
        imp_type = []
        
        all_bonds = self.ana.all_bonds[0]
        
        for improper in self._itype:
            idxs = atoms_to_indexes(self.atoms, improper.symbols)
            
            new_imps = [(cc, *ii) for cc, ii in enumerate(all_bonds)
                       if len(ii) == 3 and cc in idxs]
            
            imp_list.extend(new_imps)
            imp_type.extend(len(new_imps)*[improper.id])
            
        return [imp_list, imp_type]
    
    def copy(self):
        return copy.deepcopy(self)
    
    def get_atom_types(self, return_index=False):
        
        atom_types = []
        atom_type_idx = []
        for atom in self.atoms:
            for stype in self._stype:
                if stype.symbol == atom.symbol:
                    atom_types.append(stype.label)
                    atom_type_idx.append(stype.id)
                    break
        
        if return_index:
            return atom_types, np.array(atom_type_idx, dtype=int)
        
        return atom_types
    
    def estimate_sphere_radius(self):
        """
        Estimate the radius of the sphere containing the given points.

        Parameters:
        points (numpy.ndarray): A 2D array of shape (n, 3) where n is the number of points.

        Returns:
        float: The estimated radius of the sphere.
        """
        
        points = self.atoms.get_positions()
        
        # Calculate the centroid of the points
        centroid = np.mean(points, axis=0)

        # Calculate the distances from the centroid to each point
        distances = np.linalg.norm(points - centroid, axis=1)

        # The radius of the sphere is the maximum distance from the centroid to any point
        radius = np.max(distances)

        return radius
    
    
    def view(self):
        ase.visualize.view(self.atoms)
        return
