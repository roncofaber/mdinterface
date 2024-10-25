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
import pyinterface.utils.auxiliary as aux
from pyinterface.io.read import read_lammps_data_file

import copy
import numpy as np
#%%

class Specie(object):
    
    def __init__(self, atoms=None, charges=None, bonds=None, angles=None,
                 dihedrals=None, impropers=None, lj={}, atom_types=None, cutoff=1.0,
                 name=None, lammps_data=None):
        
        if lammps_data is not None:
            atoms, atom_types, bonds, angles, dihedrals = read_lammps_data_file(lammps_data)
            
        # read atoms
        atoms = self._read_atoms(atoms, charges)
        
        # set name
        if name is None:
            name = atoms.get_chemical_formula()
        if len(name) > 4:
            print("ATTENTION: resname for Specie could be misleading")
        self.resname = name
        
        # set up atoms
        self.set_atoms(atoms, cutoff)
        
        # set up internal topology attributes
        self._btype = as_list(bonds)
        self._atype = as_list(angles)
        self._dtype = as_list(dihedrals)
        self._itype = as_list(impropers)
        self._stype = self._setup_atom_types(lj, atom_types)
        
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
    
    # function to setup atom types
    def _setup_atom_types(self, lj, atom_types):
        
        if atom_types is None:
            # use function to retrieve IDs
            type_ids, types_map = aux.find_atom_types(self.atoms, max_depth=1)
            
            stype = []
            for atom_type in types_map:
                
                atom_symbol = types_map[atom_type][0]
                atom_neighs = "".join(types_map[atom_type][1])
                
                label = "{}_{}".format(atom_symbol, atom_neighs)
                types_map[atom_type] = label
                
                if label in lj:
                    eps, sig = lj[label]
                elif atom_type in lj:
                    eps, sig = lj[atom_type]
                elif atom_symbol in lj:
                    eps, sig = lj[atom_symbol]
                else:
                    eps, sig = None,None
                
                atom = Atom(atom_symbol, label=label,
                            eps=eps, sig=sig)
                stype.append(atom)
                
        else:
            stype = atom_types
            type_ids = list(range(len(self.atoms)))
            types_map = {}
            for cc, atom_type in enumerate(atom_types):
                types_map[cc] = atom_type.label
        
        self._atom_types = types_map
        # add tag to atoms
        self.atoms.set_tags(type_ids)
        
        return as_list(stype)
    
    def set_atoms(self, atoms, cutoff=1.0):
        
        self._atoms = atoms
        
        self._graph = aux.molecule_to_graph(atoms)
        
        return
    
    def _initialize_topology(self):
        formula = self.atoms.get_chemical_formula()
        
        for attributes in ["_btype", "_atype", "_dtype", "_itype", "_stype"]:
            attr_type = []
            for attr in self.__getattribute__(attributes):
                attr.set_formula(formula)
                attr.set_resname(self.resname)
                if attr.id is None or attr.id in attr_type:
                    idx = find_smallest_missing(attr_type, start=1)
                    attr.set_id(idx)
                else:
                    idx = attr.id
                attr_type.append(idx)
        
        return

    
    def add_topology(self, attribute):
        # Determine the type of the attribute
        attribute_type = attribute.__class__.__name__
        
        # Get the corresponding attribute list
        if attribute_type == "Bond":
            attribute_list = self._btype
        elif attribute_type == "Angle":
            attribute_list = self._atype
        elif attribute_type == "Dihedral":
            attribute_list = self._dtype
        elif attribute_type == "Improper":
            attribute_list = self._itype
        elif attribute_type == "Atom":
            attribute_list = self._stype
        else:
            raise ValueError("Invalid topology attribute type.")
        
        # Check if the attribute already exists based on symbols
        attribute_symbols = attribute.symbols
        for existing_attribute in attribute_list:
            if existing_attribute.symbols == attribute_symbols:
                print(f"{attribute_type} with symbols {attribute_symbols} already exists and will not be added again.")
                return
        
        # Add the attribute to the list
        attribute_list.append(attribute)
        
        # Reinitialize topology to ensure consistency
        self._initialize_topology()
        
        return


    # covnert to mda.Universe
    def to_universe(self, charges=True, layered=False):
        
        # empty top object
        top = mda.core.topology.Topology(n_atoms=len(self.atoms))
        
        # empty universe
        uni = mda.Universe(top, self.atoms.get_positions())
        
        # add some stuff
        uni.add_TopologyAttr("masses", self.atoms.get_masses())
        uni.add_TopologyAttr("resnames", [self.resname])
        
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
    def graph(self):
        return self._graph
        
    @property
    def bonds(self):
        # Find all unique paths of length 1 (which corresponds to bonds)
        bonds = aux.find_unique_paths_of_length(self.graph, 1)
        tags = self.atoms.get_tags()
        
        bond_list = []
        bond_type = []
        
        for i, j in bonds:
            a1 = self._atom_types[tags[i]]
            a2 = self._atom_types[tags[j]]
            
            s1 = a1.split("_")[0]
            s2 = a2.split("_")[0]
            
            atom_pair_bond = None
            symb_pair_bond = None
            
            for bond in self._btype:
                bond_symbols = bond.symbols
                atom_pair = (a1, a2)
                symbol_pair = (s1, s2)
                
                if aux.same_rev_check(atom_pair, bond_symbols):
                    atom_pair_bond = bond
                    break  # Break the loop as atom pair takes priority
                
                if aux.same_rev_check(symbol_pair, bond_symbols):
                    symb_pair_bond = bond
            
            if atom_pair_bond:
                bond_list.append([i, j])
                bond_type.append(atom_pair_bond.id)
            elif symb_pair_bond:
                bond_list.append([i, j])
                bond_type.append(symb_pair_bond.id)
        
        return [np.array(bond_list).tolist(), bond_type]
    
    @property
    def angles(self):
        # Find all unique paths of length 2 (which corresponds to angles)
        angles = aux.find_unique_paths_of_length(self.graph, 2)
        tags = self.atoms.get_tags()
        
        angle_list = []
        angle_type = []
        
        for i, j, k in angles:
            a1 = self._atom_types[tags[i]]
            a2 = self._atom_types[tags[j]]
            a3 = self._atom_types[tags[k]]
            
            s1 = a1.split("_")[0]
            s2 = a2.split("_")[0]
            s3 = a3.split("_")[0]
            
            atom_triplet_bond = None
            symb_triplet_bond = None
            
            for angle in self._atype:
                angle_symbols = angle.symbols
                atom_triplet = (a1, a2, a3)
                symbol_triplet = (s1, s2, s3)
                
                if aux.same_rev_check(atom_triplet, angle_symbols):
                    atom_triplet_bond = angle
                    break  # Break the loop as atom triplet takes priority
                
                if aux.same_rev_check(symbol_triplet, angle_symbols):
                    symb_triplet_bond = angle
            
            if atom_triplet_bond:
                angle_list.append([i, j, k])
                angle_type.append(atom_triplet_bond.id)
            elif symb_triplet_bond:
                angle_list.append([i, j, k])
                angle_type.append(symb_triplet_bond.id)
        
        return [np.array(angle_list).tolist(), angle_type]
    
    @property
    def dihedrals(self):

        dihedrals = aux.find_unique_paths_of_length(self.graph, 3)
        tags = self.atoms.get_tags()
        
        dihedral_list = []
        dihedral_type = []
        for i, j, k, l in dihedrals:
            a1 = self._atom_types[tags[i]]
            a2 = self._atom_types[tags[j]]
            a3 = self._atom_types[tags[k]]
            a4 = self._atom_types[tags[l]]
            
            s1 = a1.split("_")[0]
            s2 = a2.split("_")[0]
            s3 = a3.split("_")[0]
            s4 = a4.split("_")[0]
            
            atom_quadruplet_bond = None
            symb_quadruplet_bond = None
            
            for dihedral in self._dtype:
                dihedral_symbols = dihedral.symbols
                atom_quadruplet = (a1, a2, a3, a4)
                symbol_quadruplet = (s1, s2, s3, s4)
                
                if aux.same_rev_check(atom_quadruplet, dihedral_symbols):
                    atom_quadruplet_bond = dihedral
                    break  # Break the loop as atom quadruplet takes priority
                
                if aux.same_rev_check(symbol_quadruplet, dihedral_symbols):
                    symb_quadruplet_bond = dihedral
            
            if atom_quadruplet_bond:
                dihedral_list.append([i, j, k, l])
                dihedral_type.append(atom_quadruplet_bond.id)
            elif symb_quadruplet_bond:
                dihedral_list.append([i, j, k, l])
                dihedral_type.append(symb_quadruplet_bond.id)
        
        return [np.array(dihedral_list).tolist(), dihedral_type]
    
    @property
    def impropers(self):
        
        imp_list = []
        imp_type = []
        
        # all_bonds = self.ana.all_bonds[0]
        all_bonds = [list(self.graph.neighbors(ii))for ii in range(len(self.atoms))]
        
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
        for tag in self.atoms.get_tags():
            atom_types.append(self._stype[tag].extended_label)
            
        if return_index:
            return atom_types, self.atoms.get_tags()
        
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

    def suggest_missing_interactions(self, type="all"):
        # Get current bonds, angles, and dihedrals
        current_bonds = set(tuple(bond) for bond in self.bonds[0])
        current_angles = set(tuple(angle) for angle in self.angles[0])
        current_dihedrals = set(tuple(dihedral) for dihedral in self.dihedrals[0])
        
        # Initialize sets to store suggested interactions
        suggested_bonds = set()
        suggested_angles = set()
        suggested_dihedrals = set()
        suggested_lj = {}
        
        # Get all possible bonds, angles, and dihedrals based on the graph
        possible_bonds = aux.find_unique_paths_of_length(self.graph, 1)
        possible_angles = aux.find_unique_paths_of_length(self.graph, 2)
        possible_dihedrals = aux.find_unique_paths_of_length(self.graph, 3)
        
        tags = self.atoms.get_tags()
        
        # Check for missing bonds
        for bond in possible_bonds:
            bond_tuple = tuple(bond)
            bond_tuple_rev = tuple(reversed(bond))
            if bond_tuple not in current_bonds and bond_tuple_rev not in current_bonds:
                a1, a2 = bond
                bond_type = f"{self._atom_types[tags[a1]]}-{self._atom_types[tags[a2]]}"
                suggested_bonds.add(bond_type)
        
        # Check for missing angles
        for angle in possible_angles:
            angle_tuple = tuple(angle)
            angle_tuple_rev = tuple(reversed(angle))
            if angle_tuple not in current_angles and angle_tuple_rev not in current_angles:
                a1, a2, a3 = angle
                angle_type = f"{self._atom_types[tags[a1]]}-{self._atom_types[tags[a2]]}-{self._atom_types[tags[a3]]}"
                suggested_angles.add(angle_type)
        
        # Check for missing dihedrals
        for dihedral in possible_dihedrals:
            dihedral_tuple = tuple(dihedral)
            dihedral_tuple_rev = tuple(reversed(dihedral))
            if dihedral_tuple not in current_dihedrals and dihedral_tuple_rev not in current_dihedrals:
                a1, a2, a3, a4 = dihedral
                dihedral_type = f"{self._atom_types[tags[a1]]}-{self._atom_types[tags[a2]]}-{self._atom_types[tags[a3]]}-{self._atom_types[tags[a4]]}"
                suggested_dihedrals.add(dihedral_type)
        
        # Check for missing LJ parameters
        for atom in self._stype:
            if atom.eps is None or atom.sig is None:
                suggested_lj[atom.label] = (None, None)  # Placeholder for missing LJ parameters
        
        suggestion = {
            "bonds": list(suggested_bonds),
            "angles": list(suggested_angles),
            "dihedrals": list(suggested_dihedrals),
            "lj": suggested_lj
            }
        
        if type == "all":
            return suggestion
        
        return suggestion[type]
