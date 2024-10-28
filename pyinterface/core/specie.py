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
import pyinterface.utils.map as pmap

import copy
import numpy as np
#%%

class Specie(object):
    
    def __init__(self, atoms=None, charges=None, atom_types=None, bonds=None,
                 angles=None, dihedrals=None, impropers=None, lj={}, cutoff=1.0,
                 name=None, lammps_data=None):
        
        # if file provided, read it
        if lammps_data is not None:
            atoms, atom_types, bonds, angles, dihedrals = read_lammps_data_file(lammps_data)

        # read/update atoms atoms
        atoms = self._read_atoms(atoms, charges)
        
        # assign name
        if name is None:
            name = atoms.get_chemical_formula()
        if len(name) > 4:
            print("ATTENTION: resname for Specie could be misleading")
        self.resname = name
        
        # set up atoms and generate graph of specie
        self.set_atoms(atoms, cutoff)
        
        # read atom_types from LJ
        if atom_types is None:
            atom_types = self._atom_types_from_lj(lj)
            
        # set up internal topology attributes
        self._setup_topology(atom_types, bonds, angles, dihedrals, impropers)
        
        # initialize topology info
        self._update_topology()

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
    
    def _setup_topology(self, atoms, bonds, angles, dihedrals, impropers):
        
        # map list of inputs
        atoms_list, atom_map, atom_ids = pmap.map_atoms(as_list(atoms))
        bonds_list, bond_map, bond_ids = pmap.map_bonds(as_list(bonds))
        angles_list, angle_map, angle_ids = pmap.map_angles(as_list(angles))
        dihedrals_list, dihedral_map, dihedral_ids = pmap.map_dihedrals(as_list(dihedrals))
        impropers_list, improper_map, improper_ids = pmap.map_impropers(as_list(impropers))
        
        self._btype = bonds_list
        self._atype = angles_list
        self._dtype = dihedrals_list
        self._itype = impropers_list
        self._stype = atoms_list
        
        self._smap = atom_map
        self._bmap = bond_map
        self._amap = angle_map
        self._dmap = dihedral_map
        self._imap = improper_map
        
        self._sids = atom_ids
        self._bids = bond_ids
        self._aids = angle_ids
        self._dids = dihedral_ids
        self._iids = improper_ids
        
        # self._atom_types = types_map
        # self.atoms.set_tags(atom_type_ids)
        return
    
    # function to setup atom types
    def _atom_types_from_lj(self, lj):
        
        # use function to retrieve IDs
        atom_type_ids, types_map = aux.find_atom_types(self.atoms, max_depth=1)
        
        atom_types = []
        for atom_id in atom_type_ids:
            
            atom_symbol = types_map[atom_id][0]
            atom_neighs = "".join(types_map[atom_id][1])
            
            label = "{}_{}".format(atom_symbol, atom_neighs)
            # types_map[atom_type] = label
            
            if label in lj:
                eps, sig = lj[label]
            elif atom_symbol in lj:
                eps, sig = lj[atom_symbol]
            else:
                eps, sig = None,None
            
            atom = Atom(atom_symbol, label=label, eps=eps, sig=sig)
            atom_types.append(atom)
                
        return atom_types
    
    def set_atoms(self, atoms, cutoff=1.0):
        
        self._atoms = atoms
        
        self._graph = aux.molecule_to_graph(atoms)
        
        return
    
    def _update_topology(self):
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
        self._update_topology()
        
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
            att_types = self.type2id(att, types)
            uni._add_topology_objects(att, attribute, types=att_types)
        
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
        bond_tags = self._bmap
        bond_list = []
        bond_type = []
        
        for i, j in bonds:
            a1 = self._sids[i]
            a2 = self._sids[j]
            
            s1 = a1.split("_")[0]
            s2 = a2.split("_")[0]
            
            if (a1, a2) in bond_tags:
                bond_list.append([i, j])
                bond_type.append(bond_tags[(a1, a2)])
            elif (a2, a1) in bond_tags:
                bond_list.append([j, i])
                bond_type.append(bond_tags[(a2, a1)])
            elif (s1, s2) in bond_tags:
                bond_list.append([i, j])
                bond_type.append(bond_tags[(s1, s2)])
            elif (s2, s1) in bond_tags:
                bond_list.append([j, i])
                bond_type.append(bond_tags[(s2, s1)])
        
        bond_list = np.array(bond_list, dtype=int)
        bond_type = np.array(bond_type, dtype=int)
        
        return [bond_list.tolist(), bond_type.tolist()]
    
    @property
    def angles(self):
        # Find all unique paths of length 2 (which corresponds to angles)
        angles = aux.find_unique_paths_of_length(self.graph, 2)
        angle_tags = self._amap
        angle_list = []
        angle_type = []
        
        for i, j, k in angles:
            a1 = self._sids[i]
            a2 = self._sids[j]
            a3 = self._sids[k]
            
            s1 = a1.split("_")[0]
            s2 = a2.split("_")[0]
            s3 = a3.split("_")[0]
            
            if (a1, a2, a3) in angle_tags:
                angle_list.append([i, j, k])
                angle_type.append(angle_tags[(a1, a2, a3)])
            elif (a3, a2, a1) in angle_tags:
                angle_list.append([k, j, i])
                angle_type.append(angle_tags[(a3, a2, a1)])
            elif (s1, s2, s3) in angle_tags:
                angle_list.append([i, j, k])
                angle_type.append(angle_tags[(s1, s2, s3)])
            elif (s3, s2, s1) in angle_tags:
                angle_list.append([k, j, i])
                angle_type.append(angle_tags[(s3, s2, s1)])

            
        angle_list = np.array(angle_list, dtype=int)
        angle_type = np.array(angle_type, dtype=int)
        
        return [angle_list.tolist(), angle_type.tolist()]
    
    @property
    def dihedrals(self):
        # Find all unique paths of length 3 (which corresponds to dihedrals)
        dihedrals = aux.find_unique_paths_of_length(self.graph, 3)
        dihedral_tags = self._dmap
        dihedral_list = []
        dihedral_type = []
        
        for i, j, k, l in dihedrals:
            a1 = self._sids[i]
            a2 = self._sids[j]
            a3 = self._sids[k]
            a4 = self._sids[l]
            
            s1 = a1.split("_")[0]
            s2 = a2.split("_")[0]
            s3 = a3.split("_")[0]
            s4 = a4.split("_")[0]
            
            if (a1, a2, a3, a4) in dihedral_tags:
                dihedral_list.append([i, j, k, l])
                dihedral_type.append(dihedral_tags[(a1, a2, a3, a4)])
            elif (a4, a3, a2, a1) in dihedral_tags:
                dihedral_list.append([l, k, j, i])
                dihedral_type.append(dihedral_tags[(a4, a3, a2, a1)])
            elif (s1, s2, s3, s4) in dihedral_tags:
                dihedral_list.append([i, j, k, l])
                dihedral_type.append(dihedral_tags[(s1, s2, s3, s4)])
            elif (s4, s3, s2, s1) in dihedral_tags:
                dihedral_list.append([l, k, j, i])
                dihedral_type.append(dihedral_tags[(s4, s3, s2, s1)])
        
        dihedral_list = np.array(dihedral_list, dtype=int)
        dihedral_type = np.array(dihedral_type, dtype=int)
        
        return [dihedral_list.tolist(), dihedral_type.tolist()]
    
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
        
        type_indexes = np.array([self._smap[ii] for ii in self._sids])
        atom_types   = np.array([self._stype[ii].extended_label for ii in type_indexes])
        # atom_types = np.array([ii + "_" + self.resname for ii in self._sids])
        # atom_types = np.array([atom.extended_label for atom in self._stype])
            
        if return_index:
            return atom_types, type_indexes
        
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

    def type2id(self, attribute, types):
        
        # Get the corresponding attribute list
        if attribute == "bonds":
            attribute_list = self._btype
        elif attribute == "angles":
            attribute_list = self._atype
        elif attribute == "dihedrals":
            attribute_list = self._dtype
        elif attribute == "impropers":
            attribute_list = self._itype
        elif attribute == "atoms":
            attribute_list = self._stype
        else:
            raise ValueError("Invalid topology attribute type.")
        
        return [attribute_list[idx].id for idx in types]
    
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
