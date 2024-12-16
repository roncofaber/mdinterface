#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 13:58:43 2024

@author: roncofaber
"""

import ase
import ase.build
import ase.visualize
from ase.data.colors import jmol_colors
import MDAnalysis as mda

from mdinterface.utils.auxiliary import as_list, find_smallest_missing, atoms_to_indexes
    
from mdinterface.core.topology import Atom
import mdinterface.utils.auxiliary as aux
from mdinterface.io.read import read_lammps_data_file
import mdinterface.utils.map as pmap

import copy
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt

#%%

class Specie(object):
    
    def __init__(self, atoms=None, charges=None, atom_types=None, bonds=None,
                 angles=None, dihedrals=None, impropers=None, lj={}, cutoff=1.0,
                 name=None, lammps_data=None, fix_missing=False, chg_scaling=1.0):
        
        # if file provided, read it
        if lammps_data is not None:
            atoms, atom_types, bonds, angles, dihedrals, impropers = read_lammps_data_file(lammps_data)
        else:
            # read/update atoms atoms
            atoms, atom_types = self._read_atoms(atoms, charges, chg_scaling=chg_scaling)
        
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
        self._setup_topology(atom_types, bonds, angles, dihedrals, impropers,
                             fix_missing=fix_missing)
        
        # initialize topology info
        self._update_topology()

        return
    
    # read atoms to return ase.Atoms
    @staticmethod
    def _read_atoms(atoms, charges, chg_scaling=1.0):
        
        # initialize atoms obj
        if isinstance(atoms, str):
            try:
                atoms = ase.io.read(atoms)
            except:
                atoms = ase.build.molecule(atoms)
        elif isinstance(atoms, ase.Atoms):
            atoms = atoms.copy()
        
        # assign charges
        if charges is not None:
            charges = as_list(charges)
            if len(charges) == 1:
                charges = len(atoms)*charges
        else:
            charges = atoms.get_initial_charges()
            
        # rescale charges if needed
        charges = chg_scaling*np.asarray(charges)
        atoms.set_initial_charges(charges)
        
        # see if stype is already present
        if "stype" in atoms.arrays:
            stype = atoms.arrays["stype"]
        else:
            stype = None

        return atoms, stype
    
    def _setup_topology(self, atoms, bonds, angles, dihedrals, impropers,
                        fix_missing=False):
        
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
        
        if fix_missing: #ugly but seems to work
            
            mss_bnd = pmap.generate_missing_interactions(self, "bonds")
            mss_ang = pmap.generate_missing_interactions(self, "angles")
            mss_dih = pmap.generate_missing_interactions(self, "dihedrals")
            mss_imp = pmap.generate_missing_interactions(self, "impropers")
            
            # map list of inputs
            # atoms_list, atom_map, atom_ids = pmap.map_atoms(as_list(atoms))
            bonds_list, bond_map, bond_ids = pmap.map_bonds(as_list(bonds) + mss_bnd)
            angles_list, angle_map, angle_ids = pmap.map_angles(as_list(angles) + mss_ang)
            dihedrals_list, dihedral_map, dihedral_ids = pmap.map_dihedrals(as_list(dihedrals) + mss_dih)
            impropers_list, improper_map, improper_ids = pmap.map_impropers(as_list(impropers) + mss_imp)
            
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
        
        self._graph = aux.molecule_to_graph(atoms, cutoff_scale=cutoff)
        
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
            att_types = self._type2id(att, types)
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
    
    def _find_interactions(self, path_length, tag_map, impropers=False):
        
        if not impropers:
            paths = aux.find_unique_paths_of_length(self.graph, path_length)
        else:
            paths = aux.find_improper_idxs(self.graph)
            
        interaction_list = []
        interaction_type = []
        
        for indices in paths:
            atoms = [self._sids[idx] for idx in indices]
            symbols = [atom.split("_")[0] for atom in atoms]
            
            if tuple(atoms) in tag_map:
                interaction_list.append(indices)
                interaction_type.append(tag_map[tuple(atoms)])
            elif tuple(atoms[::-1]) in tag_map:
                interaction_list.append(indices[::-1])
                interaction_type.append(tag_map[tuple(atoms[::-1])])
            elif tuple(symbols) in tag_map:
                interaction_list.append(indices)
                interaction_type.append(tag_map[tuple(symbols)])
            elif tuple(symbols[::-1]) in tag_map:
                interaction_list.append(indices[::-1])
                interaction_type.append(tag_map[tuple(symbols[::-1])])
        
        interaction_list = np.array(interaction_list, dtype=int)
        interaction_type = np.array(interaction_type, dtype=int)
        
        return [interaction_list.tolist(), interaction_type.tolist()]
    
    @property
    def bonds(self):
        return self._find_interactions(1, self._bmap)
    
    @property
    def angles(self):
        return self._find_interactions(2, self._amap)
    
    @property
    def dihedrals(self):
        return self._find_interactions(3, self._dmap)
    
    @property
    def impropers(self):
        return self._find_interactions(3, self._imap, impropers=True)
        
    @property
    def _sids(self):
        return self.atoms.arrays["sids"]

    @_sids.setter
    def _sids(self, atom_ids):
        self.atoms.arrays["sids"] = atom_ids
            
    def copy(self):
        return copy.deepcopy(self)
    
    def get_atom_types(self, return_index=False):
        
        type_indexes = np.array([self._smap[ii] for ii in self._sids])
        atom_types   = np.array([self._stype[ii].extended_label for ii in type_indexes])
            
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

    def _type2id(self, attribute, types):
        
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
    
    def suggest_missing_interactions(self, stype="all"):
        
    
        # Check for missing interactions
        missing_bonds = pmap.find_missing_bonds(self)
        missing_angles = pmap.find_missing_angles(self)
        missing_dihedrals = pmap.find_missing_dihedrals(self)
        missing_impropers = []#pmap.find_missing_impropers(self) #TODO change
        
        suggestions = {
            "bonds": missing_bonds,
            "angles": missing_angles,
            "dihedrals": missing_dihedrals,
            "impropers": missing_impropers
        }
        
        if stype == "all":
            return suggestions
        
        return suggestions[stype]

    def __repr__(self):
        return f"{self.__class__.__name__}({self.resname})"

    def find_relevant_distances(self, Nmax, Nmin=0, centers=None, Ninv=0):
        
        # get list of relevant nodes
        if centers is None:
            relevant_nodes = self.graph.nodes()
        else:
            relevant_nodes = aux.as_list(centers)
        
        # Set to store unique pairs
        unique_pairs = set()

        # Iterate over all nodes in the graph
        for node1 in relevant_nodes:
            # Get the shortest path lengths from node node to all other reachable nodes
            shortest_paths = nx.single_source_shortest_path_length(self.graph, node1)
            
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

    def plot_graph(self, **kwargs):
        
        colors = [jmol_colors[a.number] for a in self.atoms]
        
        fig, ax = plt.subplots()
        nx.draw(self.graph, with_labels=True, node_color=colors,
                node_size=1000, edge_color='black', linewidths=2, font_size=15,
                edgecolors="black", ax=ax, width=2)
        
        
        plt.show()
        
        return
