#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 13:58:43 2024

@author: roncofaber
"""

import ase
import ase.build
import ase.visualize
from ase.data import vdw_radii
from ase.data.colors import jmol_colors
import MDAnalysis as mda

import mdinterface.utils.map as pmap
from mdinterface.core.topology import Atom
from mdinterface.io.read import read_lammps_data_file
from mdinterface.utils.auxiliary import as_list, find_smallest_missing, round_list_to_sum
from mdinterface.externals import run_ligpargen, run_OBChargeModel, calculate_RESP_charges,\
    relax_structure, run_aimd
    
from mdinterface.utils.graphs import molecule_to_graph, find_unique_paths_of_length,\
    find_improper_idxs,find_relevant_distances, find_atom_types
from mdinterface.utils.draw import draw_bond_markers

import copy
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt

#%%

class Specie(object):
    
    def __init__(self, atoms=None, charges=None, atom_types=None, bonds=None,
                 angles=None, dihedrals=None, impropers=None, lj={}, cutoff=1.0,
                 name=None, lammps_data=None, fix_missing=False, chg_scaling=1.0,
                 pbc=False, ligpargen=False, tot_charge=None, prune_z=False):
        
        # store int. variables
        self.cutoff = cutoff
        
        # if file provided, read it
        if lammps_data is not None:
            atoms, atom_types, bonds, angles, dihedrals, impropers = read_lammps_data_file(lammps_data)
        else:
            # read/update atoms atoms
            atoms, atom_types_tmp = self._read_atoms(atoms, charges,
                                                 chg_scaling=chg_scaling, pbc=pbc)
            if atom_types is None:
                atom_types = atom_types_tmp
        
        # setup nominal charge
        if tot_charge is None:
            if not "nominal_charge" in atoms.arrays:
                atoms.set_array("nominal_charge", np.array(len(atoms)*[0]))
            tot_charge = int(np.sum(atoms.arrays["nominal_charge"]))
                
        # run ligpargen to calculate parameters
        if ligpargen:
            atoms, atom_types, bonds, angles, dihedrals, impropers = run_ligpargen(
                atoms, charge=tot_charge)
        
        # assign name
        if name is None:
            name = atoms.get_chemical_formula()
        if len(name) > 4:
            print("ATTENTION: resname for Specie could be misleading")
            print("You wan the name to be 4 characters long otherwise packmol")
            print("will have problems...")
        self.resname = name[:4]
        
        # set up atoms and generate graph of specie
        self.set_atoms(atoms, cutoff, prune_z=prune_z)
        
        # read atom_types from LJ
        if atom_types is None:
            atom_types = self._atom_types_from_lj(lj)
            
        # set up internal topology attributes
        self._setup_topology(atom_types, bonds, angles, dihedrals, impropers)
        
        if fix_missing:
            self._fix_missing_interactions()
            
        # initialize topology indexing
        self._update_topology_indexing()
        
        # set tot charge
        self._tot_charge = tot_charge

        return
    
    # read atoms to return ase.Atoms
    @staticmethod
    def _read_atoms(atoms, charges, chg_scaling=1.0, pbc=False):
        
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
    
    def _setup_topology(self, atoms, bonds, angles, dihedrals, impropers):
        """
        Setup topology from scratch - this is run in the beginning
        """
        
        # define atoms
        atoms_list, atom_map, atom_ids = pmap.map_atoms(as_list(atoms))
        
        self._stype = atoms_list
        self._smap = atom_map
        self._sids = atom_ids
        
        # initialize topology values
        self._old_bonds = copy.deepcopy(as_list(bonds))
        self._old_angles = copy.deepcopy(as_list(angles))
        self._old_dihedrals = copy.deepcopy(as_list(dihedrals))
        self._old_impropers = copy.deepcopy(as_list(impropers))
        
        # map list of inputs
        self._update_topology_mappings()
        
        return
    
    def _update_topology_mappings(self):
        """
        Update topology mappings.
        """
        # Remap topology lists
        bonds_list, bond_map, bond_ids = pmap.map_bonds(self._old_bonds)
        angles_list, angle_map, angle_ids = pmap.map_angles(self._old_angles)
        dihedrals_list, dihedral_map, dihedral_ids = pmap.map_dihedrals(self._old_dihedrals)
        impropers_list, improper_map, improper_ids = pmap.map_impropers(self._old_impropers)

        self._btype = bonds_list
        self._atype = angles_list
        self._dtype = dihedrals_list
        self._itype = impropers_list

        self._bmap = bond_map
        self._amap = angle_map
        self._dmap = dihedral_map
        self._imap = improper_map

        self._bids = bond_ids
        self._aids = angle_ids
        self._dids = dihedral_ids
        self._iids = improper_ids

        # Update topology to ensure consistency
        self._update_topology_indexing()

        return
    
    def _update_topology_indexing(self):
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

    def _add_to_topology(self, bonds=[], angles=[], dihedrals=[], impropers=[]):
        
        # check for uniqueness
        new_bonds = []
        new_angles = []
        new_dihedrals = []
        new_impropers = []
        for bond in copy.deepcopy(bonds):
            if not any(bond.__eq_strict__(bb) for bb in self._old_bonds):
                new_bonds.append(bond)
        
        for angle in copy.deepcopy(angles):
            if not any(angle.__eq_strict__(aa) for aa in self._old_angles):
                new_angles.append(angle)
                
        for dihedral in copy.deepcopy(dihedrals):
            if not any(dihedral.__eq_strict__(dd) for dd in self._old_dihedrals):
                new_dihedrals.append(dihedral)
                
        for improper in copy.deepcopy(impropers):
            if not any(improper.__eq_strict__(ii) for ii in self._old_impropers):
                new_impropers.append(improper)
        
        # add to topology
        self._old_bonds     += new_bonds
        self._old_angles    += new_angles
        self._old_dihedrals += new_dihedrals
        self._old_impropers += new_impropers
        
        # update topology mapping
        self._update_topology_mappings()
        
        return
    
    def _cleanup_topology(self):
        """
        Remove topology interactions that reference deleted atoms.

        This method should be called after atoms are deleted from the system
        to ensure topology consistency.

        Parameters:
        deleted_atom_ids (list): List of atom IDs (from self._sids) that were deleted
        """
        
        atom_ids = self._sids

        # Clean bonds
        valid_bonds = []
        for bond in self._old_bonds:
            if all(sim in atom_ids for sim in bond.symbols):
                valid_bonds.append(bond)
        self._old_bonds = valid_bonds

        # Clean angles
        valid_angles = []
        for angle in self._old_angles:
            if all(sim in atom_ids for sim in angle.symbols):
                valid_angles.append(angle)
        self._old_angles = valid_angles

        # Clean dihedrals
        valid_dihedrals = []
        for dihedral in self._old_dihedrals:
            if all(sim in atom_ids for sim in dihedral.symbols):
                valid_dihedrals.append(dihedral)
        self._old_dihedrals = valid_dihedrals

        # Clean impropers
        valid_impropers = []
        for improper in self._old_impropers:
            if all(sim in atom_ids for sim in improper.symbols):
                valid_impropers.append(improper)
        self._old_impropers = valid_impropers

        # Update topology mappings
        self._update_topology_mappings()

        return

    def _fix_missing_interactions(self):
        
        mss_bnd = pmap.generate_missing_interactions(self, "bonds")
        mss_ang = pmap.generate_missing_interactions(self, "angles")
        mss_dih = pmap.generate_missing_interactions(self, "dihedrals")
        mss_imp = pmap.generate_missing_interactions(self, "impropers")
        
        self._add_to_topology(mss_bnd, mss_ang, mss_dih, mss_imp)
        
        return
    
    # function to setup atom types
    def _atom_types_from_lj(self, lj):
        
        # use function to retrieve IDs
        atom_type_ids, types_map = find_atom_types(self.graph, max_depth=1)
        
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
    
    def round_charges(self, nround=7):
        
        rcharges = round_list_to_sum(self.charges, round(sum(self.charges), nround), nround)
        self.atoms.set_initial_charges(rcharges)
        
        return
    
    def set_atoms(self, atoms, cutoff=1.0, prune_z=False):
        
        # make sure no mess
        atoms = atoms.copy()
        
        if prune_z:
            zmin = atoms.get_positions(wrap=True)[:,2].min()
            zmax = atoms.get_positions(wrap=True)[:,2].max()
            atoms.translate([0,0,-zmin])
            atoms.cell[2][2] = zmax-zmin
            atoms.pbc[2] = False
        
        if not any(atoms.pbc):
            atoms.set_center_of_mass([0,0,0])
        
        self._atoms = atoms
        self._graph = molecule_to_graph(atoms, cutoff_scale=cutoff)

        return
    
    # covnert to mda.Universe
    def to_universe(self, charges=True, layered=False, match_cell=False, xydim=None):
        
        # empty top object
        top = mda.core.topology.Topology(n_atoms=len(self.atoms))
        
        # empty universe
        uni = mda.Universe(top, self.atoms.get_positions())
        
        # add some stuff
        uni.add_TopologyAttr("masses", self.atoms.get_masses())
        uni.add_TopologyAttr("resnames", [self.resname])
        # add atom symbols
        atom_symbols = self.atoms.get_chemical_symbols()  # Assuming this method exists
        uni.add_TopologyAttr("names", atom_symbols)
        
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
        
        # match the cell!
        if match_cell:
            assert xydim is not None
            zz = self.atoms.get_cell()[2][2]
            tatoms = self.atoms.copy()
            
            tatoms.set_cell([xydim[0], xydim[1], zz], scale_atoms=True)
            
            uni.dimensions = tatoms.cell.cellpar()
            uni.atoms.positions = tatoms.get_positions()
            
        # if has cell info, pass them along
        elif self.atoms.get_cell():
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
    
    def _find_interactions(self, path_length, tag_map, impropers=False):
        
        if not impropers:
            paths = find_unique_paths_of_length(self.graph, path_length)
        else:
            paths = find_improper_idxs(self.graph)
            
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
    def atoms(self):
        return self._atoms
    
    @property
    def graph(self):
        return self._graph
    
    @property
    def bonds(self):
        if not self._bmap: return [[], []]
        return self._find_interactions(1, self._bmap)
    
    @property
    def angles(self):
        if not self._amap: return [[], []]
        return self._find_interactions(2, self._amap)
    
    @property
    def dihedrals(self):
        if not self._dmap: return [[], []]
        return self._find_interactions(3, self._dmap)
    
    @property
    def impropers(self):
        if not self._imap: return [[], []]
        return self._find_interactions(3, self._imap, impropers=True)
        
    @property
    def _sids(self):
        return self.atoms.arrays["sids"]

    @_sids.setter
    def _sids(self, atom_ids):
        self.atoms.arrays["sids"] = atom_ids
            
    def copy(self):
        return copy.deepcopy(self)
    
    @property
    def charges(self):
        return self.atoms.get_initial_charges()
    
    def get_atom_types(self, return_index=False):
        
        type_indexes = np.array([self._smap[ii] for ii in self._sids])
        atom_types   = np.array([self._stype[ii].extended_label for ii in type_indexes])
            
        if return_index:
            return atom_types, type_indexes
        
        return atom_types
    
    # method to estimate the volume of the specie
    def estimate_specie_volume(self, probe_radius=0):
        
        try:
            from libarvo import molecular_vs
        except:
            raise ImportError("libarvo NOT found. Install it.")

        centers = self.atoms.get_positions()
        radii = [vdw_radii[ii] for ii in self.atoms.get_atomic_numbers()]

        volume, surface = molecular_vs(centers, radii, probe_radius)

        return volume
    
    # method to estimate the radius of the specie if it were a sphere
    def estimate_specie_radius(self, probe_radius=0):
        volume = self.estimate_specie_volume(probe_radius=probe_radius)
        return (3 * volume / (4 * np.pi)) ** (1/3)
    
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
        
        unique_pairs_list = find_relevant_distances(self.graph, Nmax, Nmin=Nmin,
                                                    centers=centers, Ninv=Ninv)
        return unique_pairs_list

    def plot_graph(self, show_bonds=False):
        
        colors = [jmol_colors[a.number] for a in self.atoms]
        
        # node_pos = nx.spring_layout(self.graph, k=1.5/np.sqrt(self.graph.order()))
        node_pos = nx.kamada_kawai_layout(self.graph)
        
        fig, ax = plt.subplots()
        nx.draw(self.graph, pos=node_pos, with_labels=True, node_color=colors,
                node_size=1000, edge_color='black', linewidths=2, font_size=15,
                edgecolors="black", ax=ax, width=2)
        
        if show_bonds:
            draw_bond_markers(ax, self, node_pos, jmol_colors)
    
        plt.show()
        
        return

    def update_positions(self, positions=None, cellpar=None, atoms=None, prune_z=False):
        
        if atoms is not None:
            positions = atoms.get_positions()
            cellpar   = atoms.get_cell()
        
        atoms = self._atoms.copy()
        
        if positions is not None:
            atoms.set_positions(positions)
        if cellpar is not None:
            atoms.set_cell(cellpar)
            
        self.set_atoms(atoms, cutoff=self.cutoff, prune_z=prune_z)
        
        return
    
    def estimate_charges(self, method="obabel", charge=None, assign=False, **respargs):
        
        if charge is None:
            charge = 0
            print("YOU WANT TO CHECK CHARGE BUD.")
        
        if method == "obabel":
            charges = run_OBChargeModel(self.atoms)
        elif method == "ligpargen":
            system, _, _, _, _, _ = run_ligpargen(self.atoms, charge=charge)
            charges = system.get_initial_charges()
        elif method == "resp":
            charges, atoms = calculate_RESP_charges(self, **respargs)
            
        if assign:
            self.atoms.set_positions(atoms.get_positions())
            self.atoms.set_initial_charges(charges)
        
        return charges

    def estimate_OPLSAA_parameters(self, charge=None):

        if charge is None:
            charge = 0

        system, atoms, bonds, angles, dihedrals, impropers = run_ligpargen(
            self.atoms, charge=charge)

        return system, atoms, bonds, angles, dihedrals, impropers

    def relax_structure(self, optimizer='FIRE', fmax=0.05,
                       steps=200, update_positions=True, trajectory=None,
                       logfile=None, **kwargs):
        """
        Relax the Specie structure using ASE or UMA optimization.

        Parameters
        ----------
        method : str, default 'ase'
            Relaxation method to use ('ase', 'uma')
        optimizer : str, default 'BFGS'
            Optimizer for ASE method ('BFGS', 'LBFGS', 'FIRE')
        fmax : float, default 0.05
            Maximum force threshold for convergence (eV/Å)
        steps : int, default 200
            Maximum number of optimization steps
        update_positions : bool, default True
            Whether to update the Specie's atomic positions after relaxation
        trajectory : str, optional
            Path to save optimization trajectory
        logfile : str, optional
            Path to save optimization log
        **kwargs
            Additional arguments passed to the relaxation function

        Returns
        -------
        relaxed_atoms : ase.Atoms
            The relaxed atomic structure
        converged : bool
            Whether the optimization converged
        """
        
        atoms_to_relax = self.atoms
        atoms_to_relax.info["charge"] = self._tot_charge
        atoms_to_relax.info["spin"]   = 0 #TODO
        
        relaxed_atoms = relax_structure(
            atoms_to_relax, optimizer=optimizer, fmax=fmax, steps=steps,
            trajectory=trajectory, logfile=logfile, **kwargs)

        # Update positions if requested
        if update_positions:
            self.update_positions(atoms=relaxed_atoms, prune_z=False)

        return
    
    def run_aimd(self, timestep=0.5, temperature_K=300, friction=0.001,
              steps=1000, update_positions=True, trajectory=None,
              logfile=None, **kwargs):
        """
        Run AIMD for the Specie structure using ASE and FAIRChem.
    
        Parameters
        ----------
        timestep : float, default 0.1
            Time step for the simulation in fs.
        temperature_K : float, default 300
            Target temperature for the Langevin dynamics in Kelvin.
        friction : float, default 0.001
            Frictional damping coefficient in 1/fs.
        steps : int, default 1000
            Number of time steps to run the AIMD.
        update_positions : bool, default True
            Whether to update the Specie's atomic positions after AIMD.
        trajectory : str, optional
            Path to save the MD trajectory.
        logfile : str, optional
            Path to save the MD log.
        **kwargs
            Additional arguments passed to the Langevin integrator.
    
        Returns
        -------
        None
        """
        
        atoms_to_simulate = self.atoms
        atoms_to_simulate.info["charge"] = self._tot_charge
        atoms_to_simulate.info["spin"] = 0  # Update this if needed
    
        # Call the previously defined function to perform AIMD
        run_aimd(atoms_to_simulate, timestep=timestep, temperature_K=temperature_K,
                 friction=friction, steps=steps, trajectory=trajectory,
                 logfile=logfile, **kwargs)
    
        # Update positions if requested
        if update_positions:
            self.update_positions(atoms=atoms_to_simulate, prune_z=False)
        
        return
    
    def _find_rings(self, max_ring_size=8):
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
            for cycle in nx.simple_cycles(self.graph):
                if len(cycle) <= max_ring_size:
                    rings.append(sorted(cycle))
        except:
            # Fallback: use minimum cycle basis for undirected graphs
            try:
                cycle_basis = nx.minimum_cycle_basis(self.graph)
                for cycle in cycle_basis:
                    if len(cycle) <= max_ring_size:
                        rings.append(sorted(cycle))
            except:
                # If all else fails, return empty list
                rings = []

        return rings

    def _get_rings_containing_atoms(self, atom_indices, max_ring_size=8):
        """
        Find all rings that contain any of the specified atoms.

        Parameters:
        atom_indices (list): List of atom indices to check
        max_ring_size (int): Maximum ring size to detect

        Returns:
        list: List of rings containing any of the specified atoms
        """
        all_rings = self._find_rings(max_ring_size)
        relevant_rings = []

        atom_set = set(atom_indices)
        for ring in all_rings:
            if atom_set.intersection(set(ring)):
                relevant_rings.append(ring)

        return relevant_rings

