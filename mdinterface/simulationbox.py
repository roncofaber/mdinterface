#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 15:14:41 2023

@author: roncoroni
"""

from mdinterface.utils.auxiliary import label_to_element, as_list, find_smallest_missing
from mdinterface.io.lammpswriter import DATAWriter
from mdinterface.build.box import make_interface_slab, make_solvent_box, add_component

import ase

import numpy as np

import shutil
import warnings
warnings.filterwarnings('ignore')

#%%

default_params = {
    "interface": {
        "nlayers": 1,
    },
    "solvent": {
        "rho": None,
        "zdim": None,
        "nions": None,
        "concentration": None,
        "conmodel": None,
        "ion_pos": None
    },
    "miderface": {
        "nlayers": 1,
    },
    "enderface": {
        "nlayers": 1,
    },
    "vacuum": {
        "zdim": 0,
    }
}

class SimulationBox():
    
    def __init__(self, solvent=None, solute=None, interface=None,
                 enderface=None, miderface=None):
        
        # start species
        self._setup_species(solvent, solute, interface, enderface, miderface)
        
        # check interface indexing
        self._make_sandwich()
        
        # fix indexes of topology elements
        self._update_topology_indexes()
        
        return
    
    def _setup_species(self, solvent, solute, interface, enderface, miderface):
        
        self._solvent   = None
        self._solute    = None
        self._interface = None
        self._enderface = None
        self._miderface = None
        
        # assign variables
        if solvent is not None:
            self._solvent   = solvent.copy()
        if solute is not None:
            self._solute    = [ii.copy() for ii in as_list(solute)]
        if interface is not None:
            self._interface = interface.copy()
        if enderface is not None:
            self._enderface = enderface.copy() # use this to make a good sandwich!
        if miderface is not None:
            self._miderface = miderface.copy() # add some yummy stuffing!
        
        return
    
    def _make_sandwich(self):
        
        if self._interface is not None:
            for atom in self._interface._stype:
                atom.set_label(atom.label + "_i")
                
        if self._enderface is not None:
            for atom in self._enderface._stype:
                atom.set_label(atom.label + "_e")
                
        if self._miderface is not None:
            for atom in self._miderface._stype:
                atom.set_label(atom.label + "_m")   
                
        return
                    
    def _update_topology_indexes(self):
        
        nitems = {
            "_btype" : [],
            "_atype" : [],
            "_dtype" : [],
            "_itype" : [],
            }
        
        for attribute in nitems:
            for specie in self._species:
                for attr in specie.__getattribute__(attribute):
                    if attr.id not in nitems[attribute]:
                        nitems[attribute].append(attr.id)
                    else:
                        idx = find_smallest_missing(nitems[attribute], start=1)
                        attr.set_id(idx)
                        nitems[attribute].append(attr.id)
                        
        #resort atom types by alph order
        atom_types = []
        for specie in self._species:
            atom_types.extend([stype.extended_label for stype in specie._stype])
        atom_types.sort()
        
        for specie in self._species:
            for stype in specie._stype:
                idx = np.argwhere(stype.extended_label == np.array(atom_types))[0][0]
                stype.set_id(idx+1)
        
        return
    
    @staticmethod
    def _get_size_from_slab(slab):
        
        xsize = [1,0,0]@slab.atoms.cell@[1,0,0]
        ysize = [0,1,0]@slab.atoms.cell@[0,1,0]
        slab_depth = [0,0,1]@slab.atoms.cell@[0,0,1]
        
        return xsize, ysize, slab_depth
    
    # main driver to generate a simulation box given instructions
    def make_simulation_box(self, xysize, layering, padding=1.5, to_ase=False,
                            write_data=False, filename="data.lammps",
                            center_electrode=False, layered=False, hijack=None,
                            match_cell=False, remove_charges=False):
        
        # define approximate cross_section
        assert len(xysize) == 2, "'xysize' should have length of 2 [xsize, ysize]"
        xsize, ysize = xysize
        
        # make slabs of solid stuff to find actual cross_section
        slabs = []
        for layer in layering:
            layer_type = layer["type"]
            for k, v in default_params[layer_type].items():
                layer.setdefault(k, v)
            
            if layer_type in ["interface", "enderface", "miderface"]:
                tslab = getattr(self, "_" + layer["type"])
                tslab = make_interface_slab(tslab, xsize, ysize, layers=layer['nlayers'])
                slabs.append(tslab)
                layer["slab"] = tslab
        
        # define ACTUAL cell boundaries
        xsize, ysize = self._define_cell_boundaries(xsize, ysize, slabs)
        
        # build cake, layer by layer - starting from nothing
        system = None
        zdim   = 0
        for layer in layering:
            
            layer_type = layer["type"]
            
            if layer_type == "solvent":
                zsize    = layer["zdim"]
                solv_rho = layer["rho"]
                nions    = layer["nions"]
                concentration = layer["concentration"]
                conmodel = layer["conmodel"]
                ion_pos  = layer["ion_pos"]
                
                # make solvent box
                solvent = make_solvent_box(self.species, self.solvent, self._solute,
                                           [xsize, ysize, zsize], solv_rho,
                                           nions, concentration, conmodel, ion_pos)

                # add component
                system, zdim = add_component(system, solvent, zdim, padding=padding)
                
            elif layer_type in ["interface", "enderface", "miderface"]:
                tslab = layer["slab"].to_universe(layered=layered, match_cell=match_cell, xydim=[xsize, ysize])
                # add the slab
                system, zdim = add_component(system, tslab, zdim, padding=padding)
            
            elif layer_type == "vacuum":
                zdim += layer["zdim"]
            
        # update my system's dimensions        
        system.dimensions = [xsize, ysize, zdim] + [90, 90, 90] #TODO not like this
        
        # move of half unit along z
        if center_electrode:
            system.atoms.translate([0,0,zdim/2])
            _ = system.atoms.wrap()
            
        # give ase atoms to override positions
        if hijack is not None:
            system.dimensions = hijack.get_cell_lengths_and_angles()
            system.atoms.positions = hijack.get_positions()
        
        # write data file
        if write_data:
            self.write_lammps_file(system, filename=filename, remove_charges=remove_charges)
        
        # convert to ase, or not
        if to_ase:
            return self.to_ase(system)
        return system
    
    def write_lammps_file(self, system,  write_coeff=True, filename="data.lammps",
                          remove_charges=False):
        
        if remove_charges:
            system.del_TopologyAttr("charges")
        
        # first write data file
        with DATAWriter(filename) as dt:
            dt.write(system.atoms)
            
        # now write coeff where they belong
        if write_coeff:
            
            temp_file = 'tmp_data.lammps'
            
            with open(filename, 'r') as ffile, open(temp_file, 'w') as tfile:
                for ln, fl in enumerate(ffile):
                    if fl.startswith("Atoms"):
                        
                        # write coefficients
                        self.write_coefficients(system, fout=tfile)
                        
                    tfile.write(fl)
            
            shutil.move(temp_file, filename)
        
        return
    
    def write_coefficients(self, system, fname="tmp.coeff", fout=None):
        
        remember_to_close = False
        if fout is None:
            fout = open(fname, "w")
            remember_to_close = True
            
            
        fout.write("Pair Coeffs\n\n")
        
        idx = 1
        for cc, atom in enumerate(self.get_sorted_attribute("atoms")):

            if atom.extended_label not in np.unique(system.atoms.types):
                continue
            
            eps = atom.eps if atom.eps is not None else 0
            sig = atom.sig if atom.sig is not None else 0
            
            fout.write("{:>5}    {:>12.8f}    {:>12.8f}  # {}\n".format(
                idx, eps, sig, atom.extended_label))
            idx += 1
            
        if self.get_sorted_attribute("bonds"):
            fout.write("\n")
            fout.write("Bond Coeffs\n\n")
        
        for bond in self.get_sorted_attribute("bonds"):
            
            if bond.id not in np.array(system.bonds.types(), dtype=int):
                continue
            
            kr = bond.kr if bond.kr is not None else 0
            r0 = bond.r0 if bond.r0 is not None else 0
            
            btype = "{}-{}".format(*bond.symbols)
        
            fout.write("{:>5}    {:>10.6f}    {:>10.6f}  #  {:<5} | {}\n".format(
                bond.id, kr, r0, btype, bond.resname))
        
        if self.get_sorted_attribute("angles"):
            fout.write("\n")
            fout.write("Angle Coeffs\n\n")
        
        for angle in self.get_sorted_attribute("angles"):
            
            if angle.id not in np.array(system.angles.types(), dtype=int):
                continue
            
            kr     = angle.kr if angle.kr is not None else 0
            theta0 = angle.theta0 if angle.theta0 is not None else 0
            
            atype = "{}-{}-{}".format(*angle.symbols)
        
            fout.write("{:>5}    {:>10.6f}    {:>10.6f}  #  {:<8} | {}\n".format(
                angle.id, kr, theta0, atype, angle.resname))
        
        if self.get_sorted_attribute("dihedrals"):
            fout.write("\n")
            fout.write("Dihedral Coeffs\n\n")
        
        for dihedral in self.get_sorted_attribute("dihedrals"):
            
            if dihedral.id not in np.array(system.dihedrals.types(), dtype=int):
                continue
            
            dihedral.write(fout)
            # atype = "{}-{}-{}-{}".format(*dihedral.symbols)
            # value = "{:>7.4f}  {:>7.4f}  {:>7.4f}  {:>7.4f}  {:>7.4f}".format(*dihedral.values)
        
            # fout.write("{:>5}    {}  #  {:<8} | {}\n".format(dihedral.id, value, atype, dihedral.resname))
        
        if self.get_sorted_attribute("impropers"):
            fout.write("\n")
            fout.write("Improper Coeffs\n\n")
        
        for improper in self.get_sorted_attribute("impropers"):
            
            if improper.id not in np.array(system.impropers.types(), dtype=int):
                continue
            
            atype = "{}".format(*improper.symbols)
            value = "{:>7.4f}    {:>2d}    {:>2d}".format(*improper.values)
        
            fout.write("{:>5}    {}  #  {:<2} | {}\n".format(improper.id, value, atype, improper.resname))
        
        fout.write("\n")
        
        if remember_to_close:
            fout.close()
        
        return
    
    # convert to ase.Atoms
    @staticmethod
    def to_ase(system):
        
        if system is None:
            return ase.Atoms()
        
        positions = system.atoms.positions
        
        masses = system.atoms.masses
        labels = system.atoms.types
        
        symbols = [label_to_element(lab, mas) for lab, mas in zip(labels, masses)]
        
        ase_system = ase.Atoms(symbols=symbols, positions=positions)
        
        if system.dimensions is not None:
            ase_system.set_cell(system.dimensions)
            ase_system.set_pbc(True)
        
        try:
            if system.atoms.charges is not None:
                ase_system.set_initial_charges(system.atoms.charges)
        except:
            pass
        
        return ase_system

    @property
    def solvent(self):
        return self._solvent.to_universe() if self._solvent else None
    
    @property
    def solute(self):
        return [ii.to_universe() for ii in self._solute] if self._solute else None
    
    @property
    def interface(self):
        return self._interface.to_universe() if self._interface else None
    
    @property
    def enderface(self):
        return self._enderface.to_universe() if self._enderface else None
    
    @property
    def miderface(self):
        return self._miderface.to_universe() if self._miderface else None
    
    @property
    def species(self):
        return [ii.to_universe() for ii in self._species if ii is not None]
    
    @property
    def _species(self):
        all_species = np.concatenate((
            as_list(self._solvent), 
            as_list(self._solute),
            as_list(self._interface), 
            as_list(self._enderface),
            as_list(self._miderface)
        )).tolist()
        return all_species
    
    def get_sorted_attribute(self, attribute):
        
        if attribute.lower() in "bonds":
            attribute = "_btype"
        elif attribute.lower() in "angles":
            attribute = "_atype"
        elif attribute.lower() in "dihedrals":
            attribute = "_dtype"
        elif attribute.lower() in "impropers":
            attribute = "_itype"
        elif attribute.lower() in "atoms":
            attribute = "_stype"
        
        indexes    = []
        attributes = []
        for specie in self._species:
            for attr in specie.__getattribute__(attribute):
                indexes.append(attr.id)
                attributes.append(attr)
                
        return [attributes[ii] for ii in np.argsort(indexes)]
    
    def _define_cell_boundaries(self, xsize, ysize, slabs):
        
        if not slabs:
            return xsize, ysize
        
        xsize_t, ysize_t = np.NaN, np.NaN
        
        for tslab in slabs:
            xi, yi, _ = self._get_size_from_slab(tslab)
            xsize_t = np.nanmax([xi, xsize_t])
            ysize_t = np.nanmax([yi, ysize_t])
        
        return xsize_t, ysize_t
