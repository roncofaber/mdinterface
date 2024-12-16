#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 15:14:41 2023

@author: roncoroni
"""

from mdinterface.utils.auxiliary import label_to_element, as_list, find_smallest_missing
from mdinterface.io.lammpswriter import DATAWriter
from mdinterface.io.packmol import header, box_place, fix_place

import MDAnalysis as mda

import ase
from ase import units
import numpy as np

import subprocess
import shutil

import warnings
warnings.filterwarnings('ignore')

#%%

class SimulationBox():
    
    def __init__(self, solvent=None, solute=None, interface=None,
                 enderface=None):
        
        # start species
        self._setup_species(solvent, solute, interface, enderface)
        
        # check interface indexing
        self._make_sandwich()
        
        # fix indexes of topology elements
        self._update_topology_indexes()
        
        return
    
    
    def _setup_species(self, solvent, solute, interface, enderface):
        
        self._solvent   = None
        self._solute    = None
        self._interface = None
        self._enderface = None
        
        # assign variables
        if solvent is not None:
            self._solvent   = solvent.copy()
        if solute is not None:
            self._solute    = [ii.copy() for ii in as_list(solute)]
        if interface is not None:
            self._interface = interface.copy()
        if enderface is not None:
            self._enderface = enderface.copy() # use this to make a good sandwich!
        
        return
    
    def _make_sandwich(self):
        
        if self._interface is not None:
            for atom in self._interface._stype:
                atom.set_label(atom.label + "_i")
                
        if self._enderface is not None:
            for atom in self._enderface._stype:
                atom.set_label(atom.label + "_e")
                
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
    

    # populate a box with solvent and ions
    @staticmethod
    def populate_box(volume, instructions, input_file="input_packmol.in",
                     output_file="system.pdb"):

        if not instructions:
            return None

        # check volume
        assert len(volume) == 3, "Check volume!"
        
        # generate box boundaries with 1 AA padding
        box = np.concatenate(([1,1,1], np.asarray(volume)-1)).tolist()
        
        tmp_files = ["packmol.log", "input_packmol.in", "system.pdb"]
        with open(input_file, "w") as fout:
            
            fout.write(header.format(output_file, np.random.randint(100000)))
            
            for cc, instruction in enumerate(instructions):
                
                # unpack instructions
                mol = instruction[0]
                rep = instruction[1]
                typ = instruction[2]
                
                if isinstance(rep, int):
                    if not rep:
                        continue
                
                if typ == "box": # normal add
                    fout.write(box_place.format(cc, rep, " ".join(map(str, box))))
                
                elif typ == "fixed": # coordinate -> fixed point
                    
                    fout.write(fix_place.format(cc, *rep))
                    # make temp box to fill
                    # tbox     = box.copy()
                    # tbox[2]  = rep - 1
                    # tbox[-1] = rep + 1
                
                    # fout.write(box_place.format(cc, 1, " ".join(map(str, tbox))))
                
                else:
                    raise "Wrong instructions"
                
                # write tmp pdb file and store info
                mol.atoms.write("mol_{}.pdb".format(cc))
                tmp_files.append("mol_{}.pdb".format(cc))
                
        # run packmol
        try:
            subprocess.run(['packmol < {} > packmol.log'.format(input_file)],
                           shell=True, check=True, text=True)

        except:
            print("WARNING: packmol might not have worked, check system.")
        
        try:
            universe = mda.Universe(output_file)
        except:
            universe = None

        # remove temp mol files and packmol files
        subprocess.call(['rm'] + tmp_files)
        
        return universe

    # generate a slab from a unit cell
    @staticmethod
    def make_interface_slab(interface_uc, xsize, ysize, layers=1):
        
        if layers == 0 or interface_uc is None:
            return None
        
        xrep = int(np.round(xsize/interface_uc.atoms.get_cell()[0][0]))
        yrep = int(np.round(ysize/interface_uc.atoms.get_cell()[1][1]))
        
        slab = interface_uc.copy()
        
        if not np.isclose(np.dot(slab.atoms.cell[0], [1,0,0]), slab.atoms.cell[0][0]):
            xrep +=1
            print("WARNING: check interface if pattern matches")
        
        if not np.isclose(np.dot(slab.atoms.cell[1], [0,1,0]), slab.atoms.cell[1][1]):
            yrep +=1
            print("WARNING: check interface if pattern matches")
            
        slab.repeat((xrep, yrep, 1), make_cubic=True)
        
        if layers > 1: # helps with indexing
            slab.repeat([1,1,layers])
        
        slab.atoms.center()
        # slab.atoms.rattle()
        
        return slab
    
    @staticmethod #THANKS CHATGPT
    def populate_with_ions(ions, nions, volume, ion_pos=False):
        volume = np.array(volume)
        ion_coords = []
        instructions = []
    
        to_center = False
        
        if ion_pos == "left":
            volume[2] /= 2
        elif ion_pos == "center":
            if not isinstance(nions, int) and len(nions) != 1:
                raise ValueError("Center positioning can only be used with a single ion")
            to_center = True
        elif ion_pos == "box":
            for cc, ion in enumerate(ions):
                nrep = nions if isinstance(nions, int) else nions[cc]
                instructions.append((ion.to_universe(), nrep, "box"))
            return instructions
    
        for cc, ion in enumerate(ions):
            ion_radius = ion.estimate_sphere_radius()
            nrep = nions[cc] if isinstance(nions, list) else nions
    
            for _ in range(nrep):
                max_attempts = 100  # Limit placement attempts to avoid infinite loop
                for _ in range(max_attempts):
                    if to_center:
                        new_coord = volume / 2
                    else:
                        new_coord = ion_radius + 1 + np.random.rand(3) * (volume - 2 * (ion_radius + 1))
    
                    if ion_coords:
                        distances = np.linalg.norm(ion_coords - new_coord, axis=1)
                        if np.all(distances > 3):
                            break
                    else:
                        break
                else:
                    print(f"Warning: Failed to place ion {ion} after {max_attempts} attempts")
    
                ion_coords.append(new_coord)
                instructions.append((ion.to_universe(), new_coord, "fixed"))
    
        return instructions

    
    def make_solvent_box(self, solvent, ions, volume, density, nions, concentration,
                         conmodel, ion_pos):
        
        # make sure info is sound
        assert not( nions is not None and concentration is not None),\
            "'nions' and 'concentration' cannot both be not None"
        
        # convert concentration to number of ions
        if concentration is not None:
            nions = int(concentration*np.prod(volume)*units.mol/((units.m/10)**3))
        
        # define instructions for packmol
        instructions = []
        
        # populate according to continuum model
        if conmodel is not None:
            z_positions = conmodel.discretize_profile(volume)
            
            for cc, ion in enumerate(ions):
                for z_pos in z_positions[cc]:
                    
                    # radius = conmodel.species[cc].radius
                    instructions.append([ion, z_pos, "fixed"])
        
        # populate according to fixed number of ions
        elif nions is not None and ions is not None:
            
            ion_instr = self.populate_with_ions(ions, nions, volume, ion_pos=ion_pos)
            instructions.extend(ion_instr)
        
        # add solvent
        if solvent is not None:
            solvent_volume   = 1e-24*np.prod(volume)
            mass = solvent.atoms.masses.sum()
    
            # number of solvent molecules
            nummols = int(units.mol*density*(1.0/mass)*solvent_volume)
            
            instructions.append([solvent, nummols, "box"])
        
        # generate universe file
        universe = self.populate_box(volume, instructions)
        
        if universe is None:
            return None
        
        # Create a dictionary for quick lookup of species by residue name
        species_dict = {specie.residues.resnames[0]: specie for specie in self.species}
        
        alist = []
        for res in universe.residues:
            resname = res.resname
            if resname in species_dict:
                nmol = species_dict[resname].copy()
                nmol.atoms.positions = res.atoms.positions
                alist.append(nmol.atoms)
        
        solution = mda.Merge(*alist)
        solution.dimensions = volume + [90,90,90]
        
        return solution
    
    @staticmethod
    def _get_size_from_slab(slab):
        
        xsize = [1,0,0]@slab.atoms.cell@[1,0,0]
        ysize = [0,1,0]@slab.atoms.cell@[0,1,0]
        slab_depth = [0,0,1]@slab.atoms.cell@[0,0,1]
        
        
        return xsize, ysize, slab_depth
    
    def make_simulation_box(self, solvent_vol, solvent_rho, nions=None,
                            concentration=None, conmodel=None, layers=1,
                            padding=1.5, to_ase=False, mirror=False,
                            write_data=False, filename="data.lammps",
                            center_electrode=False, vacuum=None, layered=False,
                            ion_pos=None, hijack=None):
        
        # solvent volume
        xsize, ysize, zsize = solvent_vol
        
        # make slab
        islab = self.make_interface_slab(self._interface, xsize, ysize, layers=layers)
        eslab = self.make_interface_slab(self._enderface, xsize, ysize, layers=layers)
        
        xi, yi, sdi, xe, ye, sde = 0, 0, 0, 0, 0, 0
        # update the volume with multiples of UC
        if islab is not None:
            xi, yi, sdi = self._get_size_from_slab(islab)
            islab = islab.to_universe(layered=layered)
        if eslab is not None:
            xe, ye, sde = self._get_size_from_slab(eslab)
            eslab = eslab.to_universe(layered=layered)
            
        if eslab is not None and islab is not None: # check they have same size
            assert xi == xe
            assert yi == ye
        else:
            # slab_depth = 0
            padding    = 0
        
        if eslab is not None or islab is not None:
            xsize = np.maximum(xi, xe)
            ysize = np.maximum(yi, ye)
        
        # make solvent box
        solvent = self.make_solvent_box(self.solvent, self._solute,
                                        [xsize, ysize, zsize], solvent_rho,
                                        nions, concentration, conmodel, ion_pos)
        
        # now build system
        system = None
        zdim   = 0
        
        if islab is not None:
            system = islab
            zdim   += sdi
        
        if solvent is not None:
            if system is None:
                system = solvent
                zdim   += zsize
            else:
                solvent.atoms.translate([0, 0, sdi + padding])
                system = mda.Merge(system.atoms, solvent.atoms)
                zdim   += zsize + padding
        
        if eslab is not None:
            if system is None:
                system = eslab
                zdim   += sde
            else:
                if solvent is None:
                    eslab.atoms.translate([0, 0, zdim])
                    system = mda.Merge(system.atoms, eslab.atoms)
                    zdim   += sde
                else:
                    eslab.atoms.translate([0, 0, zdim + padding])
                    system = mda.Merge(system.atoms, eslab.atoms)
                    zdim   += sde + padding
                
        
        system.dimensions = [xsize, ysize, zdim] + [90, 90, 90] #TODO not like this
        
        if vacuum is not None:
            system.dimensions[2] += vacuum
            system.atoms.translate([0,0,+vacuum/2])
            zdim += vacuum
        
        if center_electrode:
            system.atoms.translate([0,0,zdim/2])
            _ = system.atoms.wrap()
            
        # give ase atoms to override positions
        if hijack is not None:
            system.dimensions = hijack.get_cell_lengths_and_angles()
            system.atoms.positions = hijack.get_positions()
        
        if write_data:
            self.write_lammps_file(system, filename=filename)
        
        if to_ase:
            return self.to_ase(system)
        
        return system
    
    def write_lammps_file(self, system,  write_coeff=True, filename="data.lammps"):
        
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
            
        if system.atoms.charges is not None:
            ase_system.set_initial_charges(system.atoms.charges)
        
        return ase_system

    @property
    def solvent(self):
        if self._solvent is None:
            return None
        return self._solvent.to_universe()
    
    @property
    def solute(self):
        if self._solute is None:
            return None
        return [ii.to_universe() for ii in self._solute]
    
    @property
    def interface(self):
        if self._interface is None:
            return None
        return self._interface.to_universe()
    
    @property
    def enderface(self):
        if self._enderface is None:
            return None
        return self._enderface.to_universe()
    
    @property
    def species(self):
        # merge all species in system
        all_species = np.concatenate((as_list(self.solvent), as_list(self.solute),
                                      as_list(self.interface), as_list(self.enderface)))
        return [ii for ii in all_species if ii is not None]
    
    @property
    def _species(self):
        
        return np.concatenate((as_list(self._solvent), as_list(self._solute),
                                      as_list(self._interface), as_list(self._enderface)))
    
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
    
