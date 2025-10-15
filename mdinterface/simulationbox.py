#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 15:14:41 2023

@author: roncoroni
"""

from typing import List, Dict, Union, Optional, Tuple, Any
from mdinterface.utils.auxiliary import label_to_element, as_list, find_smallest_missing
from mdinterface.io.lammpswriter import DATAWriter, write_lammps_coefficients
from mdinterface.build.box import make_interface_slab, make_solvent_box, add_component

import ase
import MDAnalysis as mda

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
        "nsolvent": None,  # Number of solvent molecules (alternative to rho)
        "zdim": None,
        "nspecies": None,  # Preferred parameter name
        "nions": None,     # Deprecated but maintained for backward compatibility
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
    
    def __init__(self, solvent: Optional[Any] = None, solute: Optional[Union[Any, List[Any]]] = None,
                 interface: Optional[Any] = None, enderface: Optional[Any] = None,
                 miderface: Optional[Any] = None) -> None:
        
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
                    
    @staticmethod
    def _get_size_from_slab(slab):
        
        xsize = [1,0,0]@slab.atoms.cell@[1,0,0]
        ysize = [0,1,0]@slab.atoms.cell@[0,1,0]
        slab_depth = [0,0,1]@slab.atoms.cell@[0,0,1]
        
        return xsize, ysize, slab_depth

    def _validate_xysize(self, xysize: Union[Tuple[float, float], List[float]]) -> Tuple[float, float]:
        """
        Validate and normalize XY size input.

        Parameters:
        -----------
        xysize : Union[Tuple[float, float], List[float]]
            Cross-sectional dimensions

        Returns:
        --------
        Tuple[float, float]
            Validated (xsize, ysize) tuple

        Raises:
        -------
        TypeError
            If xysize is not a list, tuple, or array
        ValueError
            If xysize doesn't have exactly 2 elements or contains non-positive values
        """
        if not isinstance(xysize, (list, tuple, np.ndarray)):
            raise TypeError(
                f"xysize must be a list, tuple, or numpy array, got {type(xysize).__name__}. "
                f"Expected format: [x_size, y_size] where both values are positive numbers."
            )

        if len(xysize) != 2:
            raise ValueError(
                f"xysize must have exactly 2 elements [x_size, y_size], got {len(xysize)} elements: {xysize}. "
                f"Example: xysize=[20.0, 20.0] for a 20x20 Angstrom cross-section."
            )

        try:
            xsize, ysize = float(xysize[0]), float(xysize[1])
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"xysize elements must be numeric, got {xysize}. "
                f"Both x_size and y_size must be convertible to float. Error: {e}"
            )

        if xsize <= 0 or ysize <= 0:
            raise ValueError(
                f"xysize values must be positive, got x_size={xsize}, y_size={ysize}. "
                f"Both dimensions must be greater than 0 Angstroms."
            )

        return xsize, ysize

    def _validate_layering(self, layering: List[Dict[str, Any]]) -> None:
        """
        Validate layering configuration before processing.

        Parameters:
        -----------
        layering : List[Dict[str, Any]]
            List of layer dictionaries

        Raises:
        -------
        TypeError
            If layering is not a list or contains non-dictionary elements
        ValueError
            If layer configurations are invalid
        """
        if not isinstance(layering, list):
            raise TypeError(
                f"layering must be a list of layer dictionaries, got {type(layering).__name__}. "
                f"Expected format: [{{\"type\": \"solvent\", ...}}, {{\"type\": \"interface\", ...}}]"
            )

        if not layering:
            raise ValueError(
                "layering cannot be empty. At least one layer must be specified. "
                "Example: layering=[{\"type\": \"solvent\", \"zdim\": 25, \"rho\": 1.0}]"
            )

        valid_types = {"interface", "enderface", "miderface", "solvent", "vacuum"}

        for i, layer in enumerate(layering):
            if not isinstance(layer, dict):
                raise TypeError(
                    f"Layer {i} must be a dictionary, got {type(layer).__name__}: {layer}. "
                    f"Each layer must have a 'type' key and appropriate parameters."
                )

            if "type" not in layer:
                raise ValueError(
                    f"Layer {i} missing required 'type' key: {layer}. "
                    f"Valid types: {valid_types}"
                )

            layer_type = layer["type"]
            if layer_type not in valid_types:
                raise ValueError(
                    f"Layer {i} has invalid type '{layer_type}'. "
                    f"Valid types: {valid_types}. "
                    f"Check spelling and ensure the layer type is supported."
                )

            # Type-specific validation
            if layer_type == "solvent":
                required_keys = {"zdim"}
                missing_keys = required_keys - layer.keys()
                if missing_keys:
                    raise ValueError(
                        f"Solvent layer {i} missing required keys: {missing_keys}. "
                        f"Solvent layers must specify 'zdim' (thickness in Angstroms). "
                        f"Example: {{\"type\": \"solvent\", \"zdim\": 25, \"rho\": 1.0}}"
                    )

                # Validate rho vs nsolvent mutual exclusivity
                rho = layer.get("rho")
                nsolvent = layer.get("nsolvent")
                if rho is not None and nsolvent is not None:
                    import warnings
                    warnings.warn(
                        f"Solvent layer {i}: Both 'rho' (density) and 'nsolvent' (number of molecules) are specified. "
                        f"Using 'nsolvent' and ignoring 'rho'. "
                        f"To avoid this warning, specify only one parameter.",
                        UserWarning, stacklevel=3
                    )
                elif rho is None and nsolvent is None:
                    raise ValueError(
                        f"Solvent layer {i} must specify either 'rho' (density in g/cm³) or 'nsolvent' (number of molecules). "
                        f"Examples: {{\"rho\": 1.0}} or {{\"nsolvent\": 500}}"
                    )


            elif layer_type in ["interface", "enderface", "miderface"]:
                if "nlayers" in layer and not isinstance(layer["nlayers"], int):
                    raise ValueError(
                        f"{layer_type.capitalize()} layer {i} 'nlayers' must be an integer, "
                        f"got {type(layer['nlayers']).__name__}: {layer['nlayers']}"
                    )

    def make_simulation_box(
        self,
        xysize: Union[Tuple[float, float], List[float]],
        layering: List[Dict[str, Any]],
        padding: float = 1.5,
        to_ase: bool = False,
        write_data: bool = False,
        filename: str = "data.lammps",
        center_electrode: bool = False,
        layered: bool = False,
        hijack: Optional[ase.Atoms] = None,
        match_cell: bool = False,
        atom_style: str = "full",
        write_coeff: bool = True
    ) -> Union[mda.Universe, ase.Atoms]:
        """
        Generate a simulation box from layered components.

        This is the main driver function that assembles a complete simulation system
        by stacking layers (interfaces, solvents, vacuum) according to the provided
        configuration.

        Parameters:
        -----------
        xysize : Tuple[float, float] or List[float]
            Cross-sectional dimensions [x_size, y_size] in Angstroms
        layering : List[Dict[str, Any]]
            List of layer configurations. Each layer must have a 'type' key.
            Valid types: 'interface', 'enderface', 'miderface', 'solvent', 'vacuum'
        padding : float, default=1.5
            Spacing between layers in Angstroms
        to_ase : bool, default=False
            If True, return ase.Atoms object; otherwise return MDAnalysis.Universe
        write_data : bool, default=False
            If True, write LAMMPS data file
        filename : str, default="data.lammps"
            Name of the LAMMPS data file to write
        center_electrode : bool, default=False
            If True, shift system by 50% along Z to center the first interface
        layered : bool, default=False
            Assign different molecule indexes to each interface layer for LAMMPS
        hijack : Optional[ase.Atoms], default=None
            Override positions with provided ase.Atoms object (use with caution)
        match_cell : bool, default=False
            Deform slabs to match XY dimensions (use with care)
        atom_style : str, default="full"
            LAMMPS atom style format ('full' or 'atomic')
        write_coeff : bool, default=True
            Whether to write force field coefficients in LAMMPS data file

        Returns:
        --------
        Union[mda.Universe, ase.Atoms]
            The assembled simulation system

        Raises:
        -------
        TypeError
            If xysize or layering have incorrect types
        ValueError
            If xysize dimensions are invalid or layering configuration is malformed
        """
        
        # Validate and parse input parameters
        xsize, ysize = self._validate_xysize(xysize)
        self._validate_layering(layering)
        
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
                nsolvent = layer["nsolvent"]

                # Handle backward compatibility: nions -> nspecies
                nspecies = layer.get("nspecies")
                nions_legacy = layer.get("nions")

                if nspecies is not None and nions_legacy is not None:
                    raise ValueError("Cannot specify both 'nspecies' and 'nions'. Use 'nspecies' (preferred) or 'nions' (deprecated).")
                elif nspecies is not None:
                    nions = nspecies  # Use the preferred parameter
                elif nions_legacy is not None:
                    nions = nions_legacy  # Use legacy parameter for backward compatibility
                    import warnings
                    warnings.warn("Parameter 'nions' is deprecated. Use 'nspecies' instead.",
                                DeprecationWarning, stacklevel=2)
                else:
                    nions = None

                concentration = layer["concentration"] # in molar
                conmodel = layer["conmodel"]
                ion_pos  = layer["ion_pos"]

                # make solvent box
                solvent = make_solvent_box(self.species, self.solvent, self._solute,
                                           [xsize, ysize, zsize], solv_rho,
                                           nions, concentration, conmodel, ion_pos, nsolvent)

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
            self.write_lammps_file(system, filename=filename,
                                   atom_style=atom_style, write_coeff=write_coeff)
        
        # convert to ase, or not
        if to_ase:
            return self.to_ase(system)
        return system
    
    def write_lammps_file(self, system: mda.Universe, write_coeff: bool = True,
                          filename: str = "data.lammps", atom_style: str = "full") -> None:
        
        # just make sure we are not messing things up
        system = system.copy()
        
        # remove coefficients
        if not write_coeff:
            for attribute in ["bonds", "angles", "dihedrals", "impropers"]:
                try:
                    system.del_TopologyAttr(attribute)
                except:
                    pass
        
        # first write data file
        with DATAWriter(filename) as dt:
            dt.write(system.atoms, atom_style=atom_style)
            
        # now write coeff where they belong
        if write_coeff:
            
            temp_file = 'tmp_data.lammps'
            
            with open(filename, 'r') as ffile, open(temp_file, 'w') as tfile:
                for ln, fl in enumerate(ffile):
                    if fl.startswith("Atoms"):

                        # write coefficients using refactored function
                        sorted_attrs = {
                            "atoms": self.get_sorted_attribute("atoms"),
                            "bonds": self.get_sorted_attribute("bonds"),
                            "angles": self.get_sorted_attribute("angles"),
                            "dihedrals": self.get_sorted_attribute("dihedrals"),
                            "impropers": self.get_sorted_attribute("impropers")
                        }
                        write_lammps_coefficients(system, sorted_attrs, fout=tfile)
                        
                    tfile.write(fl)
            
            shutil.move(temp_file, filename)
        
        return

    # convert to ase.Atoms
    @staticmethod
    def to_ase(system: Optional[mda.Universe]) -> ase.Atoms:
        
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
        
        xsize_t, ysize_t = np.nan, np.nan
        
        for tslab in slabs:
            xi, yi, _ = self._get_size_from_slab(tslab)
            xsize_t = np.nanmax([xi, xsize_t])
            ysize_t = np.nanmax([yi, ysize_t])
        
        return xsize_t, ysize_t
