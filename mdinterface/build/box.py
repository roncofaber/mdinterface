#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 10:23:57 2025

@author: roncofaber
"""

from typing import List, Optional, Union, Tuple, Dict, Any
from mdinterface.io.packmol import header, box_place, fix_place
from mdinterface.build.continuum2sim import discretize_concentration

import MDAnalysis as mda
from ase import units

import numpy as np

import subprocess

#%%

def _validate_solvent_box_parameters(
    nions: Optional[Union[int, List[int]]],
    concentration: Optional[float],
    conmodel: Optional[Dict[int, Tuple[List[float], List[float]]]],
    ions: Optional[List[Any]],
    solvent: Optional[Any],
    density: Optional[float],
    nsolvent: Optional[int] = None
) -> None:
    """
    Validate parameter combinations for make_solvent_box function.

    Raises appropriate errors for invalid parameter combinations.
    """

    # Basic mutual exclusivity check
    if nions is not None and concentration is not None:
        raise ValueError("Cannot specify both 'nions' and 'concentration'. Use one or the other.")

    # If concentration model is provided, ions must be provided
    if conmodel is not None and (ions is None or len(ions) == 0):
        raise ValueError("When using 'conmodel', 'ions' must be provided and non-empty.")

    # If nions is a list, it must match the number of ion species
    if isinstance(nions, (list, tuple)) and ions is not None:
        if len(nions) != len(ions):
            raise ValueError(f"Length of 'nions' ({len(nions)}) must match number of ion species ({len(ions)}).")

    # Validate density vs nsolvent mutual exclusivity
    if density is not None and nsolvent is not None:
        import warnings
        warnings.warn(
            "Both 'density' and 'nsolvent' are specified. Using 'nsolvent' and ignoring 'density'.",
            UserWarning, stacklevel=4
        )

    # If solvent density or nsolvent is provided but no solvent, warn the user
    if (density is not None or nsolvent is not None) and solvent is None:
        import warnings
        warnings.warn("Density or nsolvent specified but no solvent provided. Will be ignored.",
                     UserWarning, stacklevel=4)

    # If no solvent and no ions, nothing to do
    if solvent is None and (ions is None or len(ions) == 0):
        import warnings
        warnings.warn("No solvent or ions specified. Empty box will be created.",
                     UserWarning, stacklevel=3)

def make_solvent_box(
    species: List[Any],
    solvent: Optional[Any],
    ions: Optional[List[Any]],
    volume: List[float],
    density: Optional[float],
    nions: Optional[Union[int, List[int]]],
    concentration: Optional[float],
    conmodel: Optional[Dict[int, Tuple[List[float], List[float]]]],
    ion_pos: Optional[str],
    nsolvent: Optional[int] = None
) -> Optional[mda.Universe]:
    """
    Build a solvent box with optional ionic species.

    This function creates a simulation box containing solvent molecules and ionic species
    using PACKMOL for molecular packing. It supports various placement strategies and
    concentration models.

    Parameters:
    -----------
    species : list
        List of all available species in the simulation
    solvent : object or None
        Solvent molecule object (e.g., Water). If None, only ions are placed.
    ions : list or None
        List of ionic species to add to the box
    volume : list
        Box dimensions [x, y, z] in Angstroms
    density : float or None
        Solvent density in g/cm³. Ignored if solvent is None or if nsolvent is specified.
    nions : int, list, or None
        Number of each ionic species. Can be:
        - int: Same number for all ion types
        - list: Different number for each ion type (must match len(ions))
        - None: No ions added
    concentration : float or None
        Ionic concentration in Molar. Alternative to nions.
        Cannot be used simultaneously with nions.
    conmodel : dict or None
        Advanced concentration model for spatially varying concentrations.
        Format: {ion_index: (z_coords, concentration_profile)}
    ion_pos : str or None
        Ion placement strategy:
        - "random": Random placement (default)
        - "center": Place all ions at box center
        - "box": Use PACKMOL box placement
        - "left": Constrain to left half of box
        - None: Defaults to "random"
    nsolvent : int or None
        Number of solvent molecules to place. If specified, takes precedence over density.
        Cannot be used simultaneously with density.

    Returns:
    --------
    MDAnalysis.Universe or None
        Merged universe containing solvent and ions, or None if no components

    Examples:
    ---------
    # Simple water box with NaCl
    make_solvent_box(species, water, [na, cl], [20, 20, 20], 1.0, [5, 5], None, None, "random")

    # Concentration-based approach
    make_solvent_box(species, water, [na, cl], [20, 20, 20], 1.0, None, 0.1, None, "random")

    # Complex polymer solution
    make_solvent_box(species, None, [polymer, hydronium, water], [50, 50, 50], None, [10, 50, 200], None, None, "box")
    """

    # Validate parameter combinations
    _validate_solvent_box_parameters(nions, concentration, conmodel, ions, solvent, density, nsolvent)

    # Legacy parameter compatibility is handled in the calling function (simulationbox.py)

    # convert concentration to number of ions
    if concentration is not None:
        nions = int(concentration*np.prod(volume)*units.mol/((units.m/10)**3))
    
    # define instructions for packmol
    instructions = []
     
    # populate ions
    if (conmodel is not None) or (nions is not None and ions is not None):
        ion_instr = populate_with_ions(ions, nions, volume, ion_pos=ion_pos,
                                       conmodel=conmodel)
        instructions.extend(ion_instr)
    
    # add solvent
    if solvent is not None:
        if nsolvent is not None:
            # Use directly specified number of solvent molecules
            nummols = nsolvent
        elif density is not None:
            # Calculate number of solvent molecules from density
            solvent_volume = 1e-24*np.prod(volume)
            mass = solvent.atoms.masses.sum()
            nummols = int(units.mol*density*(1.0/mass)*solvent_volume)
        else:
            # This should be caught by validation, but defensive programming
            raise ValueError("Either 'density' or 'nsolvent' must be specified for solvent placement")

        instructions.append([solvent, nummols, "box"])

    # generate universe file
    universe = populate_box(volume, instructions)
    
    if universe is None:
        return None
    
    # Create a dictionary for quick lookup of species by residue name
    species_dict = {specie.residues.resnames[0]: specie for specie in species}
    
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
    
def populate_box(
    volume: List[float],
    instructions: List[Tuple[Any, Union[int, List[float]], str]],
    input_file: str = "input_packmol.in",
    output_file: str = "system.pdb"
) -> Optional[mda.Universe]:

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
                
            elif typ == "zfixed": # small bin to use in CM
                tbox = box.copy()
                tbox[2] = tbox[2] - mol.estimate_specie_radius()
                tbox[5] = tbox[5] + mol.estimate_specie_radius()
                fout.write(box_place.format(cc, 1, " ".join(map(str, tbox))))
                mol = mol.to_universe()
            
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

#THANKS CHATGPT (but mostly me tbh)
def populate_with_ions(ions, nions, volume, ion_pos=False, conmodel=None):
    
    def place_ion(ion, volume, ion_coords, ion_radii, zpos=None, max_attempts=100):
        ion_radius = ion.estimate_specie_radius()
        for _ in range(max_attempts):
            new_coord = ion_radius + 1 + np.random.rand(3) * (volume - 2 * (ion_radius + 1))
            if zpos is not None:
                new_coord[2] = zpos
                
            if not ion_coords or np.all(np.linalg.norm(ion_coords - new_coord, axis=1) >= np.array(ion_radii) + ion_radius + 1):
                return new_coord
        print(f"Warning: Failed to place ion {ion} after {max_attempts} attempts")
        return None

    def place_ions_conmodel(ions, conmodel, volume, max_attempts=100):
        instructions = []
        ion_coords = []
        ion_radii = []
        for cc, ion in enumerate(ions):
            z_coords, conc_profile = conmodel[cc]
            z_positions = discretize_concentration(ion, conc_profile, z_coords, volume)
            for zpos in z_positions:
                new_coord = place_ion(ion, volume, ion_coords, ion_radii, zpos=zpos,
                                      max_attempts=max_attempts)
                if new_coord is not None:
                    ion_coords.append(new_coord)
                    ion_radii.append(ion.estimate_specie_radius())
                    instructions.append((ion.to_universe(), new_coord, "fixed"))
        return instructions

    def place_ions_random(ions, nions, volume, to_center, max_attempts=100):
        instructions = []
        ion_coords = []
        ion_radii = []
        for cc, ion in enumerate(ions):
            nrep = nions if isinstance(nions, int) else nions[cc]
            for _ in range(nrep):
                new_coord = volume / 2 if to_center else place_ion(ion, volume,
                                                                   ion_coords,
                                                                   ion_radii,
                                                                   max_attempts=max_attempts)
                if new_coord is not None:
                    ion_coords.append(new_coord)
                    ion_radii.append(ion.estimate_specie_radius())
                    instructions.append((ion.to_universe(), new_coord, "fixed"))
        return instructions

    max_attempts = 100  # Limit placement attempts to avoid infinite loop
    volume = np.array(volume)

    if conmodel is not None:
        assert len(conmodel) == len(ions), "Need one profile per specie"
        return place_ions_conmodel(ions, conmodel, volume, max_attempts=max_attempts)

    if ion_pos == "box":
        return [(ion.to_universe(), nions if isinstance(nions, int) else nions[cc], "box") for cc, ion in enumerate(ions)]

    to_center = ion_pos == "center"
    if ion_pos == "left":
        volume[2] /= 2

    return place_ions_random(ions, nions, volume, to_center, max_attempts=max_attempts)


# add a component to the system
def add_component(system, component, zdim, padding=0):
    
    # nothing to add here
    if component is None:
        return system, zdim
    
    # ohh, let's lego the shit out of this
    component = component.copy()
    
    # component: "look at me, I am the system now."
    if system is None:
        component.atoms.translate([0, 0, zdim])
        system = component
        zdim += component.dimensions[2]
    
    # make space and add it to the pile
    else:
        component.atoms.translate([0, 0, zdim + padding])
        system = mda.Merge(system.atoms, component.atoms)
        zdim += component.dimensions[2] + padding
    return system, zdim
