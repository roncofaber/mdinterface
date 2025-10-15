#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 15:03:47 2023

@author: roncoroni
"""

from typing import Dict, List, Optional, Union, TextIO, Any
import numpy as np

from MDAnalysis.lib import util, mdamath
from MDAnalysis.core.groups import requires
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.coordinates import base
import MDAnalysis as mda


btype_sections = {'bond':'Bonds', 'angle':'Angles',
                  'dihedral':'Dihedrals', 'improper':'Impropers'}


class DATAWriter(base.WriterBase):
    """Write out the current time step as a LAMMPS DATA file.

    This writer supports the sections Atoms, Masses, Velocities, Bonds,
    Angles, Dihedrals, and Impropers. This writer will write the header
    and these sections (if applicable). Atoms section is written in the
    "full" sub-style if charges are available or "molecular" sub-style
    if they are not. Molecule id is set to 0 for all atoms.

    Note
    ----
    This writer assumes "conventional" or "real" LAMMPS units where length
    is measured in Angstroms and velocity is measured in Angstroms per
    femtosecond. To write in different units, specify `lengthunit`

    If atom types are not already positive integers, the user must set them
    to be positive integers, because the writer will not automatically
    assign new types.

    To preserve numerical atom types when writing a selection, the Masses
    section will have entries for each atom type up to the maximum atom type.
    If the universe does not contain atoms of some type in
    {1, ... max(atom_types)}, then the mass for that type will be set to 1.

    In order to write bonds, each selected bond type must be explicitly set to
    an integer >= 1.

    """
    format = 'DATA'

    def __init__(self, filename, convert_units=True, **kwargs):
        """Set up a DATAWriter

        Parameters
        ----------
        filename : str
            output filename
        convert_units : bool, optional
            units are converted to the MDAnalysis base format; [``True``]
        """
        self.filename = util.filename(filename, ext='data', keep=True)

        self.convert_units = convert_units

        self.units = {'time': 'fs', 'length': 'Angstrom'}
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['velocity'] = kwargs.pop('velocityunit',
                                 self.units['length']+'/'+self.units['time'])

    def _write_atoms(self, atoms, data, atom_style):
        self.f.write('\n')
        self.f.write('Atoms\n')
        self.f.write('\n')

        try:
            charges = atoms.charges
            has_charges = True
        except (NoDataError, AttributeError):
            has_charges = False

        indices = atoms.indices + 1
        i_l = len(str(indices.max()))
        
        unique_types, types = np.unique(atoms.types, return_inverse=True)
        types += 1
        t_l = len(str(len(unique_types)))

        moltags = atoms.resindices
        m_l = len(str(moltags.max()))+1
    
        if self.convert_units:
            coordinates = self.convert_pos_to_native(atoms.positions, inplace=False)
            
        b_l = len(str(int(coordinates.max()))) + 1
        b_t = b_l + 6

        if atom_style == "atomic":
            for index, atype, coords in zip(indices, types, coordinates):
                x, y, z = coords
                self.f.write(f"{index:{i_l}d}  {atype:{t_l}d}"
                             f"  {x:> {b_t}.5f} {y:> {b_t}.5f} {z:> {b_t}.5f}\n")
        elif atom_style == "full":
            
            if not has_charges:
                for index, moltag, atype, coords in zip(indices, moltags, types,
                        coordinates):
                    x, y, z = coords
                    self.f.write(f"{index:{i_l}d} {moltag:{m_l}d} {atype:{t_l}d}"
                                 f" {x:> {b_t}.5f} {y:> {b_t}.5f} {z:> {b_t}.5f}\n")
            elif has_charges:
                for index, moltag, atype, charge, coords in zip(indices, moltags,
                                                                types, charges, coordinates):
                    x, y, z = coords
                    self.f.write(f"{index:{i_l}d} {moltag:{m_l}d} {atype:{t_l}d}"
                                 f"  {charge: .7f}"
                                 f" {x:> {b_t}.5f} {y:> {b_t}.5f} {z:> {b_t}.5f}\n")

    def _write_velocities(self, atoms):
        self.f.write('\n')
        self.f.write('Velocities\n')
        self.f.write('\n')
        indices = atoms.indices + 1
        velocities = self.convert_velocities_to_native(atoms.velocities,
                                                       inplace=False)
        for index, vel in zip(indices, velocities):
            self.f.write('{i:d} {x:f} {y:f} {z:f}\n'.format(i=index, x=vel[0],
                y=vel[1], z=vel[2]))

    def _write_masses(self, atoms):
        # self.f.write('\n')
        self.f.write('Masses\n')
        self.f.write('\n')
        
        for cc, atype in enumerate(np.unique(atoms.types)):
            
            # search entire universe for mass info, not just writing selection
            masses = set(atoms.select_atoms('type {}'.format(atype)).masses)
            
            if len(masses) == 0:   
                mass = 1.0  
            else:
                mass = masses.pop()
                
            if masses:
                raise ValueError('LAMMPS DATAWriter: to write data file, '+
                        'atoms with same type must have same mass')
        
        
            self.f.write('{:d} {:> 8.3f}   # {}\n'.format(cc+1, mass, atype))

    def _write_bonds(self, bonds):
        self.f.write('\n')
        self.f.write('{}\n'.format(btype_sections[bonds.btype]))
        self.f.write('\n')
        for bond, i in zip(bonds, range(1, len(bonds)+1)):
            try:
                self.f.write('{:d} {:d} '.format(i, int(bond.type))+\
                        ' '.join((bond.atoms.indices + 1).astype(str))+'\n')
            except TypeError:
                errmsg = (f"LAMMPS DATAWriter: Trying to write bond, but bond "
                          f"type {bond.type} is not numerical.")
                raise TypeError(errmsg) from None

    def _write_dimensions(self, dimensions):
        """Convert dimensions to triclinic vectors, convert lengths to native
        units and then write the dimensions section
        """
        if self.convert_units:
            triv = self.convert_pos_to_native(mdamath.triclinic_vectors(
                                              dimensions),inplace=False)
        self.f.write('\n')
        self.f.write('{:f} {:f} xlo xhi\n'.format(0., triv[0][0]))
        self.f.write('{:f} {:f} ylo yhi\n'.format(0., triv[1][1]))
        self.f.write('{:f} {:f} zlo zhi\n'.format(0., triv[2][2]))
        if any([triv[1][0], triv[2][0], triv[2][1]]):
            self.f.write('{xy:f} {xz:f} {yz:f} xy xz yz\n'.format(
                xy=triv[1][0], xz=triv[2][0], yz=triv[2][1]))
        self.f.write('\n')

    @requires('types', 'masses')
    def write(self, selection, frame=None, atom_style='full'):
        """Write selection at current trajectory frame to file.

        The sections for Atoms, Masses, Velocities, Bonds, Angles,
        Dihedrals, and Impropers (if these are defined) are
        written. The Atoms section is written in the "full" sub-style
        if charges are available or "molecular" sub-style if they are
        not. Molecule id in atoms section is set to to 0.

        No other sections are written to the DATA file.
        As of this writing, other sections are not parsed into the topology
        by the :class:`DATAReader`.

        Note
        ----
        If the selection includes a partial fragment, then only the bonds,
        angles, etc. whose atoms are contained within the selection will be
        included.

        Parameters
        ----------
        selection : AtomGroup or Universe
            MDAnalysis AtomGroup (selection or Universe.atoms) or also Universe
        frame : int (optional)
            optionally move to frame number `frame`

        """
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]
        else:
            frame = u.trajectory.ts.frame

        # make sure to use atoms (Issue 46)
        atoms = selection.atoms

        # check that types can be converted to ints if they aren't ints already
        # try:
        #     atoms.types.astype(np.int32)
        # except ValueError:
        #     errmsg = ("LAMMPS.DATAWriter: atom types must be convertible to "
        #               "integers")
        #     raise ValueError(errmsg) from None

        try:
            velocities = atoms.velocities
        except (NoDataError, AttributeError):
            has_velocities = False
        else:
            has_velocities = True

        features = {}
        with util.openany(self.filename, 'wt') as self.f:
            self.f.write('LAMMPS data file via MDAnalysis and mdinterface\n')
            self.f.write('\n')
            self.f.write('{:>12d}  atoms\n'.format(len(atoms)))

            attrs = [('bond', 'bonds'), ('angle', 'angles'),
                ('dihedral', 'dihedrals'), ('improper', 'impropers')]

            for btype, attr_name in attrs:
                try:
                    features[btype] = atoms.__getattribute__(attr_name)
                    self.f.write('{:>12d}  {}\n'.format(len(features[btype]),
                                                        attr_name))
                    features[btype] = features[btype].atomgroup_intersection(
                                        atoms, strict=True)
                except:
                    pass

            self.f.write('\n')
            
            ntypes = len(set(atoms.types))
            self.f.write('{:>12d}  atom types\n'.format(ntypes))

            for btype, attr in features.items():
                self.f.write('{:>12d}  {} types\n'.format(len(attr.types()),
                                                          btype))

            self._write_dimensions(atoms.dimensions)

            self._write_masses(atoms)
            self._write_atoms(atoms, u.trajectory.ts.data, atom_style)
            for attr in features.values():
                if attr is None or len(attr) == 0:
                    continue
                self._write_bonds(attr)

            if has_velocities:
                self._write_velocities(atoms)


def write_lammps_coefficients(
    system: mda.Universe,
    sorted_attributes: Dict[str, List[Any]],
    fout: Optional[TextIO] = None,
    fname: str = "coefficients.lammps"
) -> None:
    """
    Write LAMMPS force field coefficients to a file.

    This function writes pair, bond, angle, dihedral, and improper coefficients
    in LAMMPS format. It's designed to be used with SimulationBox systems that
    contain force field parameters.

    Parameters:
    -----------
    system : MDAnalysis.Universe
        The system containing atoms and topology for which coefficients are written
    sorted_attributes : dict
        Dictionary containing sorted attributes from SimulationBox.get_sorted_attribute():
        - "atoms": sorted atom types
        - "bonds": sorted bond types
        - "angles": sorted angle types
        - "dihedrals": sorted dihedral types
        - "impropers": sorted improper types
    fout : file object, optional
        Open file object to write to. If None, creates new file with `fname`
    fname : str, optional
        Filename to create if `fout` is None. Default: "coefficients.lammps"

    Returns:
    --------
    None

    Examples:
    ---------
    # Basic usage with SimulationBox
    sorted_attrs = {
        "atoms": simbox.get_sorted_attribute("atoms"),
        "bonds": simbox.get_sorted_attribute("bonds"),
        "angles": simbox.get_sorted_attribute("angles"),
        "dihedrals": simbox.get_sorted_attribute("dihedrals"),
        "impropers": simbox.get_sorted_attribute("impropers")
    }
    write_lammps_coefficients(system, sorted_attrs, fname="my_coeffs.lammps")

    # Writing to an existing file object
    with open("data.lammps", "w") as f:
        write_lammps_coefficients(system, sorted_attrs, fout=f)
    """

    remember_to_close = False
    if fout is None:
        fout = open(fname, "w")
        remember_to_close = True

    # Write Pair Coefficients
    fout.write("Pair Coeffs\n\n")

    idx = 1
    atom_types = sorted_attributes.get("atoms", [])
    for atom in atom_types:
        if atom.extended_label not in np.unique(system.atoms.types):
            continue

        eps = atom.eps if atom.eps is not None else 0
        sig = atom.sig if atom.sig is not None else 0

        fout.write("{:>5}    {:>12.8f}    {:>12.8f}  # {}\n".format(
            idx, eps, sig, atom.extended_label))
        idx += 1

    # Write Bond Coefficients
    bond_types = sorted_attributes.get("bonds", [])
    if bond_types:
        fout.write("\n")
        fout.write("Bond Coeffs\n\n")

    for bond in bond_types:
        if bond.id not in np.array(system.bonds.types(), dtype=int):
            continue

        kr = bond.kr if bond.kr is not None else 0
        r0 = bond.r0 if bond.r0 is not None else 0

        btype = "{}-{}".format(*bond.symbols)

        fout.write("{:>5}    {:>10.6f}    {:>10.6f}  #  {:<5} | {}\n".format(
            bond.id, kr, r0, btype, bond.resname))

    # Write Angle Coefficients
    angle_types = sorted_attributes.get("angles", [])
    if angle_types:
        fout.write("\n")
        fout.write("Angle Coeffs\n\n")

    for angle in angle_types:
        if angle.id not in np.array(system.angles.types(), dtype=int):
            continue

        kr = angle.kr if angle.kr is not None else 0
        theta0 = angle.theta0 if angle.theta0 is not None else 0

        atype = "{}-{}-{}".format(*angle.symbols)

        fout.write("{:>5}    {:>10.6f}    {:>10.6f}  #  {:<8} | {}\n".format(
            angle.id, kr, theta0, atype, angle.resname))

    # Write Dihedral Coefficients
    dihedral_types = sorted_attributes.get("dihedrals", [])
    if dihedral_types:
        fout.write("\n")
        fout.write("Dihedral Coeffs\n\n")

    for dihedral in dihedral_types:
        if dihedral.id not in np.array(system.dihedrals.types(), dtype=int):
            continue

        dihedral.write(fout)

    # Write Improper Coefficients
    improper_types = sorted_attributes.get("impropers", [])
    if improper_types and system.impropers.types():
        fout.write("\n")
        fout.write("Improper Coeffs\n\n")
        
    for improper in improper_types:
        if improper.id not in np.array(system.impropers.types(), dtype=int):
            continue

        atype = "{}".format(*improper.symbols)
        value = "{:>7.4f}    {:>2d}    {:>2d}".format(*improper.values)

        fout.write("{:>5}    {}  #  {:<2} | {}\n".format(improper.id, value, atype, improper.resname))

    fout.write("\n")

    if remember_to_close:
        fout.close()
