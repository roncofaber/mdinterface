#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 11:55:55 2024

@author: roncofaber
"""

import networkx as nx
import numpy as np
import ase
from io import StringIO

from mdinterface.utils.auxiliary import label_to_element
from collections import deque
from ase.io.lammpsrun import construct_cell, lammps_data_to_ase_atoms
from ase import Atoms
#%%

def read_lammps_dump_text(fileobj, **kwargs):
    """Process cleartext lammps dumpfiles

    :param fileobj: filestream providing the trajectory data
    :param index: integer or slice object (default: get the last timestep)
    :returns: list of Atoms objects
    :rtype: list
    """
    # Load all dumped timesteps into memory simultaneously
    lines = deque(fileobj.readlines())

    n_atoms = 0

    # avoid references before assignment in case of incorrect file structure
    cell, celldisp, pbc = None, None, False

    while len(lines) > n_atoms:
        line = lines.popleft()

        if "ITEM: TIMESTEP" in line:
            n_atoms = 0
            line = lines.popleft()
            # !TODO: pyflakes complains about this line -> do something
            # ntimestep = int(line.split()[0])  # NOQA

        if "ITEM: NUMBER OF ATOMS" in line:
            line = lines.popleft()
            n_atoms = int(line.split()[0])

        if "ITEM: BOX BOUNDS" in line:
            # save labels behind "ITEM: BOX BOUNDS" in triclinic case
            # (>=lammps-7Jul09)
            tilt_items = line.split()[3:]
            celldatarows = [lines.popleft() for _ in range(3)]
            celldata = np.loadtxt(celldatarows)
            diagdisp = celldata[:, :2].reshape(6, 1).flatten()

            # determine cell tilt (triclinic case!)
            if len(celldata[0]) > 2:
                # for >=lammps-7Jul09 use labels behind "ITEM: BOX BOUNDS"
                # to assign tilt (vector) elements ...
                offdiag = celldata[:, 2]
                # ... otherwise assume default order in 3rd column
                # (if the latter was present)
                if len(tilt_items) >= 3:
                    sort_index = [tilt_items.index(i)
                                  for i in ["xy", "xz", "yz"]]
                    offdiag = offdiag[sort_index]
            else:
                offdiag = (0.0,) * 3

            cell, celldisp = construct_cell(diagdisp, offdiag)

            # Handle pbc conditions
            if len(tilt_items) == 3:
                pbc_items = tilt_items
            elif len(tilt_items) > 3:
                pbc_items = tilt_items[3:6]
            else:
                pbc_items = ["f", "f", "f"]
            pbc = ["p" in d.lower() for d in pbc_items]

        if "ITEM: ATOMS" in line:
            colnames = line.split()[2:]
            datarows = [lines.popleft() for _ in range(n_atoms)]
            data = np.loadtxt(datarows, dtype=str)
            
            out_atoms = lammps_data_to_ase_atoms(
                data=data,
                colnames=colnames,
                cell=cell,
                celldisp=celldisp,
                atomsobj=Atoms,
                pbc=pbc,
                **kwargs
            )
            
    return out_atoms


def traj2chunks(filename, every=1):
    """Generator to read trajectory file in chunks (memory efficient)."""
    with open(filename, "r") as fin:
        lines = []
        for cc, line in enumerate(fin):
            lines.append(line)
            if line.startswith("ITEM: TIMESTEP") and len(lines) > 1:
                
                if not cc%every:
                    yield "".join(lines[:-1])
                lines = [line]
                
        if lines:  # Handle remaining lines
            yield "".join(lines)

def dump2ase(lines, specorder=None):
    # Your processing logic for a single item goes here
    # atoms = ase.io.read(StringIO("".join(lines)), format="lammps-dump-text",
                        # parallel=False, specorder=specorder)
    
    atoms = read_lammps_dump_text(StringIO("".join(lines)), specorder=specorder)
    return atoms

# use networking algo to merge SeaUrchin coordination environments
def merge_clusters(cluster_idxs, tar_idxs=None):

    def to_edges(l):
        """
            treat `l` as a Graph and returns it's edges
            to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
        """
        it = iter(l)
        last = next(it)

        for current in it:
            yield last, current
            last = current

    clusters_graph = nx.Graph()

    if tar_idxs is not None:
        clusters_graph.add_nodes_from(tar_idxs)

    for idx in cluster_idxs:
        clusters_graph.add_edges_from(to_edges(idx))

    return [list(ii) for ii in nx.connected_components(clusters_graph)]

# # use a data.lammps file, BEST for lammps trajectories!
def read_data_file(dta_file, type2sym=False):

    def nonblank_lines(f):
        for l in f:
            line = l.rstrip()
            if line:
                yield line

    with open(dta_file) as fin:

        read_types, read_atoms, read_bonds = False, False, False
        labels, atoms, bonds, charges = [], [], [], []
        
        for line in nonblank_lines(fin):
            if line.endswith("atoms"):
                n_atoms = int(line.strip().split()[0])

            elif line.endswith("atom types"):
                n_types = int(line.strip().split()[0])

            elif line.startswith("Masses"):
                read_types = True
                continue

            elif line.startswith("Atoms"):
                read_atoms = True
                continue

            elif line.startswith("Bonds"):
                read_bonds = True
                continue

            elif line.startswith("Angles"):
                break

            elif read_types:
                parts = line.split()
                labels.append(label_to_element(parts[-1], float(parts[1])))
                if len(labels) >= n_types:
                    read_types = False

            elif read_atoms:
                parts = line.split()
                atoms.append(int(parts[2]))
                charges.append(float(parts[3]))
                if len(atoms) >= n_atoms:
                    read_atoms = False

            elif read_bonds:
                cline = line.split()
                bonds.append([int(cline[ii])-1 for ii in [2,3]])
                
    if not atoms:
        raise RuntimeError("No atoms found in file -- check again!")

    # make list of atomic symbols
    symbols = np.array([labels[ii-1] for ii in atoms])

    # use bonding info to calculate connectivity
    connectivity = [[] for _ in range(len(atoms))]
    for bond in bonds:
        for atom in bond:
            connectivity[atom].extend([ii for ii in bond if ii != atom])

    # remove duplicate entries
    connectivity = [list(set(con)) for con in connectivity]

    # generate atomic types depending on their environment (bonded to same elements)
    connected_elements = ["".join(sorted([symbols[ii] for ii in con]))
                          for con in connectivity]

    atom_types = np.unique(connected_elements)
    ato_typ = []
    for atom_type in connected_elements:
        ato_typ.append(np.where(atom_type == atom_types)[0][0])

    # create list of atoms belonging to same molecule
    list_of_molecules   = merge_clusters(bonds)

    # generate bonding data
    frame = ase.Atoms(symbols)

    # TODO soooooo messy but works thanks God

    # generate unique molecule ID number to recognize it later
    mol_idx = []
    counter = 0
    for cc in range(len(frame)):

        index = np.where([cc in ii for ii in list_of_molecules])[0].tolist()

        # not in a molecule
        if index == []:
            mol_idx.append(counter + len(list_of_molecules))
            counter += 1
        else:
            mol_idx.append(index[0])

    mol_idx = np.array(mol_idx)

    # generate the molecule names
    mol_names = []
    mol_typ = np.array(len(mol_idx)*[-1])
    for  mol_id in set(mol_idx):

        indexes = np.where(mol_idx == mol_id)[0]
        cmol = frame[indexes]

        chem_formula = cmol.get_chemical_formula()

        if chem_formula not in mol_names:
            mol_names.append(chem_formula)

        idnr = np.where([chem_formula == mol_type for mol_type in mol_names])[0]
        mol_typ[indexes] = idnr
        
    if type2sym:
        return labels , symbols, np.array(mol_idx), np.array(mol_typ), mol_names,\
            np.array(ato_typ), np.array(connectivity, dtype=object), np.array(charges)
        
    return symbols, np.array(mol_idx), np.array(mol_typ), mol_names,\
        np.array(ato_typ), np.array(connectivity, dtype=object), np.array(charges)
