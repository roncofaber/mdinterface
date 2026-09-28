#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Polymer snippet utilities for local topology refinement.

A snippet is a small sub-molecule cut from a polymer chain around a junction
atom, capped with terminating atoms (typically H), and passed to LigParGen to
obtain accurate OPLS-AA parameters for the chain interior.  The resulting
topology is then remapped back onto the full polymer indices.
"""

import numpy as np
import copy

def make_snippet(polymer, center, Nmax, ending="H", preserve_rings=True):

    from rdkit import Chem
    from mdinterface.core.chemistry import capped_molecule

    mol = polymer.to_rdkit()
    selected = {center}
    selected.update(np.concatenate(polymer.find_relevant_distances(Nmax, centers=center)).tolist())
    changed = True
    while changed:
        previous = set(selected)
        if preserve_rings:
            for ring in mol.GetRingInfo().AtomRings():
                if selected.intersection(ring):
                    selected.update(ring)
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a not in selected and b not in selected:
                continue
            charged = mol.GetAtomWithIdx(a).GetFormalCharge() or mol.GetAtomWithIdx(b).GetFormalCharge()
            terminal = mol.GetAtomWithIdx(a).GetDegree() == 1 or mol.GetAtomWithIdx(b).GetDegree() == 1
            if charged or terminal or bond.GetBondType() != Chem.BondType.SINGLE:
                selected.update((a, b))
        changed = selected != previous
    return capped_molecule(polymer.atoms, mol, selected, ending=ending)


def remap_snippet_topology(original_idxs, sn_atoms, sn_atypes, sn_bonds,
                           sn_angles, sn_dihedrals, sn_impropers, local_idxs):
    
    # Create filtered topology objects
    sn_bonds_out     = []
    sn_angles_out    = []
    sn_dihedrals_out = [] 
    sn_impropers_out = []
    
    # Prepare remapping IDs based on original IDs from snippet indices
    remap_ids = {}
    for cc, atype in enumerate(sn_atypes):
        remap_ids[atype.label] = str(original_idxs[cc])
    
    # Update bonds
    for bond in copy.deepcopy(sn_bonds):
        a1 = remap_ids[bond._a1]
        a2 = remap_ids[bond._a2]
        
        if all([ii in local_idxs for ii in [a1, a2]]):
            bond.update(a1=a1, a2=a2)
            sn_bonds_out.append(bond)
    
    # Update angles
    for angle in copy.deepcopy(sn_angles):
        a1 = remap_ids[angle._a1]
        a2 = remap_ids[angle._a2]
        a3 = remap_ids[angle._a3]
        if all([ii in local_idxs for ii in [a1, a2, a3]]):
            angle.update(a1=a1, a2=a2, a3=a3)
            sn_angles_out.append(angle)
    
    # Update dihedrals
    for dihedral in copy.deepcopy(sn_dihedrals):
        a1 = remap_ids[dihedral._a1]
        a2 = remap_ids[dihedral._a2]
        a3 = remap_ids[dihedral._a3]
        a4 = remap_ids[dihedral._a4]
        if all([ii in local_idxs for ii in [a1, a2, a3, a4]]):
            dihedral.update(a1=a1, a2=a2, a3=a3, a4=a4)
            sn_dihedrals_out.append(dihedral)
    
    # Update impropers
    for improper in copy.deepcopy(sn_impropers):
        a1 = remap_ids[improper._a1]
        a2 = remap_ids[improper._a2]
        a3 = remap_ids[improper._a3]
        a4 = remap_ids[improper._a4]
        if all([ii in local_idxs for ii in [a1, a2, a3, a4]]):
            # K, d, n = improper.values  # Assuming values consist of K, d, n parameters
            improper.update(a1=a1, a2=a2, a3=a3, a4=a4)
            sn_impropers_out.append(improper)
        
    return sn_bonds_out, sn_angles_out, sn_dihedrals_out, sn_impropers_out
