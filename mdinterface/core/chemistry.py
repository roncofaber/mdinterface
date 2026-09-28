"""RDKit chemical graphs with atom ordering shared by ASE coordinates."""

import numpy as np
import networkx as nx
from ase import Atoms
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdDistGeom


_MOLECULE_KEY = "mdinterface_molecule"


def _with_coordinates(mol, atoms):
    mol = Chem.Mol(mol)
    if [atom.GetAtomicNum() for atom in mol.GetAtoms()] != atoms.numbers.tolist():
        raise ValueError("Chemical graph and ASE atom ordering differ.")
    conformer = Chem.Conformer(len(atoms))
    for index, position in enumerate(atoms.positions):
        conformer.SetAtomPosition(index, position)
    mol.RemoveAllConformers()
    mol.AddConformer(conformer)
    return mol


def store_molecule(atoms, mol):
    mol = _with_coordinates(mol, atoms)
    used = set()
    next_map = max((atom.GetAtomMapNum() for atom in mol.GetAtoms()), default=0) + 1
    for atom in mol.GetAtoms():
        number = atom.GetAtomMapNum()
        if number and number in used:
            raise ValueError("Atom-map numbers must be unique within a molecule.")
        if not number:
            number = next_map
            next_map += 1
            atom.SetAtomMapNum(number)
        used.add(number)
    atoms.info[_MOLECULE_KEY] = mol.ToBinary()
    atoms.set_array("nominal_charge", np.array([a.GetFormalCharge() for a in mol.GetAtoms()], dtype=int))
    atoms.set_array("atom_map", np.array([a.GetAtomMapNum() for a in mol.GetAtoms()], dtype=int))


def stored_molecule(atoms):
    if _MOLECULE_KEY not in atoms.info:
        return None
    mol = Chem.Mol(atoms.info[_MOLECULE_KEY])
    if "atom_map" in atoms.arrays:
        original = [a.GetAtomMapNum() for a in mol.GetAtoms()]
        current = atoms.arrays["atom_map"].tolist()
        if len(current) != len(original) or len(set(current)) != len(current) or set(current) != set(original):
            raise ValueError("ASE atom maps no longer match the stored chemical graph.")
        if current != original:
            order = {number: index for index, number in enumerate(original)}
            mol = Chem.RenumberAtoms(mol, [order[number] for number in current])
    mol = _with_coordinates(mol, atoms)
    formal = np.array([a.GetFormalCharge() for a in mol.GetAtoms()])
    if "nominal_charge" in atoms.arrays and not np.array_equal(formal, atoms.arrays["nominal_charge"]):
        raise ValueError("nominal_charge conflicts with the stored RDKit formal charges.")
    return mol


def perceive_molecule(atoms, charge=None):
    mol = stored_molecule(atoms)
    if mol is not None:
        if charge is not None and Chem.GetFormalCharge(mol) != charge:
            raise ValueError("Total charge conflicts with the chemical graph.")
        return mol
    if any(atoms.pbc):
        raise ValueError("RDKit bond perception requires a nonperiodic molecular structure.")
    if charge is None:
        charge = int(np.sum(atoms.arrays.get("nominal_charge", [0])))
    rw = Chem.RWMol()
    for number in atoms.numbers:
        rw.AddAtom(Chem.Atom(int(number)))
    mol = _with_coordinates(rw.GetMol(), atoms)
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=int(charge))
        Chem.SanitizeMol(mol)
    except (ValueError, RuntimeError) as exc:
        raise ValueError("Cannot determine a valid molecular graph; provide charged SMILES or an RDKit molecule with explicit bonds.") from exc
    if "nominal_charge" in atoms.arrays and np.any(atoms.arrays["nominal_charge"]):
        formal = [a.GetFormalCharge() for a in mol.GetAtoms()]
        if not np.array_equal(formal, atoms.arrays["nominal_charge"]):
            raise ValueError("nominal_charge sites disagree with the perceived chemical graph; provide an explicit RDKit molecule.")
    return mol


def atoms_from_molecule(mol, seed=0):
    mol = Chem.Mol(mol)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addCoords=True)
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Use a separate Specie for each disconnected SMILES/RDKit fragment.")
    if not mol.GetNumConformers():
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = seed
        if rdDistGeom.EmbedMolecule(mol, params) != 0:
            raise ValueError("RDKit could not generate coordinates for this molecule.")
    atoms = Atoms([a.GetAtomicNum() for a in mol.GetAtoms()], positions=mol.GetConformer().GetPositions())
    store_molecule(atoms, mol)
    return atoms


def graph_from_molecule(mol):
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        graph.add_node(atom.GetIdx(), element=atom.GetSymbol(), formal_charge=atom.GetFormalCharge(), atom_map=atom.GetAtomMapNum())
    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_order=bond.GetBondTypeAsDouble())
    return graph


def capped_molecule(atoms, mol, selected, ending="H"):
    from ase.data import covalent_radii, atomic_numbers

    if ending not in {"H", "F", "Cl", "Br", "I"}:
        raise ValueError("Snippet caps must be neutral monovalent atoms.")
    selected = sorted(set(selected))
    retained = set(selected)
    cuts = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (a in retained) == (b in retained):
            continue
        if bond.GetBondType() != Chem.BondType.SINGLE:
            raise ValueError("Snippet boundaries may only cut single bonds.")
        inside, outside = (a, b) if a in retained else (b, a)
        if mol.GetAtomWithIdx(inside).GetFormalCharge() or mol.GetAtomWithIdx(outside).GetFormalCharge():
            raise ValueError("Snippet boundary cuts a formal-charge site; increase snippet_radius.")
        cuts.append((inside, outside))
    rw = Chem.RWMol(mol)
    for index in reversed(range(len(atoms))):
        if index not in retained:
            rw.RemoveAtom(index)
    result = atoms[selected]
    mapping = list(selected)
    for inside, outside in cuts:
        cap = Chem.Atom(ending)
        cap.SetNoImplicit(True)
        new_index = rw.AddAtom(cap)
        rw.AddBond(selected.index(inside), new_index, Chem.BondType.SINGLE)
        direction = atoms.positions[outside] - atoms.positions[inside]
        distance = np.linalg.norm(direction)
        if distance == 0:
            raise ValueError("Cannot cap a zero-length bond.")
        length = covalent_radii[atoms.numbers[inside]] + covalent_radii[atomic_numbers[ending]]
        result += Atoms(ending, positions=[atoms.positions[inside] + direction * length / distance])
        mapping.append(outside)
    capped = rw.GetMol()
    Chem.SanitizeMol(capped)
    store_molecule(result, capped)
    return result, np.array(mapping, dtype=int)


def minimize_molecule(mol, max_iterations):
    from rdkit.Chem import rdForceFieldHelpers

    if not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer.")
    if not rdForceFieldHelpers.MMFFHasAllMoleculeParams(mol):
        raise ValueError("MMFF94 parameters are unavailable for this molecule's chemistry.")
    properties = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)
    forcefield = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, properties)
    if forcefield is None or forcefield.Minimize(maxIts=int(max_iterations)) != 0:
        raise RuntimeError("MMFF94 minimization did not converge; increase max_iterations or use another conformer.")
    energy = float(forcefield.CalcEnergy())
    if not np.isfinite(energy) or not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise RuntimeError("RDKit produced nonfinite coordinates or energy.")
    return energy
