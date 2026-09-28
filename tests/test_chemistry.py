import numpy as np
import pytest
import networkx as nx
from rdkit import Chem
from ase import Atoms

from mdinterface import Specie, Polymer
from mdinterface.core.chemistry import stored_molecule, capped_molecule
from mdinterface.build.snippets import make_snippet


def test_smiles_preserves_chemistry_and_coordinates():
    specie = Specie(smiles='C[N+](C)(C)CCO', seed=19)
    assert len(specie.atoms) == 21
    mol = specie.to_rdkit()
    assert Chem.GetFormalCharge(mol) == 1
    nitrogen = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'N')
    assert specie.atoms.arrays['nominal_charge'][nitrogen] == 1
    assert specie._tot_charge == 1
    np.testing.assert_allclose(mol.GetConformer().GetPositions(), specie.atoms.positions)
    assert len(set(specie.atoms.arrays['atom_map'])) == len(specie.atoms)
    assert specie.to_smiles() == 'C[N+](C)(C)CCO'
    mol.GetAtomWithIdx(nitrogen).SetFormalCharge(0)
    assert Chem.GetFormalCharge(specie.to_rdkit()) == 1


def test_aromatic_stereo_and_seed_roundtrip():
    smiles = 'C[C@H](O)c1ccccc1'
    first = Specie(smiles=smiles, seed=17)
    second = Specie(smiles=smiles, seed=17)
    assert first.to_smiles() == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    np.testing.assert_allclose(first.atoms.positions, second.atoms.positions)
    assert any(data['bond_order'] == 1.5 for _, _, data in first.graph.edges(data=True))
    before = set(first.graph.edges)
    first.update_positions(positions=first.atoms.positions * 10)
    assert set(first.graph.edges) == before


def test_rdkit_input_not_mutated():
    mol = Chem.MolFromSmiles('[CH3:7][OH:9]')
    specie = Specie(mol)
    assert mol.GetNumAtoms() == 2
    assert mol.GetNumConformers() == 0
    assert specie.atoms.arrays['atom_map'][:2].tolist() == [7, 9]
    assert ':7' in specie.to_smiles(mapped=True)


@pytest.mark.parametrize('kwargs', [
    {'smiles': 'not smiles'},
    {'smiles': 'O', 'atoms': Atoms('H')},
    {'smiles': '[Na+]', 'tot_charge': 0},
    {'smiles': '[Na+].[Cl-]'},
])
def test_invalid_chemical_inputs_raise(kwargs):
    with pytest.raises(ValueError):
        Specie(**kwargs)


def test_ase_names_keep_their_meaning():
    assert len(Specie('CO').atoms) == 2
    assert len(Specie(smiles='CO').atoms) == 6


def test_formal_charge_cannot_be_overwritten_silently():
    specie = Specie(smiles='C[N+](C)(C)C')
    specie.atoms.arrays['nominal_charge'][:] = 0
    with pytest.raises(ValueError, match='conflicts'):
        specie.to_rdkit()


def ethane_monomer():
    monomer = Specie(smiles='CC', charges=0.0)
    mol = monomer.to_rdkit()
    marks = np.zeros(len(monomer.atoms), dtype=int)
    for carbon, mark in [(0, 1), (1, 2)]:
        hydrogen = next(a.GetIdx() for a in mol.GetAtomWithIdx(carbon).GetNeighbors() if a.GetAtomicNum() == 1)
        marks[hydrogen] = mark
    monomer.atoms.set_array('polymerize', marks)
    return monomer


def test_polymer_junctions_are_explicit_and_survive_motion():
    monomer = ethane_monomer()
    chain = Polymer(monomer, nrep=3)
    assert chain.to_smiles() == 'CCCCCC'
    assert nx.is_connected(chain.graph)
    assert len(chain.junction_bonds) == 2
    for a, b in chain.junction_bonds:
        assert chain.atoms.get_distance(a, b) == pytest.approx(1.52)
    maps = chain.atoms.arrays['atom_map']
    assert len(set(maps)) == len(maps)
    assert maps[0] == monomer.atoms.arrays['atom_map'][0]
    pairs = chain.junction_bonds
    chain.update_positions(positions=chain.atoms.positions * 5)
    assert chain.junction_bonds == pairs
    assert nx.is_connected(chain.graph)


def test_snippet_preserves_charged_group_and_neutral_caps():
    specie = Specie(smiles='CCCC[N+](C)(C)CCCC')
    mol = specie.to_rdkit()
    nitrogen = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
    snippet, indices = make_snippet(specie, nitrogen, 1)
    capped = stored_molecule(snippet)
    assert Chem.GetFormalCharge(capped) == 1
    assert sum(a.GetAtomicNum() == 7 for a in capped.GetAtoms()) == 1
    for index, original in enumerate(indices):
        if snippet[index].number != specie.atoms[original].number:
            assert snippet.arrays['nominal_charge'][index] == 0
    Chem.SanitizeMol(capped)


def test_capping_rejects_charged_boundary():
    specie = Specie(smiles='CC[N+](C)(C)C')
    with pytest.raises(ValueError, match='formal-charge site'):
        capped_molecule(specie.atoms, specie.to_rdkit(), [0, 1])


def test_refinement_does_not_reuse_stale_charge_cache(monkeypatch):
    import mdinterface.core.polymer as module
    from mdinterface.core.topology import Atom
    chain = Polymer(ethane_monomer(), nrep=2)
    calls = []

    def parameterize(atoms, charge, is_snippet):
        calls.append(charge)
        result = atoms.copy()
        result.set_initial_charges(np.full(len(atoms), 0.02))
        types = [Atom(a.symbol, label=f'{a.symbol}_{i}', eps=0.1, sig=1) for i, a in enumerate(atoms)]
        return result, types, [], [], [], []

    monkeypatch.setattr(module, 'run_ligpargen', parameterize)
    report = chain.refine_junctions(charge_correction="none")
    assert chain.charges.sum() > 0
    assert report['initial_charge'] == 0
    assert report['refined_charge'] == pytest.approx(chain.charges.sum())
    assert report['residual'] == pytest.approx(chain.charges.sum())
    assert report['correction_per_atom'] == 0
    assert report['junctions'] == 1
    report = chain.refine_junctions(charge_correction="uniform")
    assert report['final_charge'] == 0.0
    assert report['correction_per_atom'] * len(chain.atoms) == pytest.approx(-report['residual'])
    assert calls == [0, 0]
    assert chain.charges.sum() == pytest.approx(0, abs=1e-12)


def test_ligpargen_charges_are_applied(monkeypatch):
    import mdinterface.externals.ligpargen as module
    from mdinterface.core.topology import Atom

    def parameterize(atoms, charge):
        result = atoms.copy()
        result.set_initial_charges([1.0])
        return result, [Atom('Na', label='Na', eps=0.1, sig=2)], [], [], [], []

    monkeypatch.setattr(module, 'run_ligpargen', parameterize)
    specie = Specie(smiles='[Na+]', ligpargen=True)
    assert specie.charges.tolist() == [1.0]


def test_charged_ase_monomer_wrong_formal_site_rejected():
    from mdinterface.core.chemistry import perceive_molecule
    specie = Specie(smiles='C[N+](C)(C)C')
    atoms = specie.atoms.copy()
    atoms.info.clear()
    atoms.arrays['nominal_charge'][:] = 0
    atoms.arrays['nominal_charge'][0] = 1
    with pytest.raises(ValueError, match='sites disagree'):
        perceive_molecule(atoms)


def test_segment_keeps_chemical_graph_and_neutral_caps():
    from mdinterface.externals.ligpargen import _make_capped_segment
    specie = Specie(smiles='CCCC')
    mol = specie.to_rdkit()
    selected = {0, 1}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and any(n.GetIdx() in {0, 1} for n in atom.GetNeighbors()):
            selected.add(atom.GetIdx())
    capped, n_real = _make_capped_segment(specie, sorted(selected), [(1, 2)])
    result = stored_molecule(capped)
    assert n_real == len(selected)
    assert Chem.MolToSmiles(Chem.RemoveHs(result)).count('C') == 2
    assert capped.arrays['nominal_charge'][-1] == 0


@pytest.mark.integration
def test_smiles_packing_keeps_molecule_counts():
    from mdinterface import SimCell
    specie = Specie(smiles='CCO', name='ETOH')
    box = SimCell(xysize=[20, 20], verbose=False)
    box.add_solvent(specie, nsolvent=4, zdim=20)
    box.build()
    assert len(box.universe.atoms) == 4 * len(specie.atoms)
    assert len(box.universe.residues) == 4


def test_atom_maps_preserve_chemistry_after_ase_reordering():
    specie = Specie(smiles='CC(=O)[O-]', name='ACET')
    order = list(reversed(range(len(specie.atoms))))
    reordered = Specie(specie.atoms[order])
    assert reordered.to_smiles() == specie.to_smiles()
    mol = reordered.to_rdkit()
    np.testing.assert_allclose(mol.GetConformer().GetPositions(), reordered.atoms.positions)
    assert [a.GetAtomMapNum() for a in mol.GetAtoms()] == specie.atoms.arrays['atom_map'][order].tolist()


def test_atom_map_corruption_rejected():
    specie = Specie(smiles='CCO')
    specie.atoms.arrays['atom_map'][1] = specie.atoms.arrays['atom_map'][0]
    with pytest.raises(ValueError, match='atom maps'):
        specie.to_rdkit()


def test_polymer_conformer_preserves_chemistry_and_arrays():
    first = Polymer(ethane_monomer(), nrep=3)
    second = Polymer(ethane_monomer(), nrep=3)
    first.atoms.positions *= 10
    arrays = {key: value.copy() for key, value in first.atoms.arrays.items() if key != 'positions'}
    smiles = first.to_smiles(mapped=True)
    pairs = first.junction_bonds
    energy = first.generate_conformer(minimize=True, seed=42)
    assert np.isfinite(energy)
    assert second.generate_conformer(minimize=True, seed=42) == pytest.approx(energy)
    np.testing.assert_allclose(first.atoms.positions, second.atoms.positions)
    for key, value in arrays.items():
        np.testing.assert_array_equal(first.atoms.arrays[key], value)
    assert first.to_smiles(mapped=True) == smiles
    assert first.junction_bonds == pairs
    for a, b in pairs:
        assert 1.4 < first.atoms.get_distance(a, b) < 1.7


@pytest.mark.parametrize('failure', ['embedding', 'parameters', 'convergence'])
def test_polymer_conformer_failure_preserves_coordinates(monkeypatch, failure):
    from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
    chain = Polymer(ethane_monomer(), nrep=3)
    positions = chain.atoms.positions.copy()
    if failure == 'embedding':
        monkeypatch.setattr(rdDistGeom, 'EmbedMolecule', lambda *args: -1)
    elif failure == 'parameters':
        monkeypatch.setattr(rdForceFieldHelpers, 'MMFFHasAllMoleculeParams', lambda mol: False)
    else:
        class Unconverged:
            def Minimize(self, **kwargs):
                return 1
        monkeypatch.setattr(rdForceFieldHelpers, 'MMFFGetMoleculeForceField', lambda *args: Unconverged())
    with pytest.raises((ValueError, RuntimeError)):
        chain.generate_conformer(minimize=True)
    np.testing.assert_array_equal(chain.atoms.positions, positions)


@pytest.mark.parametrize('kwargs', [{'seed': -1}, {'max_iterations': 0}])
def test_polymer_conformer_invalid_arguments(kwargs):
    chain = Polymer(ethane_monomer(), nrep=2)
    with pytest.raises(ValueError):
        chain.generate_conformer(**kwargs)


def test_charged_polymer_conformer_preserves_charge_sites():
    monomer = Specie(smiles='C[N+](C)(C)CC')
    mol = monomer.to_rdkit()
    marks = np.zeros(len(monomer.atoms), dtype=int)
    for carbon, mark in [(0, 1), (5, 2)]:
        hydrogen = next(a.GetIdx() for a in mol.GetAtomWithIdx(carbon).GetNeighbors() if a.GetAtomicNum() == 1)
        marks[hydrogen] = mark
    monomer.atoms.set_array('polymerize', marks)
    chain = Polymer(monomer, nrep=2)
    chain.atoms.set_initial_charges(np.linspace(-0.1, 0.2, len(chain.atoms)))
    formal = chain.atoms.arrays['nominal_charge'].copy()
    partial = chain.charges.copy()
    chain.generate_conformer(minimize=True, seed=42, max_iterations=2000)
    np.testing.assert_array_equal(chain.atoms.arrays['nominal_charge'], formal)
    np.testing.assert_array_equal(chain.charges, partial)
    assert Chem.GetFormalCharge(chain.to_rdkit()) == 2


@pytest.mark.parametrize('kind', ['periodic', 'constrained'])
def test_polymer_conformer_rejects_unsupported_boundary_conditions(kind):
    from ase.constraints import FixAtoms
    chain = Polymer(ethane_monomer(), nrep=2)
    if kind == 'periodic':
        chain.atoms.set_pbc(True)
    else:
        chain.atoms.set_constraint(FixAtoms(indices=[0]))
    with pytest.raises(ValueError, match='unconstrained, nonperiodic'):
        chain.generate_conformer(minimize=True)


def test_specie_minimization_uses_existing_conformer(monkeypatch):
    from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
    specie = Specie(smiles='CCCC', seed=42)
    mol = specie.to_rdkit()
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, rdForceFieldHelpers.MMFFGetMoleculeProperties(mol))
    initial = ff.CalcEnergy()
    def unexpected_embed(*args):
        raise AssertionError('Minimization must not embed another conformer')
    monkeypatch.setattr(rdDistGeom, 'EmbedMolecule', unexpected_embed)
    assert specie.minimize_geometry() <= initial
    assert specie.to_smiles() == 'CCCC'


def test_conformer_can_embed_without_mmff(monkeypatch):
    from rdkit.Chem import rdForceFieldHelpers
    specie = Specie(smiles='CCO')
    monkeypatch.setattr(rdForceFieldHelpers, 'MMFFHasAllMoleculeParams', lambda mol: False)
    assert specie.generate_conformer(seed=42) is None


def test_junction_failure_does_not_modify_original(monkeypatch):
    import pickle
    import mdinterface.core.polymer as module
    from mdinterface.core.topology import Atom
    chain = Polymer(ethane_monomer(), nrep=3)
    before = pickle.dumps({key: value for key, value in chain.__dict__.items() if key != "_graph"})
    graph = chain.graph.copy()
    atoms = chain.atoms
    calls = []
    def parameterize(atoms, charge, is_snippet):
        calls.append(charge)
        if len(calls) == 2:
            raise RuntimeError('second junction failed')
        result = atoms.copy()
        result.set_initial_charges(np.full(len(atoms), 0.02))
        types = [Atom(a.symbol, label=f'{a.symbol}_{i}', eps=0.123, sig=2.5) for i, a in enumerate(atoms)]
        return result, types, [], [], [], []
    monkeypatch.setattr(module, 'run_ligpargen', parameterize)
    with pytest.raises(RuntimeError, match='second junction'):
        chain.refine_junctions(charge_correction='uniform')
    assert len(calls) == 2
    assert chain.atoms is atoms
    assert pickle.dumps({key: value for key, value in chain.__dict__.items() if key != "_graph"}) == before
    assert nx.utils.graphs_equal(chain.graph, graph)


@pytest.mark.parametrize('kwargs', [
    {'snippet_radius': 3}, {'snippet_radius': 4.5},
    {'charge_correction': True}, {'cap_element': 'C'},
])
def test_junction_options_fail_before_parameterization(monkeypatch, kwargs):
    import mdinterface.core.polymer as module
    chain = Polymer(ethane_monomer(), nrep=2)
    def unexpected(*args, **kwargs):
        raise AssertionError('Invalid options must fail before backend execution')
    monkeypatch.setattr(module, 'run_ligpargen', unexpected)
    with pytest.raises(ValueError):
        chain.refine_junctions(**kwargs)


@pytest.mark.parametrize('method', ['ligpargen', 'obabel', 'resp', 'opls'])
def test_charge_helpers_use_stored_ionic_charge(monkeypatch, method):
    import mdinterface.core.specie as module
    specie = Specie(smiles='C[N+](C)(C)C')
    seen = []
    expected = np.full(len(specie.atoms), 1 / len(specie.atoms))
    def ligpargen(atoms, charge):
        seen.append(charge)
        result = atoms.copy()
        result.set_initial_charges(expected)
        return result, [], [], [], [], []
    def obabel(atoms, charge):
        seen.append(charge)
        return expected
    def resp(specie, charge, **kwargs):
        seen.append(charge)
        return expected, specie.atoms.copy()
    monkeypatch.setattr(module, 'run_ligpargen', ligpargen)
    monkeypatch.setattr(module, 'run_OBChargeModel', obabel)
    monkeypatch.setattr(module, 'calculate_RESP_charges', resp)
    call = specie.estimate_OPLSAA_parameters if method == 'opls' else lambda **kwargs: specie.estimate_charges(method=method, **kwargs)
    call()
    assert seen == [1]
    with pytest.raises(ValueError, match='conflicts'):
        call(charge=0)
    assert seen == [1]


def test_resp_receives_charge_when_constructing_molecule(monkeypatch):
    import sys
    from types import SimpleNamespace
    from mdinterface.externals.pyscf import calculate_RESP_charges
    specie = Specie(smiles='C[N+](C)(C)C')
    seen = []
    def convert(atoms, basis, charge):
        seen.append(charge)
        raise RuntimeError('conversion reached')
    monkeypatch.setitem(sys.modules, 'gpu4pyscf.pop', SimpleNamespace(esp=None))
    monkeypatch.setitem(sys.modules, 'pymbxas.build.structure', SimpleNamespace(ase_to_mole=convert, mole_to_ase=None))
    monkeypatch.setitem(sys.modules, 'pymbxas.build.input_pyscf', SimpleNamespace(make_pyscf_calculator=None))
    monkeypatch.setitem(sys.modules, 'pymbxas.md.solvers', SimpleNamespace(Geometry_optimizer=None))
    with pytest.raises(RuntimeError, match='conversion reached'):
        calculate_RESP_charges(specie)
    assert seen == [1]


def test_openbabel_preserves_ionic_total():
    pytest.importorskip('openbabel.openbabel')
    specie = Specie(smiles='C[N+](C)(C)C')
    assert specie.estimate_charges(method='obabel').sum() == pytest.approx(1.0)


def test_explicit_hydrogen_maps_survive_smiles_parsing():
    specie = Specie(smiles='[H:7][C:1]([H:8])([H:9])[C:2]([H:10])([H:11])[H:12]')
    maps = {a.GetAtomMapNum(): a.GetSymbol() for a in specie.to_rdkit().GetAtoms()}
    assert maps == {1: 'C', 2: 'C', 7: 'H', 8: 'H', 9: 'H', 10: 'H', 11: 'H', 12: 'H'}
    specie.mark_attachment_sites(head_map=7, tail_map=12)
    marks = specie.atoms.arrays['polymerize']
    assert specie.atoms.arrays['atom_map'][marks == 1].tolist() == [7]
    assert specie.atoms.arrays['atom_map'][marks == 2].tolist() == [12]


def test_repeat_scales_ionic_charge_and_keeps_original_atom_maps():
    specie = Specie(smiles='[CH3:7][N+:9](C)(C)C')
    original_maps = specie.atoms.arrays['atom_map'].copy()
    specie.atoms.set_cell([20, 20, 20])
    specie.repeat([2, 1, 1])
    mol = specie.to_rdkit()
    assert Chem.GetFormalCharge(mol) == 2
    assert specie._resolve_charge(None) == 2
    np.testing.assert_array_equal(specie.atoms.arrays['atom_map'][:len(original_maps)], original_maps)
    assert len(set(specie.atoms.arrays['atom_map'])) == len(specie.atoms)
    assert len(Chem.GetMolFrags(mol)) == 2


def test_polymer_charge_audit_normalizes_roundoff(monkeypatch):
    chain = Polymer(ethane_monomer(), nrep=2)
    charges = np.zeros(len(chain.atoms))
    charges[:3] = [0.1, 0.2, -0.3]
    chain.atoms.set_initial_charges(charges)
    monkeypatch.setattr(Polymer, '_update_connection', lambda *args, **kwargs: None)
    report = chain.refine_junctions(charge_correction='uniform')
    for key in ('initial_charge', 'refined_charge', 'residual', 'final_charge', 'correction_per_atom'):
        assert report[key] == 0.0
        assert not np.signbit(report[key])
    np.testing.assert_array_equal(chain.charges, charges)
