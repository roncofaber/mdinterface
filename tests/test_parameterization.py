import pickle

import networkx as nx
import numpy as np
import pytest
from ase import Atoms
from rdkit import Chem

from mdinterface import Specie, Polymer, SimCell
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral
from mdinterface.utils.graphs import find_unique_paths_of_length


def parameterize_mock(atoms, charge=None, is_snippet=False):
    bare = atoms.copy()
    bare.arrays.pop("stype", None)
    specie = Specie(bare)
    types = [Atom(a.symbol, label=f'{a.symbol}_{i}', eps=0.1, sig=2) for i, a in enumerate(atoms)]
    labels = [a.label for a in types]
    bonds = [Bond(labels[a], labels[b], kr=100, r0=1.5) for a, b in specie.graph.edges]
    angles = [Angle(*(labels[i] for i in ids), kr=20, theta0=109) for ids in find_unique_paths_of_length(specie.graph, 2)]
    dihedrals = [Dihedral(*(labels[i] for i in ids), 0, 0, 0, 0) for ids in find_unique_paths_of_length(specie.graph, 3)]
    result = atoms.copy()
    result.set_initial_charges(np.full(len(atoms), charge / len(atoms)))
    return result, types, bonds, angles, dihedrals, []


def snapshot(specie):
    return pickle.dumps({key: value for key, value in specie.__dict__.items() if key != '_graph'})


def test_parameterize_applies_parameters_without_moving_atoms(monkeypatch):
    monkeypatch.setattr('mdinterface.externals.ligpargen.run_ligpargen', parameterize_mock)
    specie = Specie(smiles='C[N+](C)(C)C')
    positions = specie.atoms.positions.copy()
    report = specie.parameterize()
    assert report['final_charge'] == pytest.approx(1)
    assert report['correction_per_atom'] == 0
    assert report['junctions'] == 0
    assert len(specie.bonds[0]) == specie.graph.number_of_edges()
    assert specie.validate_force_field() is None
    np.testing.assert_array_equal(specie.atoms.positions, positions)


def test_invalid_parameters_are_rejected_atomically(monkeypatch):
    specie = Specie(smiles='CC')
    before = snapshot(specie)
    def bad(*args, **kwargs):
        result = parameterize_mock(*args, **kwargs)
        result[1][0].eps = None
        return result
    monkeypatch.setattr('mdinterface.externals.ligpargen.run_ligpargen', bad)
    with pytest.raises(ValueError, match='pair parameters'):
        specie.parameterize()
    assert snapshot(specie) == before


def test_segment_limits_include_caps(monkeypatch):
    sizes = []
    def run(atoms, **kwargs):
        sizes.append(len(atoms))
        return parameterize_mock(atoms, **kwargs)
    monkeypatch.setattr('mdinterface.externals.ligpargen.run_ligpargen', run)
    specie = Specie(smiles='CCCCCCCCCCCCCCCCCCCC')
    report = specie.parameterize(segment_size=20, snippet_radius=4)
    segments = len(sizes) - report['junctions']
    assert segments == report['junctions'] + 1
    assert segments > 1
    assert max(sizes[:segments]) <= 20
    assert max(sizes) <= 200
    specie.validate_force_field()


def test_large_failure_preserves_original(monkeypatch):
    specie = Specie(smiles='CCCCCCCCCCCC')
    before = snapshot(specie)
    calls = []
    def run(atoms, **kwargs):
        calls.append(len(atoms))
        if len(calls) == 2:
            raise RuntimeError('backend failed')
        return parameterize_mock(atoms, **kwargs)
    monkeypatch.setattr('mdinterface.externals.ligpargen.run_ligpargen', run)
    with pytest.raises(RuntimeError, match='backend failed'):
        specie.parameterize(segment_size=18, snippet_radius=4)
    assert len(calls) == 2
    assert snapshot(specie) == before


def test_oversized_junction_fails_before_any_backend_call(monkeypatch):
    specie = Specie(smiles='CCCCCCCCCCCC')
    def unexpected(*args, **kwargs):
        raise AssertionError('backend should not run')
    monkeypatch.setattr('mdinterface.externals.ligpargen.run_ligpargen', unexpected)
    monkeypatch.setattr('mdinterface.build.snippets.make_snippet', lambda *args, **kwargs: (Atoms('H201'), np.arange(201)))
    with pytest.raises(ValueError, match='201 atoms'):
        specie.parameterize(segment_size=18)


def test_parameterized_copolymer_inherits_without_label_collisions(monkeypatch):
    monkeypatch.setattr('mdinterface.externals.ligpargen.run_ligpargen', parameterize_mock)
    first = Specie(smiles='[CH3:1][CH3:2]')
    second = Specie(smiles='[CH3:1]O[CH3:2]')
    for monomer in (first, second):
        monomer.parameterize()
        monomer.mark_attachment_sites(1, 2)
    before = [snapshot(first), snapshot(second)]
    chain = Polymer([first, second], sequence=[0, 1, 0])
    assert len(set(chain._sids)) == len(chain.atoms)
    assert len(chain.bonds[0]) == chain.graph.number_of_edges() - 2
    assert all(a.eps is not None for a in chain._stype)
    assert [snapshot(first), snapshot(second)] == before
    monkeypatch.setattr('mdinterface.core.polymer.run_ligpargen', parameterize_mock)
    chain.refine_junctions(snippet_radius=4)
    chain.validate_force_field()


def test_attachment_sites_follow_atom_maps_after_reordering():
    specie = Specie(smiles='[CH3:7]O[CH3:9]')
    reordered = Specie(specie.atoms[::-1])
    reordered.mark_attachment_sites(7, 9)
    mol = reordered.to_rdkit()
    for value, number in [(1, 7), (2, 9)]:
        index = np.flatnonzero(reordered.atoms.arrays['polymerize'] == value).item()
        atom = mol.GetAtomWithIdx(index)
        assert atom.GetSymbol() == 'H'
        assert atom.GetNeighbors()[0].GetAtomMapNum() == number
    before = reordered.atoms.arrays['polymerize'].copy()
    with pytest.raises(ValueError):
        reordered.mark_attachment_sites(7, 7)
    np.testing.assert_array_equal(reordered.atoms.arrays['polymerize'], before)


def test_explicit_zero_torsions_are_valid(monkeypatch):
    monkeypatch.setattr('mdinterface.externals.ligpargen.run_ligpargen', parameterize_mock)
    specie = Specie(smiles='CCCC')
    specie.parameterize()
    assert len(specie.dihedrals[0]) > 0
    assert all(tuple(d.values[:4]) == (0, 0, 0, 0) for d in specie._dtype)
    specie.validate_force_field()


def test_unparameterized_and_missing_bonds_rejected(monkeypatch):
    specie = Specie(smiles='CC')
    with pytest.raises(ValueError, match='pair parameters'):
        specie.validate_force_field()
    monkeypatch.setattr('mdinterface.externals.ligpargen.run_ligpargen', parameterize_mock)
    specie.parameterize()
    specie._old_bonds = []
    specie._update_topology_mappings()
    with pytest.raises(ValueError, match='unassigned bonds'):
        specie.validate_force_field()


def test_export_preflight_does_not_create_or_overwrite_file(tmp_path):
    specie = Specie(smiles='CC')
    box = SimCell([20, 20])
    box._all_species = [specie]
    box._universe = specie.to_universe()
    box._universe.dimensions = [20, 20, 20, 90, 90, 90]
    output = tmp_path / 'data.lammps'
    output.write_text('original')
    with pytest.raises(ValueError, match='Incomplete force field'):
        box.write_lammps(str(output))
    assert output.read_text() == 'original'


def test_export_expected_charge_is_explicit(tmp_path):
    from mdinterface.database import Ion
    specie = Ion('Na')
    box = SimCell([20, 20])
    box._all_species = [specie]
    box._universe = specie.to_universe()
    box._universe.dimensions = [20, 20, 20, 90, 90, 90]
    output = tmp_path / 'data.lammps'
    with pytest.raises(ValueError, match='expected_charge'):
        box.write_lammps(str(output), expected_charge=0)
    assert not output.exists()
    box.write_lammps(str(output), expected_charge=1)


def test_imported_connectivity_survives_coordinates_and_reordering():
    from pathlib import Path
    path = Path(__file__).parents[1] / 'examples/polymer/data/mon1/monomer_1.lammps.lmp'
    specie = Specie(lammps_data=path)
    expected = set(specie.graph.edges)
    specie.update_positions(positions=specie.atoms.positions * 10)
    assert set(specie.graph.edges) == expected
    reordered = Specie(specie.atoms[::-1])
    n = len(specie.atoms)
    assert {frozenset((n - 1 - a, n - 1 - b)) for a, b in expected} == {frozenset(edge) for edge in reordered.graph.edges}
    specie.atoms.set_cell([100, 100, 100])
    specie.repeat([2, 1, 1])
    assert nx.number_connected_components(specie.graph) == 2
    assert specie.graph.number_of_edges() == 2 * len(expected)


def test_refinement_retains_distinct_improper_assignments():
    from mdinterface.core.topology import Improper
    specie = Specie(smiles='CC')
    first = Improper('C_0', 'H_2', 'H_3', 'C_1', K=1, d=-1, n=2)
    second = Improper('C_1', 'H_5', 'H_6', 'C_0', K=1, d=-1, n=2)
    assert first == second
    assert not first.__eq_strict__(second)
    specie._add_to_topology(impropers=[first])
    specie._add_to_topology(impropers=[second, first.copy()])
    assert len(specie._old_impropers) == 2
    assert set(specie._imap) == {first.symbols, second.symbols}


@pytest.mark.parametrize('name, charge', [('Na', 1), ('Cl', -1), ('Hydronium', 1), ('Hydroxide', -1), ('Perchlorate', -1)])
def test_database_ions_have_their_formal_total(name, charge):
    from mdinterface import database
    specie = database.Ion(name, chg_scaling=0.8) if name in ('Na', 'Cl') else getattr(database, name)()
    assert specie._resolve_charge(None) == charge
    if name in ('Na', 'Cl'):
        assert specie.charges.sum() == pytest.approx(0.8 * charge)


def test_export_does_not_clobber_fixed_temporary_filename(tmp_path, monkeypatch):
    from mdinterface.database import Ion
    monkeypatch.chdir(tmp_path)
    original = tmp_path / 'tmp_data.lammps'
    original.write_text('user file')
    specie = Ion('Na')
    box = SimCell([20, 20])
    box._all_species = [specie]
    box._universe = specie.to_universe()
    box._universe.dimensions = [20, 20, 20, 90, 90, 90]
    box.write_lammps('export.data')
    assert original.read_text() == 'user file'
    assert not list(tmp_path.glob('.mdinterface-*.lammps'))


@pytest.mark.parametrize('residual', [0.0, 5e-13, -5e-13, 5e-9])
@pytest.mark.parametrize('policy', ['none', 'uniform'])
def test_charge_audit_normalizes_roundoff_without_erasing_partial_charges(monkeypatch, residual, policy):
    specie = Specie(smiles='C')
    charges = np.array([0.1, 0.2, -0.3, 0, residual])
    specie.atoms.set_initial_charges(charges)

    def parameterize(atoms, **kwargs):
        result = parameterize_mock(atoms, **kwargs)
        result[0].set_initial_charges(charges)
        return result

    monkeypatch.setattr('mdinterface.externals.ligpargen.run_ligpargen', parameterize)
    report = specie.parameterize(charge_correction=policy)
    if abs(residual) <= 1e-12:
        for key in ('initial_charge', 'refined_charge', 'residual', 'final_charge', 'correction_per_atom'):
            assert report[key] == 0.0
            assert not np.signbit(report[key])
        np.testing.assert_array_equal(specie.charges, charges)
    else:
        assert report['residual'] == pytest.approx(residual, abs=1e-16)
        if policy == 'none':
            assert report['final_charge'] != 0
            np.testing.assert_array_equal(specie.charges, charges)
        else:
            assert report['final_charge'] == 0
            assert report['correction_per_atom'] != 0
