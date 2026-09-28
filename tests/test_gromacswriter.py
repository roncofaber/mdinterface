import json
import os
from pathlib import Path
import shutil
import subprocess

import MDAnalysis as mda
import networkx as nx
import numpy as np
import pytest

from mdinterface import SimCell, Specie
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral, Improper
from mdinterface.database import Water, Ion
from mdinterface.io.gromacswriter import write_gromacs_itp
from mdinterface.utils.graphs import find_unique_paths_of_length, find_improper_idxs


def parameterized(smiles='CCCC', coefficients=(1, 2, 3, 4)):
    bare = Specie(smiles=smiles, name='TEST')
    labels = [f'{atom.symbol}_{i}' for i, atom in enumerate(bare.atoms)]
    types = [Atom(atom.symbol, label, eps=0.1, sig=2.5)
             for atom, label in zip(bare.atoms, labels)]
    bonds = [Bond(labels[a], labels[b], kr=100, r0=1.5) for a, b in bare.graph.edges]
    angles = [Angle(*(labels[i] for i in ids), kr=20, theta0=109)
              for ids in find_unique_paths_of_length(bare.graph, 2)]
    torsions = [Dihedral(*(labels[i] for i in ids), *coefficients)
                for ids in find_unique_paths_of_length(bare.graph, 3)]
    return Specie(bare.atoms, name='TEST', atom_types=types, bonds=bonds,
                  angles=angles, dihedrals=torsions)


def assembled(species):
    cell = SimCell(xysize=[50, 50])
    cell._all_species = species
    cell._update_topology_indexes()
    groups = []
    for i, specie in enumerate(species):
        universe = specie.to_universe()
        universe.atoms.positions += [5 + 6 * i, 5, 5]
        groups.append(universe.atoms)
    cell._universe = mda.Merge(*groups)
    cell._universe.dimensions = [50, 50, 50, 90, 90, 90]
    cell._universe.add_TopologyAttr('resids', np.arange(1, len(species) + 1))
    return cell


def sections(path):
    result = {}
    section = None
    for line in Path(path).read_text().splitlines():
        line = line.split(';')[0].strip()
        if line.startswith('['):
            section = line.strip('[] ')
            result.setdefault(section, [])
        elif line and section:
            result[section].append(line.split())
    return result


@pytest.mark.parametrize('smiles', ['CCCC', 'C1CCC1', 'C1CCCCC1'])
def test_pairs_use_shortest_bond_distance_even_for_zero_torsions(tmp_path, smiles):
    specie = parameterized(smiles, (0, 0, 0, 0))
    write_gromacs_itp(specie, tmp_path / 'm.itp')
    rows = sections(tmp_path / 'm.itp').get('pairs', [])
    pairs = [(int(row[0]) - 1, int(row[1]) - 1) for row in rows]
    expected = {(a, b) for a in specie.graph for b in specie.graph
                if a < b and nx.shortest_path_length(specie.graph, a, b) == 3}
    assert set(pairs) == expected
    assert len(pairs) == len(expected)
    if smiles == 'C1CCC1':
        assert (0, 3) not in expected


def test_mixed_export_has_global_types_and_coordinate_order(tmp_path):
    water, ion = Water(), Ion('Na')
    cell = assembled([water, ion, water.copy()])
    cell.write_gromacs(outdir=tmp_path)
    top = sections(tmp_path / 'system.top')
    assert top['molecules'] == [['H2O', '1'], ['Na', '1'], ['H2O', '1']]
    text = (tmp_path / 'system.top').read_text()
    assert text.index('[ atomtypes ]') < text.index('#include')
    for name in ('H2O', 'Na'):
        assert 'atomtypes' not in sections(tmp_path / f'{name}.itp')
    coords = mda.Universe(str(tmp_path / 'system.gro'))
    assert list(coords.residues.resnames) == ['H2O', 'Na', 'H2O']
    np.testing.assert_allclose(coords.atoms.positions, cell.universe.atoms.positions, atol=0.006)


@pytest.mark.parametrize('kind', ['bare', 'nonfinite', 'missing_bond', 'fifth_term', 'bad_improper'])
def test_invalid_parameters_preserve_existing_output(tmp_path, kind):
    specie = parameterized()
    if kind == 'bare':
        specie = Specie(smiles='CC')
    elif kind == 'nonfinite':
        specie.atoms.set_initial_charges(np.full(len(specie.atoms), np.nan))
    elif kind == 'missing_bond':
        specie._btype[0].kr = None
    elif kind == 'fifth_term':
        specie._dtype[0].update(A5=1)
    else:
        specie = parameterized('C=C')
        ids = find_improper_idxs(specie.graph)[0]
        labels = tuple(specie._sids[ids])
        imp = Improper(*labels, K=1, d=1, n=2)
        imp._d = 0
        specie._itype = [imp]
        specie._imap = {labels: 0}
    path = tmp_path / 'm.itp'
    path.write_text('keep')
    with pytest.raises(ValueError):
        write_gromacs_itp(specie, path)
    assert path.read_text() == 'keep'


def test_conflicting_names_fail_before_writing(tmp_path):
    cell = assembled([Water(), Water(model='charmm')])
    outdir = tmp_path / 'output'
    with pytest.raises(ValueError, match='Conflicting GROMACS molecule'):
        cell.write_gromacs(outdir=outdir)
    assert not outdir.exists()


def test_changed_universe_charges_are_not_silently_discarded(tmp_path):
    cell = assembled([Water()])
    cell.universe.atoms.charges += 0.1
    with pytest.raises(ValueError, match='does not match'):
        cell.write_gromacs(outdir=tmp_path / 'output')
    assert not (tmp_path / 'output').exists()


def test_bond_angle_and_lj_units(tmp_path):
    write_gromacs_itp(Water(), tmp_path / 'water.itp')
    data = sections(tmp_path / 'water.itp')
    assert list(map(float, data['bonds'][0][3:])) == pytest.approx([0.09572, 376560])
    assert list(map(float, data['angles'][0][4:])) == pytest.approx([104.52, 460.24])
    assert list(map(float, data['atomtypes'][0][-2:])) == pytest.approx([0.3188, 0.426768])


def test_opls_energy_scan(tmp_path):
    specie = parameterized()
    write_gromacs_itp(specie, tmp_path / 'm.itp')
    data = sections(tmp_path / 'm.itp')['dihedrals'][:4]
    phi = np.linspace(-np.pi, np.pi, 361)
    actual = sum(float(row[6]) * (1 + np.cos(int(row[7]) * phi - np.deg2rad(float(row[5]))))
                 for row in data)
    expected = 4.184 * sum(k / 2 * (1 + (-1)**(n + 1) * np.cos(n * phi))
                           for n, k in enumerate((1, 2, 3, 4), 1))
    np.testing.assert_allclose(actual, expected, atol=1e-12)


def gromacs_binary():
    if 'GMX_BINARY' in os.environ:
        return shutil.which(os.environ['GMX_BINARY'])
    return shutil.which('gmx_d') or shutil.which('gmx')


def run_tool(args, directory, stdin=None):
    result = subprocess.run(args, cwd=directory, input=stdin, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            env={**os.environ, 'OMP_NUM_THREADS': '1'}, timeout=90)
    assert result.returncode == 0, result.stdout
    return result.stdout


@pytest.mark.integration
@pytest.mark.parametrize('case', ['mixture', 'chain', 'ring'])
def test_grompp_accepts_export(tmp_path, case):
    gmx = gromacs_binary()
    if not gmx:
        pytest.skip('GROMACS is not installed')
    species = ([Water(), Ion('Na'), Water(), Ion('Cl')] if case == 'mixture'
               else [parameterized('CCCC' if case == 'chain' else 'C1CCCCC1')])
    assembled(species).write_gromacs(outdir=tmp_path)
    (tmp_path / 'check.mdp').write_text(
        'integrator = md\nnsteps = 0\ndt = 0.0001\ncutoff-scheme = Verlet\n'
        'coulombtype = Cut-off\nrcoulomb = 1\nrvdw = 1\nrlist = 1\n'
        'constraints = none\npbc = xyz\n')
    run_tool([gmx, 'grompp', '-f', 'check.mdp', '-c', 'system.gro',
              '-p', 'system.top', '-o', 'check.tpr'], tmp_path)
    assert (tmp_path / 'check.tpr').exists()


def four_atom_specie(improper=None, zero_torsion=False):
    from ase import Atoms

    positions = [[0, 0, 0], [1.4, 0.2, 0], [2.1, 1.3, 0.3], [3, 1.4, 1.5]]
    edges = [(0, 1), (1, 2), (2, 3)]
    if improper is not None:
        positions = [[0, 0, 0], [1.4, 0, 0], [-0.5, 1.3, 0], [-0.5, -0.7, 1.1]]
        edges = [(0, 1), (0, 2), (0, 3)]
    atoms = Atoms('C4', positions=positions)
    atoms.info['mdinterface_bonds'] = edges
    atoms.set_array('topology_atom_id', np.arange(4))
    labels = [f'C_{i}' for i in range(4)]
    graph = nx.Graph(edges)
    return Specie(
        atoms, name='FOUR', charges=[-0.2, 0.1, 0.15, -0.05],
        atom_types=[Atom('C', label, eps=0.2, sig=2.5) for label in labels],
        bonds=[Bond(labels[a], labels[b], kr=100, r0=1.5) for a, b in edges],
        angles=[Angle(*(labels[i] for i in ids), kr=20, theta0=109)
                for ids in find_unique_paths_of_length(graph, 2)],
        dihedrals=([Dihedral(*labels, *((0, 0, 0, 0) if zero_torsion else (1, 2, 3, 4)))]
                   if improper is None else []),
        impropers=([Improper(*labels, K=2, d=improper[0], n=improper[1])]
                   if improper is not None else []),
    )


@pytest.mark.integration
@pytest.mark.parametrize('case', ['proper', 'zero', 'improper_positive', 'improper_negative', 'improper_constant'])
def test_gromacs_lammps_single_point_energy_and_forces(tmp_path, case):
    gmx, lmp = gromacs_binary(), shutil.which('lmp')
    if not gmx or not lmp:
        pytest.skip('GROMACS and LAMMPS are required')
    improper = {'improper_positive': (1, 1), 'improper_negative': (-1, 2),
                'improper_constant': (1, 0)}.get(case)
    cell = assembled([four_atom_specie(improper, case == 'zero')])
    compare_engines(cell, tmp_path)


def compare_engines(cell, tmp_path):
    gmx, lmp = gromacs_binary(), shutil.which('lmp')
    if not gmx or not lmp:
        pytest.skip('GROMACS and LAMMPS are required')
    cell.universe.atoms.positions = np.round(cell.universe.atoms.positions, 2)
    cell.write_gromacs(outdir=tmp_path)
    cell.write_lammps(tmp_path / 'system.data')
    (tmp_path / 'check.mdp').write_text(
        'integrator = md\nnsteps = 0\ndt = 0.0001\ncutoff-scheme = Verlet\n'
        'coulombtype = Cut-off\ncoulomb-modifier = None\nvdw-modifier = None\n'
        'rcoulomb = 1\nrvdw = 1\nrlist = 1\nconstraints = none\npbc = xyz\n'
        'nstfout = 1\nnstxout = 1\nnstenergy = 1\n')
    run_tool([gmx, 'grompp', '-f', 'check.mdp', '-c', 'system.gro',
              '-p', 'system.top', '-o', 'check.tpr'], tmp_path)
    run_tool([gmx, 'mdrun', '-s', 'check.tpr', '-deffnm', 'check', '-nt', '1'], tmp_path)
    run_tool([gmx, 'energy', '-f', 'check.edr', '-o', 'energy.xvg'], tmp_path, 'Potential\n0\n')
    energy_lines = [line for line in (tmp_path / 'energy.xvg').read_text().splitlines()
                    if line and line[0] not in '#@']
    gmx_energy = float(energy_lines[-1].split()[1])
    with mda.coordinates.TRR.TRRReader(str(tmp_path / 'check.trr'), convert_units=False) as trajectory:
        gmx_forces = trajectory[0].forces.copy()
    (tmp_path / 'check.in').write_text(
        'units real\natom_style full\npair_style lj/cut/coul/cut 10\n'
        'pair_modify mix geometric\nspecial_bonds lj/coul 0 0 0.5\n'
        'bond_style harmonic\nangle_style harmonic\ndihedral_style opls\n'
        'improper_style cvff\nread_data system.data\n'
        'dump forces all custom 1 forces.dump id fx fy fz\ndump_modify forces sort id format float %.12g\n'
        'run 0\nvariable energy equal pe\nprint "ENERGY ${energy}"\n')
    output = run_tool([lmp, '-in', 'check.in'], tmp_path)
    lmp_energy = float(next(line.split()[1] for line in output.splitlines() if line.startswith('ENERGY ')))
    lines = (tmp_path / 'forces.dump').read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('ITEM: ATOMS')) + 1
    lmp_forces = np.array([[float(v) for v in line.split()[1:]] for line in lines[start:start + len(cell.universe.atoms)]])
    report = {
        'atoms': len(cell.universe.atoms),
        'molecules': len(cell.universe.residues),
        'gromacs_energy_kj_mol': gmx_energy,
        'lammps_energy_kj_mol': lmp_energy * 4.184,
        'energy_absolute_error_kj_mol': abs(gmx_energy - lmp_energy * 4.184),
        'force_max_absolute_error_kj_mol_nm': float(np.max(abs(gmx_forces - lmp_forces * 41.84))),
        'force_rms_error_kj_mol_nm': float(np.sqrt(np.mean((gmx_forces - lmp_forces * 41.84)**2))),
    }
    (tmp_path / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    assert gmx_energy == pytest.approx(lmp_energy * 4.184, abs=0.002, rel=2e-5)
    np.testing.assert_allclose(gmx_forces, lmp_forces * 41.84, atol=0.02, rtol=2e-4)


def test_invalid_second_species_does_not_create_partial_export(tmp_path):
    cell = assembled([Water(), Specie(smiles='CC')])
    with pytest.raises(ValueError, match='Incomplete force field'):
        cell.write_gromacs(outdir=tmp_path / 'output')
    assert not (tmp_path / 'output').exists()


def test_conflicting_atomtype_masses_are_rejected(tmp_path):
    water = Water()
    masses = water.atoms.get_masses()
    masses[-1] = 2
    water.atoms.set_masses(masses)
    with pytest.raises(ValueError, match='Conflicting GROMACS atom type'):
        write_gromacs_itp(water, tmp_path / 'water.itp')
    assert not (tmp_path / 'water.itp').exists()


def test_specie_wrapper_can_omit_atomtypes(tmp_path):
    Water().write_gromacs_itp(tmp_path / 'water.itp', include_atomtypes=False)
    assert 'atomtypes' not in sections(tmp_path / 'water.itp')


@pytest.mark.integration
def test_packed_mixture_exports(tmp_path):
    cell = SimCell(xysize=[25, 25])
    cell.add_solvent(Water(), zdim=25, nsolvent=10,
                     solute=[Ion('Na'), Ion('Cl')], nsolute=[1, 1])
    cell.build()
    cell.write_gromacs(outdir=tmp_path)
    counts = sections(tmp_path / 'system.top')['molecules']
    assert sum(int(count) for _, count in counts) == 12
    assert mda.Universe(str(tmp_path / 'system.gro')).atoms.n_atoms == 32


def common_specie(name):
    from ase import Atoms

    if name == 'water':
        return Water()
    if name in ('Na', 'Cl'):
        return Ion(name)
    data = json.loads((Path(__file__).parent / 'data' / 'gromacs' / f'{name}.json').read_text())
    atoms = Atoms(numbers=data['numbers'], positions=data['positions'])
    atoms.info['mdinterface_bonds'] = data['bonds_graph']
    atoms.set_array('topology_atom_id', np.arange(len(atoms)))
    return Specie(atoms, name=data['name'], charges=data['charges'],
                  atom_types=[Atom(*row) for row in data['atom_types']],
                  bonds=[Bond(*row) for row in data['bonds']],
                  angles=[Angle(*row) for row in data['angles']],
                  dihedrals=[Dihedral(*row) for row in data['dihedrals']],
                  impropers=[Improper(*row) for row in data['impropers']])


@pytest.mark.integration
@pytest.mark.parametrize('names', [
    ('water',), ('benzene',), ('ethanol',), ('methane',),
    ('water', 'water'), ('benzene', 'benzene'), ('ethanol', 'ethanol'),
    ('water', 'ethanol', 'water'), ('water', 'benzene', 'water'),
    ('water', 'Na', 'water', 'Cl', 'water'),
], ids=['water', 'benzene', 'ethanol', 'methane', 'water_dimer', 'benzene_dimer',
        'ethanol_dimer', 'water_ethanol', 'water_benzene', 'water_salt'])
def test_common_molecule_energy_and_forces(tmp_path, names):
    cell = assembled([common_specie(name) for name in names])
    if names == ('benzene', 'benzene'):
        cell.universe.residues[1].atoms.positions += [1.5, 0, 0]
    for i, residue in enumerate(cell.universe.residues):
        residue.atoms.positions += [0, 0.35 * (i % 2), 0.6 * (i % 3)]
    compare_engines(cell, tmp_path)


@pytest.mark.integration
def test_packed_electrolyte_energy_and_forces(tmp_path):
    gmx = gromacs_binary()
    if not gmx or not shutil.which('lmp'):
        pytest.skip('GROMACS and LAMMPS are required')
    version = run_tool([gmx, '--version'], tmp_path)
    if not any(line.split() == ['Precision:', 'double'] for line in version.splitlines()):
        pytest.skip('Double-precision GROMACS is required for the packed force comparison')
    cell = SimCell(xysize=[25, 25])
    cell.add_solvent(Water(), zdim=25, nsolvent=24, solute=[Ion('Na'), Ion('Cl')],
                     nsolute=[2, 2], packmol_tolerance=3.5, seed=42)
    cell.build()
    compare_engines(cell, tmp_path)
