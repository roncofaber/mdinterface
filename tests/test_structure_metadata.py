import hashlib
import json
import logging
import shutil

import MDAnalysis as mda
import numpy as np
import pytest

from mdinterface import SimCell
from mdinterface.database import Water
from mdinterface.io.lammpswriter import DATAWriter
from mdinterface.io.structure_metadata import lammps_metadata


@pytest.fixture(autouse=True)
def restore_logging():
    logger = logging.getLogger("mdinterface")
    state = logger.level, logger.propagate, list(logger.handlers)
    yield
    logger.setLevel(state[0])
    logger.propagate = state[1]
    logger.handlers[:] = state[2]


def test_exported_sparse_types_and_molecule_ids(tmp_path):
    universe = mda.Universe.empty(3, n_residues=1, atom_resindex=[0, 0, 0], trajectory=True)
    for name, values in {"types": ["O_z", "H_a", "H_a"], "elements": ["O", "H", "H"],
                         "masses": [15.999, 1.008, 1.008], "charges": [-.8, .4, .4],
                         "resnames": ["WATR"], "resids": [42]}.items():
        universe.add_TopologyAttr(name, values)
    universe.add_bonds([(0, 1), (0, 2)], types=[8, 17])
    universe.add_angles([(1, 0, 2)], types=[29])
    universe.atoms.positions = [[3, 3, 3], [4, 3, 3], [3, 4, 3]]
    universe.dimensions = [10, 11, 12, 90, 90, 80]
    data = tmp_path / "data.lammps"
    with DATAWriter(str(data)) as writer:
        writer.write(universe.atoms)
    description = lammps_metadata(data, universe, "full", "test")
    assert [atom["type_id"] for atom in description["atoms"]] == [2, 1, 1]
    assert {atom["molecule_id"] for atom in description["atoms"]} == {0}
    assert [item["type_id"] for item in description["topology"]["bonds"]] == [1, 2]
    assert description["topology"]["angles"][0]["type_id"] == 1
    assert description["groups"] == {"WATR": [1, 2, 3]}
    assert description["cell"]["tilt_A"][0] != 0
    assert description["data_file"]["sha256"] == hashlib.sha256(data.read_bytes()).hexdigest()


@pytest.mark.integration
@pytest.mark.parametrize("atomic", [False, True])
def test_simcell_optional_metadata(tmp_path, monkeypatch, atomic):
    if not shutil.which("packmol"):
        pytest.skip("PACKMOL required")
    monkeypatch.chdir(tmp_path)
    box = SimCell(xysize=[15, 15], verbose=False)
    box.add_solvent(Water(), nsolvent=4, zdim=15, seed=1234)
    box.build(padding=0)
    box.write_lammps("data.lammps", atom_style="atomic" if atomic else "full",
                     write_coeff=not atomic, metadata="structure.json")
    result = json.loads((tmp_path / "structure.json").read_text())
    assert result["schema"] == "mdinterface.structure"
    assert len(result["atoms"]) == 12
    assert result["simulation_choices"]["constraints"] is None
    if atomic:
        assert result["topology"]["bonds"] == []
        assert all(atom["molecule_id"] is None and atom["charge"] is None for atom in result["atoms"])
    else:
        exported = mda.Universe("data.lammps", format="DATA")
        np.testing.assert_array_equal([atom["type_id"] for atom in result["atoms"]], exported.atoms.types.astype(int))
        assert len(result["topology"]["bonds"]) == 8
        assert result["coefficients"]["Bond Coeffs"]
    before = (tmp_path / "data.lammps").read_bytes()
    with pytest.raises(FileExistsError): box.write_lammps("data.lammps", metadata="structure.json")
    assert (tmp_path / "data.lammps").read_bytes() == before
    with pytest.raises(ValueError, match="differ"):
        box.write_lammps("same.lammps", metadata="same.lammps")
    box.write_lammps("legacy.lammps")
    assert not (tmp_path / "legacy.json").exists()
    if not atomic:
        assert (tmp_path / "legacy.lammps").read_bytes() == before


@pytest.mark.integration
def test_mixed_electrode_electrolyte_export(tmp_path, monkeypatch):
    from mdinterface.database import Ion, Metal111

    monkeypatch.chdir(tmp_path)
    box = SimCell(xysize=[15, 15], verbose=False)
    box.add_slab(Metal111('Au'), nlayers=1)
    for _ in range(2):
        box.add_solvent(Water(), nsolvent=4, solute=[Ion('Na'), Ion('Cl')],
                        nsolute=[1, 1], zdim=15)
    box.build(center=True)
    box.write_lammps('mixed.lammps', metadata='mixed.json')
    result = json.loads((tmp_path / 'mixed.json').read_text())
    exported = mda.Universe('mixed.lammps', format='DATA')
    assert len(exported.atoms) == len(box.universe.atoms)
    assert {atom['element'] for atom in result['atoms']} == {'Au', 'O', 'H', 'Na', 'Cl'}
    np.testing.assert_allclose(exported.atoms.charges.sum(), 0, atol=1e-6)
    for name, expected in [('bonds', 16), ('angles', 8)]:
        assert len(result['topology'][name]) == len(getattr(exported, name)) == expected
        section = name[:-1].title() + ' Coeffs'
        assert {item['type_id'] for item in result['topology'][name]} == {
            item['type_id'] for item in result['coefficients'][section]
        }
    assert {atom['type_id'] for atom in result['atoms']} == {
        item['type_id'] for item in result['coefficients']['Pair Coeffs']
    }
    header = (tmp_path / 'mixed.lammps').read_text().split('Masses')[0]
    assert f'{len(exported.atoms):12d}  atoms' in header
    assert f'{16:12d}  bonds' in header
    assert f'{8:12d}  angles' in header


def test_torsion_export_matches_coefficients_and_metadata(tmp_path):
    from mdinterface.core.topology import Atom, Dihedral, Improper
    from mdinterface.io.lammpswriter import write_lammps_coefficients

    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0]*4, trajectory=True)
    for name, values in {'types': ['C_TEST']*4, 'elements': ['C']*4, 'masses': [12.011]*4,
                         'charges': [0.0]*4, 'resnames': ['TEST'], 'resids': [1]}.items():
        universe.add_TopologyAttr(name, values)
    universe.atoms.positions = [[1, 1, 1], [2, 1, 1], [2, 2, 1], [3, 2, 2]]
    universe.dimensions = [10, 10, 10, 90, 90, 90]
    universe.add_bonds([])
    universe.add_angles([])
    universe.add_dihedrals([(0, 1, 2, 3)], types=[17])
    universe.add_impropers([(1, 0, 2, 3)], types=[29])
    dihedral = Dihedral('C', 'C', 'C', 'C', A1=1, A2=2, A3=3, A4=4)
    improper = Improper('C', 'C', 'C', 'C', K=2, d=-1, n=2)
    atom = Atom('C', eps=0.1, sig=3)
    atom.set_resname('TEST')
    dihedral.set_id(17)
    improper.set_id(29)
    path = tmp_path / 'torsions.lammps'
    with DATAWriter(str(path)) as writer:
        writer.write(universe.atoms)
    with path.open('a') as output:
        output.write('\n')
        write_lammps_coefficients(universe, {'atoms': [atom],
                                           'dihedrals': [dihedral], 'impropers': [improper]}, fout=output)
    result = lammps_metadata(path, universe, 'full', 'test')
    exported = mda.Universe(str(path), format='DATA')
    for name, section, expected in [('dihedrals', 'Dihedral Coeffs', [1, 2, 3, 4]),
                                     ('impropers', 'Improper Coeffs', [2, -1, 2])]:
        assert len(getattr(exported, name)) == 1
        assert result['topology'][name][0]['type_id'] == 1
        assert result['coefficients'][section][0]['type_id'] == 1
        np.testing.assert_allclose([float(x) for x in result['coefficients'][section][0]['tokens']], expected)
        np.testing.assert_array_equal(getattr(exported, name).indices[0] + 1,
                                      result['topology'][name][0]['atom_ids'])
    header = path.read_text().split('Masses')[0]
    for kind in ('dihedral', 'improper'):
        assert f'{1:12d}  {kind}s' in header
        assert f'{1:12d}  {kind} types' in header
