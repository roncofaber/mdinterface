"""Machine-readable structural facts from the final LAMMPS export."""

import hashlib
from pathlib import Path


def lammps_metadata(filename, system, atom_style, version):
    path = Path(filename)
    names = {"Atoms", "Masses", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers",
             "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs"}
    sections = {name: [] for name in names}
    bounds, tilt, section = {}, [0.0, 0.0, 0.0], None
    for line in path.read_text().splitlines():
        body, _, annotation = line.partition("#")
        body = body.strip()
        if body in names:
            section = body
        elif body:
            fields = body.split()
            if fields[-2:] in (["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]):
                bounds[fields[-2][0]] = [float(value) for value in fields[:2]]
            elif fields[-3:] == ["xy", "xz", "yz"]:
                tilt = [float(value) for value in fields[:3]]
            elif section:
                sections[section].append((fields, annotation.strip()))
    atoms, groups, molecules = [], {}, {}
    masses = {int(row[0]): float(row[1]) for row, _ in sections["Masses"]}
    for row, _ in sections["Atoms"]:
        atom_id = int(row[0])
        original = system.atoms[atom_id - 1]
        molecule_id = int(row[1]) if atom_style == "full" else None
        type_id = int(row[2] if atom_style == "full" else row[1])
        entry = dict(id=atom_id, type_id=type_id, molecule_id=molecule_id,
                     source_type_label=str(original.type), species=str(original.resname),
                     element=str(original.element), mass=masses[type_id],
                     charge=float(row[3]) if atom_style == "full" else None,
                     position_A=[float(value) for value in row[-3:]])
        atoms.append(entry)
        groups.setdefault(entry["species"], []).append(atom_id)
        if molecule_id is not None:
            molecules.setdefault(molecule_id, []).append(atom_id)
    topology = {name.lower(): [dict(id=int(row[0]), type_id=int(row[1]), atom_ids=[int(value) for value in row[2:]])
                               for row, _ in sections[name]]
                for name in ("Bonds", "Angles", "Dihedrals", "Impropers")}
    coefficients = {name: [dict(type_id=int(row[0]), tokens=row[1:], annotation=annotation)
                          for row, annotation in sections[name]]
                    for name in sorted(names) if name.endswith("Coeffs")}
    types = [dict(id=type_id, mass=masses[type_id],
                  source_labels=sorted({atom["source_type_label"] for atom in atoms if atom["type_id"] == type_id}),
                  elements=sorted({atom["element"] for atom in atoms if atom["type_id"] == type_id}))
             for type_id in sorted(masses)]
    return dict(schema="mdinterface.structure", schema_version=1,
                generator=dict(package="mdinterface", version=version),
                data_file=dict(name=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), atom_style=atom_style),
                units=dict(length="angstrom", mass="g/mol", charge="e", coefficients="LAMMPS real"),
                cell=dict(bounds_A=bounds, tilt_A=tilt, boundary_conditions=None),
                atoms=atoms, atom_types=types, groups=groups,
                molecules=[dict(id=key, atom_ids=value) for key, value in sorted(molecules.items())],
                topology=topology, coefficients=coefficients, parameter_provenance=None,
                simulation_choices=dict(pair_style=None, mixing=None, constraints=None, ensemble=None))
