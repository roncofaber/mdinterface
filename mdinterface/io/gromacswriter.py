#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GROMACS include-topology writer for :class:`~mdinterface.core.specie.Specie`.

Writes ``[ atomtypes ]``, ``[ moleculetype ]``, ``[ atoms ]``, ``[ bonds ]``,
``[ angles ]``, and ``[ dihedrals ]`` sections in GROMACS ITP format with the
appropriate unit conversions from LAMMPS *real* units.

Unit conversions
----------------
- Energy  : kcal/mol → kJ/mol  (× 4.184)
- Length  : Å → nm             (× 0.1)
- Bond kr : kcal/mol/Å² → kJ/mol/nm²   (× 418.4)
- Angle kr: kcal/mol/rad² → kJ/mol/rad² (× 4.184)

OPLS-AA proper dihedrals are split into up to four ``funct=9`` terms.
LAMMPS cvff improper dihedrals are mapped to ``funct=4``.
"""

import logging
from itertools import groupby

import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)

#%%

# unit conversion factors
_KCAL_TO_KJ = 4.184
_ANG_TO_NM  = 0.1

# ---------------------------------------------------------------------------
# Column format strings -- header and data share the same field widths so
# that labels sit directly above their values.  Both lines start with a
# 2-character prefix ("; " for headers, "  " for data rows).  Explicit
# 2-space gaps between every column prevent adjacent fields from touching.
#
# The name/type column width (nw) is computed dynamically inside
# write_gromacs_itp from the actual longest atom type label so that short
# molecules don't get unnecessary whitespace.
# ---------------------------------------------------------------------------

def _at_hdr(nw):
    return f"; {{:<{nw}}}  {{:>6}}  {{:>9}}  {{:>8}}  {{:<5}}  {{:>12}}  {{:>14}}\n"

def _at_fmt(nw):
    return f"  {{:<{nw}}}  {{:>6d}}  {{:>9.4f}}  {{:>8.4f}}  {{:<5}}  {{:>12.6f}}  {{:>14.6f}}\n"

def _atoms_hdr(nw):
    return f"; {{:>6}}  {{:<{nw}}}  {{:>5}}  {{:<6}}  {{:<6}}  {{:>4}}  {{:>12}}  {{:>10}}\n"

def _atoms_fmt(nw):
    return f"  {{:>6d}}  {{:<{nw}}}  {{:>5d}}  {{:<6}}  {{:<6}}  {{:>4d}}  {{:>12.6f}}  {{:>10.5f}}\n"

# [ bonds ]  ai(5R)  aj(5R)  funct(5R)  b0(12R,.6f)  kb(16R,.4f)
_BONDS_HDR = "; {:>5}  {:>5}  {:>5}  {:>12}  {:>16}\n"
_BONDS_FMT = "  {:>5d}  {:>5d}  {:>5d}  {:>12.6f}  {:>16.4f}\n"

# [ angles ]  ai(5R)  aj(5R)  ak(5R)  funct(5R)  th0(12R,.4f)  cth(16R,.4f)
_ANGS_HDR  = "; {:>5}  {:>5}  {:>5}  {:>5}  {:>12}  {:>16}\n"
_ANGS_FMT  = "  {:>5d}  {:>5d}  {:>5d}  {:>5d}  {:>12.4f}  {:>16.4f}\n"

# [ dihedrals ]  ai(5R)  aj(5R)  ak(5R)  al(5R)  funct(5R)  phi0(10R,.3f)  kphi(14R,.6f)  n(3R)
_DIHS_HDR  = "; {:>5}  {:>5}  {:>5}  {:>5}  {:>5}  {:>10}  {:>14}  {:>3}\n"
_DIHS_FMT  = "  {:>5d}  {:>5d}  {:>5d}  {:>5d}  {:>5d}  {:>10.3f}  {:>14.6f}  {:>3d}\n"


def _validate_specie(specie):
    specie.validate_force_field()
    if not len(specie.atoms):
        raise ValueError("Cannot export an empty GROMACS molecule.")
    masses = specie.atoms.get_masses()
    if not np.isfinite(masses).all() or (masses <= 0).any():
        raise ValueError("GROMACS atom masses must be finite and positive.")
    for tid in set(specie.dihedrals[1]):
        if specie._dtype[tid].values[4] is not None:
            raise ValueError("GROMACS export requires four-term OPLS dihedrals; A5 is unsupported.")
    for tid in set(specie.impropers[1]):
        _, sign, multiplicity = specie._itype[tid].values
        if sign not in (-1, 1) or multiplicity not in range(7):
            raise ValueError("GROMACS export requires CVFF impropers with d = +/-1 and integer n = 0..6.")


def _atomtype_records(species):
    records = {}
    for specie in species:
        labels, indices = specie.get_atom_types(return_index=True)
        for label, index, atom in zip(labels, indices, specie.atoms):
            parameter = specie._stype[index]
            record = (atom.number, atom.mass, parameter.sig, parameter.eps)
            if label in records and records[label] != record:
                raise ValueError(f"Conflicting GROMACS atom type {label!r}.")
            records[label] = record
    return records


def _write_atomtypes(stream, records):
    nw = max(4, max(map(len, records)))
    stream.write("[ atomtypes ]\n")
    stream.write(_at_hdr(nw).format(
        "name", "atnum", "mass", "charge", "ptype", "sigma(nm)", "eps(kJ/mol)"))
    for label, (number, mass, sigma, epsilon) in records.items():
        stream.write(_at_fmt(nw).format(
            label, number, mass, 0.0, "A", sigma * _ANG_TO_NM, epsilon * _KCAL_TO_KJ))
    stream.write("\n")


def _species_signature(specie):
    interactions = []
    for name, key in (("bonds", "_btype"), ("angles", "_atype"),
                      ("dihedrals", "_dtype"), ("impropers", "_itype")):
        indices, types = getattr(specie, name)
        parameters = getattr(specie, key)
        interactions.append(tuple((tuple(ids), tuple(parameters[tid].values))
                                  for ids, tid in zip(indices, types)))
    return (tuple(specie.get_atom_types()), tuple(specie._sids),
            tuple(specie.atoms.numbers), tuple(specie.atoms.get_masses()),
            tuple(specie.charges), tuple(_atomtype_records([specie]).items()),
            tuple(interactions))


def _prepare_species(species, universe):
    unique = {}
    signatures = {}
    for specie in species:
        _validate_specie(specie)
        signature = _species_signature(specie)
        if specie.resname in signatures and signatures[specie.resname] != signature:
            raise ValueError(f"Conflicting GROMACS molecule definitions for {specie.resname!r}; use distinct species names.")
        unique[specie.resname] = specie
        signatures[specie.resname] = signature
    _atomtype_records(unique.values())
    for residue in universe.residues:
        if residue.resname not in unique:
            raise ValueError(f"No GROMACS molecule definition for {residue.resname!r}.")
        specie = unique[residue.resname]
        atoms = residue.atoms
        if (len(atoms) != len(specie.atoms)
                or not np.array_equal(atoms.types, specie.get_atom_types())
                or not np.allclose(atoms.charges, specie.charges, rtol=0, atol=1e-8)
                or not np.allclose(atoms.masses, specie.atoms.get_masses(), rtol=0, atol=1e-8)):
            raise ValueError(f"Assembled residue {residue.resname!r} does not match its GROMACS molecule definition.")
    return list(unique.values())


def write_gromacs_itp(specie, filename=None, *, include_atomtypes=True):
    """
    Write a GROMACS include topology (.itp) file for a Specie.

    .. warning::
        GROMACS output is experimental and assumes harmonic bonds and angles,
        four-term OPLS proper torsions, CVFF impropers, geometric LJ mixing,
        and 0.5 LJ/Coulomb 1-4 scaling. Constraints are not generated.

    Parameters
    ----------
    specie : Specie
        The molecular species to write.
    filename : str, optional
        Output filename. Defaults to ``{resname}.itp``.
    include_atomtypes : bool, default True
        Include atom-type definitions. For multi-species systems, set False
        and pass the species to :func:`write_gromacs_top` so all atom types
        precede all molecule definitions.

    Raises
    ------
    ValueError
        If parameters are incomplete, conflicting, or unsupported.
    """
    _validate_specie(specie)
    records = _atomtype_records([specie])
    if filename is None:
        filename = f"{specie.resname}.itp"

    masses   = specie.atoms.get_masses()
    charges  = specie.atoms.get_initial_charges()
    atom_types = specie.get_atom_types()

    bond_idxs,  bond_tids  = specie.bonds
    angle_idxs, angle_tids = specie.angles
    dih_idxs,   dih_tids   = specie.dihedrals
    imp_idxs,   imp_tids   = specie.impropers

    # column width for name/type: fit the longest label, at least "name" (4)
    nw = max(max(len(t) for t in atom_types), 4)

    ATOMS_HDR = _atoms_hdr(nw)
    ATOMS_FMT = _atoms_fmt(nw)

    with open(filename, "w") as f:
        f.write("; GROMACS ITP file generated by mdinterface\n")
        f.write(f"; Molecule: {specie.resname}\n\n")

        if include_atomtypes:
            _write_atomtypes(f, records)

        # ------------------------------------------------------------------
        # [ moleculetype ]
        # ------------------------------------------------------------------
        f.write("[ moleculetype ]\n")
        f.write("; name      nrexcl\n")
        f.write(f"  {specie.resname:<10s}  3\n\n")

        # ------------------------------------------------------------------
        # [ atoms ]
        # ------------------------------------------------------------------
        f.write("[ atoms ]\n")
        f.write(ATOMS_HDR.format(
            "nr", "type", "resnr", "residu", "atom", "cgnr", "charge", "mass"))
        for ii, (atype, charge, mass, sid) in enumerate(
                zip(atom_types, charges, masses, specie._sids)):
            atom_name = sid.split("_")[0]
            f.write(ATOMS_FMT.format(
                ii+1, atype, 1, specie.resname, atom_name, ii+1, charge, mass))
        f.write("\n")

        # ------------------------------------------------------------------
        # [ bonds ]
        # ------------------------------------------------------------------
        if bond_idxs:
            f.write("[ bonds ]\n")
            f.write(_BONDS_HDR.format("ai", "aj", "funct", "b0(nm)", "kb(kJ/mol/nm2)"))
            for (ai, aj), tid in zip(bond_idxs, bond_tids):
                bond = specie._btype[tid]
                b0 = bond.r0 * _ANG_TO_NM
                # LAMMPS: E = K*(r-r0)^2  |  GROMACS: V = (kb/2)*(r-r0)^2  => kb = 2*K
                kb = 2 * bond.kr * _KCAL_TO_KJ / _ANG_TO_NM**2
                f.write(_BONDS_FMT.format(ai+1, aj+1, 1, b0, kb))
            f.write("\n")

        graph = nx.Graph()
        graph.add_edges_from(bond_idxs)
        pairs = sorted((a, b) for a, distances in nx.all_pairs_shortest_path_length(graph, cutoff=3)
                       for b, distance in distances.items() if a < b and distance == 3)
        if pairs:
            f.write("[ pairs ]\n; ai  aj  funct\n")
            for ai, aj in pairs:
                f.write(f"  {ai + 1}  {aj + 1}  1\n")
            f.write("\n")

        # ------------------------------------------------------------------
        # [ angles ]
        # ------------------------------------------------------------------
        if angle_idxs:
            f.write("[ angles ]\n")
            f.write(_ANGS_HDR.format("ai", "aj", "ak", "funct", "th0(deg)", "cth(kJ/mol/rad2)"))
            for (ai, aj, ak), tid in zip(angle_idxs, angle_tids):
                angle = specie._atype[tid]
                th0 = angle.theta0
                # LAMMPS: E = K*(θ-θ0)^2  |  GROMACS: V = (kt/2)*(θ-θ0)^2  => kt = 2*K
                cth = 2 * angle.kr * _KCAL_TO_KJ
                f.write(_ANGS_FMT.format(ai+1, aj+1, ak+1, 1, th0, cth))
            f.write("\n")

        # ------------------------------------------------------------------
        # [ dihedrals ] -- proper (OPLS Fourier -> funct=9)
        #
        # LAMMPS: E = (K1/2)(1+cos φ) + (K2/2)(1-cos 2φ)
        #           + (K3/2)(1+cos 3φ) + (K4/2)(1-cos 4φ)
        # GROMACS funct=9: V = kphi * (1 + cos(n*phi - phi0))
        #   K1 -> n=1, phi0=0,   kphi = K1/2
        #   K2 -> n=2, phi0=180, kphi = K2/2
        #   K3 -> n=3, phi0=0,   kphi = K3/2
        #   K4 -> n=4, phi0=180, kphi = K4/2
        # ------------------------------------------------------------------
        if dih_idxs:
            f.write("[ dihedrals ]\n")
            f.write("; proper dihedrals -- OPLS Fourier series, funct=9\n")
            f.write(_DIHS_HDR.format("ai", "aj", "ak", "al", "funct", "phi0(deg)", "kphi(kJ/mol)", "n"))
            opls_map = [(0, 1, 0.0), (1, 2, 180.0),
                        (2, 3, 0.0), (3, 4, 180.0)]
            for (ai, aj, ak, al), tid in zip(dih_idxs, dih_tids):
                dih  = specie._dtype[tid]
                vals = dih.values  # [K1, K2, K3, K4, K5?]
                for vi, n, phi0 in opls_map:
                    K = vals[vi]
                    if K == 0.0:
                        continue
                    kphi = K / 2.0 * _KCAL_TO_KJ
                    f.write(_DIHS_FMT.format(ai+1, aj+1, ak+1, al+1, 9, phi0, kphi, int(n)))
            f.write("\n")

        # ------------------------------------------------------------------
        # [ dihedrals ] -- improper (LAMMPS cvff -> funct=4)
        #
        # LAMMPS: E = K(1 + d*cos(n*phi))
        # GROMACS funct=4: V = kphi(1 + cos(n*psi - psi0))
        #   d= 1 -> psi0=0
        #   d=-1 -> psi0=180
        # ------------------------------------------------------------------
        if imp_idxs:
            f.write("[ dihedrals ]\n")
            f.write("; improper dihedrals, funct=4\n")
            f.write(_DIHS_HDR.format("ai", "aj", "ak", "al", "funct", "phi0(deg)", "kphi(kJ/mol)", "n"))
            for (ai, aj, ak, al), tid in zip(imp_idxs, imp_tids):
                imp      = specie._itype[tid]
                K, d, n  = imp.values
                phi0 = 0.0 if d == 1 else 180.0
                kphi = K * _KCAL_TO_KJ
                f.write(_DIHS_FMT.format(ai+1, aj+1, ak+1, al+1, 4, phi0, kphi, int(n)))
            f.write("\n")



def write_gromacs_top(universe, itp_files, filename="system.top",
                      system_name="MD System", *, species=None):
    """
    Write a GROMACS system topology (.top) file.

    Parameters
    ----------
    universe : mda.Universe
        The assembled system (after SimCell.build()).  Used to determine
        molecule counts from ``universe.residues``.
    itp_files : list of str
        Basenames of the per-species ITP files to ``#include``.
    filename : str, optional
        Output filename. Default ``"system.top"``.
    system_name : str, optional
        Title written in the ``[ system ]`` section.
    species : sequence of Specie, optional
        Species whose atom types are written before the molecule includes.
        These ITP files must be written with ``include_atomtypes=False``.
        Also validates residue identity against the assembled universe.
    """
    records = None
    if species is not None:
        records = _atomtype_records(_prepare_species(species, universe))
    runs = [(name, sum(1 for _ in group))
            for name, group in groupby(universe.residues.resnames)]

    with open(filename, "w") as f:
        f.write("; GROMACS topology file generated by mdinterface\n\n")

        # OPLS-AA defaults: LJ (nbfunc=1), geometric combining (comb-rule=3)
        f.write("[ defaults ]\n")
        f.write("; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n")
        f.write("  1       3          yes        0.5      0.5\n\n")

        if records:
            _write_atomtypes(f, records)

        for itp in itp_files:
            f.write(f'#include "{itp}"\n')
        f.write("\n")

        f.write("[ system ]\n")
        f.write(f"{system_name}\n\n")

        f.write("[ molecules ]\n")
        f.write("; {:<14}  {:>6}\n".format("molecule", "nmols"))
        for resname, count in runs:
            f.write(f"  {resname:<14s}  {count:>6d}\n")

