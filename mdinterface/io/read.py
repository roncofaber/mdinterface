#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Readers for LAMMPS data and dump files.

Wraps ASE's LAMMPS I/O with mdinterface topology reconstruction so that
:class:`~mdinterface.core.specie.Specie` objects can be initialised directly
from LAMMPS data files, and trajectory frames can be loaded for use with the
``hijack`` mode of :meth:`~mdinterface.build.builder.SimCell.build`.
"""

import collections
import logging
from io import StringIO

import numpy as np

logger = logging.getLogger(__name__)

import ase
import ase.visualize
import ase.io.lammpsdata
import ase.io.lammpsrun
from mdinterface.core.topology import Atom, Bond, Angle, Dihedral, Improper
from mdinterface.utils.auxiliary import mass2symbol

#%%

def read_lammps_data_file(filename, pbc=False, ato_start_idx=0, is_snippet=False):
    logger.debug("Reading LAMMPS data file: %s", filename)
    system = ase.io.lammpsdata.read_lammps_data(filename)
    
    if not pbc:
        system.set_pbc(False)
        system.set_cell(None)

    symbols = system.get_chemical_symbols()
    unique_symbols = set(symbols)
    
    pair_coeff     = []
    bond_coeff     = []
    angle_coeff    = []
    dihedral_coeff = []
    improper_coeff = []

    atoms     = []
    atosym    = []
    bonds     = []
    angles    = []
    dihedrals = []
    impropers = []

    section_map = {
        "masses"     :     "masses",
        'pair coeffs':     'pair_coeff',
        'bond coeffs':     'bond_coeff',
        'angle coeffs':    'angle_coeff',
        'dihedral coeffs': 'dihedral_coeff',
        'improper coeffs': 'improper_coeff',
        'atoms':            'atoms',
        'bonds':            'bonds',
        'angles':           'angles',
        'dihedrals':        'dihedrals',
        'impropers':        'impropers'
    }

    with open(filename, "r") as file:
        to_read = None

        for line in file:
            line = line.strip()

            if not line:
                continue

            section_key = line.lower()

            if section_key in section_map:
                to_read = section_map[section_key]
                continue

            line = line.split('#')[0].split()

            if to_read == "masses":
                mass = float(line[1])
                atosym.append(mass2symbol(mass, unique_symbols))
                
            elif to_read == "pair_coeff":
                eps = float(line[1])
                sig = float(line[2])
                pair_coeff.append([eps, sig])

            elif to_read == "bond_coeff":
                kr = float(line[1])
                r0 = float(line[2])
                bond_coeff.append([kr, r0])

            elif to_read == "angle_coeff":
                kr = float(line[1])
                theta0 = float(line[2])
                angle_coeff.append([kr, theta0])

            elif to_read == "dihedral_coeff":
                values = [float(ii) for ii in line[1:]]
                dihedral_coeff.append(values)

            elif to_read == "improper_coeff":
                values = [float(line[1]), int(line[2]), int(line[3])]
                improper_coeff.append(values)

            elif to_read == "atoms":
                atoidx = ato_start_idx + (int(line[0]) - 1)
                atotyp = int(line[2]) - 1
                eps, sig = pair_coeff[atotyp]
                symbol   = atosym[atotyp]
                
                if not is_snippet:
                    atom_name = f"{symbol}_{str(atoidx).zfill(3)}"
                else:
                    atom_name = f"{symbol}_{str(atoidx).zfill(3)}sn"
                
                atom = Atom(symbol, atom_name, eps, sig, atoid=atotyp)
                atoms.append(atom)

            elif to_read == "bonds":
                idx1 = int(line[2]) - 1
                idx2 = int(line[3]) - 1
                kr, r0 = bond_coeff[int(line[1]) - 1]
                bond = Bond(atoms[idx1].label, atoms[idx2].label, kr, r0)
                bonds.append(bond)

            elif to_read == "angles":
                idx1 = int(line[2]) - 1
                idx2 = int(line[3]) - 1
                idx3 = int(line[4]) - 1
                kr, theta0 = angle_coeff[int(line[1]) - 1]
                angle = Angle(atoms[idx1].label, atoms[idx2].label, atoms[idx3].label, kr, theta0)
                angles.append(angle)

            elif to_read == "dihedrals":
                idxs = [int(ii) - 1 for ii in line[2:]]
                values = dihedral_coeff[int(line[1]) - 1]
                dihedral = Dihedral(atoms[idxs[0]].label, atoms[idxs[1]].label,
                                    atoms[idxs[2]].label, atoms[idxs[3]].label, *values)
                dihedrals.append(dihedral)

            elif to_read == "impropers":
                idxs = [int(ii) - 1 for ii in line[2:]]
                values = improper_coeff[int(line[1]) - 1]
                improper = Improper(atoms[idxs[0]].label, atoms[idxs[1]].label,
                                    atoms[idxs[2]].label, atoms[idxs[3]].label,
                                    K=values[0], d=values[1], n=values[2])
                impropers.append(improper)

    system.new_array("stype", np.array(atoms))
    logger.debug("  └─> %d atoms, %d bonds, %d angles, %d dihedrals, %d impropers",
                 len(atoms), len(bonds), len(angles), len(dihedrals), len(impropers))
    return system, atoms, bonds, angles, dihedrals, impropers


# ---------------------------------------------------------------------------
# LAMMPS dump frame reader
# ---------------------------------------------------------------------------

def _iter_lammps_dump_frames(filename):
    """Yield raw text chunks, one per frame, from a LAMMPS dump file."""
    with open(filename) as fh:
        lines = []
        for line in fh:
            if line.startswith("ITEM: TIMESTEP") and lines:
                yield "".join(lines)
                lines = []
            lines.append(line)
        if lines:
            yield "".join(lines)


def _read_last_frame_text(filename):
    """Return the raw text of the last frame by seeking backwards from EOF.

    Scans from the end of the file in 64 KiB chunks until the last
    ``ITEM: TIMESTEP`` marker is found.  This is O(frame size) rather than
    O(file size), making it orders of magnitude faster on multi-GB trajectories.
    """
    marker = b"ITEM: TIMESTEP"
    chunk_size = 1 << 16  # 64 KiB
    with open(filename, "rb") as fh:
        fh.seek(0, 2)
        pos = fh.tell()
        buf = b""
        while True:
            step = min(chunk_size, pos)
            pos -= step
            fh.seek(pos)
            buf = fh.read(step) + buf
            idx = buf.rfind(marker)
            if idx != -1:
                return buf[idx:].decode()
            if pos == 0:
                break
    return buf.decode()


def read_lammps_nth_frame(filename, frame=-1):
    """
    Read a single frame from a LAMMPS dump file without loading the full trajectory.

    Parameters
    ----------
    filename : str
        Path to the LAMMPS dump file.
    frame : int
        Frame index. ``0`` = first frame, ``-1`` = last frame (default),
        ``-2`` = second-to-last, etc.
        Reading the last frame (``frame=-1``) seeks from the end of the file
        and is O(frame size), making it fast even for multi-GB trajectories.

    Returns
    -------
    ase.Atoms
        The requested frame as an ASE Atoms object with cell and PBC set.

    Raises
    ------
    IndexError
        If the requested frame index is out of range.
    """
    logger.debug("Reading frame %d from dump: %s", frame, filename)
    if frame == -1:
        # Fast path: seek from EOF -- O(frame size) regardless of file size.
        return ase.io.lammpsrun.read_lammps_dump_text(
            StringIO(_read_last_frame_text(filename))
        )
    elif frame >= 0:
        for ii, chunk in enumerate(_iter_lammps_dump_frames(filename)):
            if ii == frame:
                return ase.io.lammpsrun.read_lammps_dump_text(StringIO(chunk))
        raise IndexError(
            f"Dump file '{filename}' has fewer than {frame + 1} frame(s)."
        )
    else:
        # Negative index: keep a rolling buffer of the last abs(frame) chunks.
        buf = collections.deque(maxlen=-frame)
        for chunk in _iter_lammps_dump_frames(filename):
            buf.append(chunk)
        if len(buf) < -frame:
            raise IndexError(
                f"Dump file '{filename}' has fewer than {-frame} frame(s)."
            )
        return ase.io.lammpsrun.read_lammps_dump_text(StringIO(buf[0]))
