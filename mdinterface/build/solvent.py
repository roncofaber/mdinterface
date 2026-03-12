#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
solvent.py — solvent box assembly.

Contains parameter validation, molecule-count computation, and the
high-level ``make_solvent_box`` entry point.  Low-level PACKMOL
orchestration lives in ``box.py`` (``populate_box``).
"""

from typing import Any, Dict, List, Optional, Tuple, Union

import warnings

import numpy as np
import MDAnalysis as mda
from ase import units

from mdinterface.build.box import populate_box
from mdinterface.build.continuum2sim import discretize_concentration

import logging

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _validate_solvent_box_parameters(
    nsolute: Optional[Union[int, List[int]]],
    concentration: Optional[float],
    conmodel: Optional[Dict[int, Tuple[List[float], List[float]]]],
    solute: Optional[List[Any]],
    solvents: List[Any],
    density: Optional[float],
    nsolvent=None,
    ratio: Optional[List[float]] = None,
) -> None:
    """Validate parameter combinations for make_solvent_box."""

    if nsolute is not None and concentration is not None:
        raise ValueError("Cannot specify both 'nsolute' and 'concentration'. Use one or the other.")

    if conmodel is not None and (solute is None or len(solute) == 0):
        raise ValueError("When using 'conmodel', 'solute' must be provided and non-empty.")

    if isinstance(nsolute, (list, tuple)) and solute is not None:
        if len(nsolute) != len(solute):
            raise ValueError(f"Length of 'nsolute' ({len(nsolute)}) must match number of solute species ({len(solute)}).")

    if ratio is not None:
        if len(ratio) != len(solvents):
            raise ValueError(f"Length of 'ratio' ({len(ratio)}) must match number of solvent species ({len(solvents)}).")
        if isinstance(nsolvent, (list, tuple)):
            raise ValueError("Cannot use 'ratio' with 'nsolvent' as a list; use ratio+density or ratio+nsolvent(int).")
        if density is None and nsolvent is None:
            raise ValueError("'ratio' requires either 'density' or 'nsolvent' (total count).")

    if isinstance(nsolvent, (list, tuple)) and solvents:
        if len(nsolvent) != len(solvents):
            raise ValueError(f"Length of 'nsolvent' ({len(nsolvent)}) must match number of solvent species ({len(solvents)}).")

    if len(solvents) > 1 and ratio is None and not isinstance(nsolvent, (list, tuple)) and density is not None and nsolvent is None:
        raise ValueError("For a solvent mixture, specify 'nsolvent' (list), 'ratio'+'density', or 'ratio'+'nsolvent'.")

    if density is not None and nsolvent is not None and ratio is None:
        warnings.warn(
            "Both 'density' and 'nsolvent' are specified. Using 'nsolvent' and ignoring 'density'.",
            UserWarning, stacklevel=4,
        )

    if (density is not None or nsolvent is not None) and not solvents:
        warnings.warn("Density or nsolvent specified but no solvent provided. Will be ignored.",
                      UserWarning, stacklevel=4)

    if not solvents and (solute is None or len(solute) == 0):
        warnings.warn("No solvent or solute specified. Empty box will be created.",
                      UserWarning, stacklevel=3)


# ---------------------------------------------------------------------------
# Solute placement — internal helpers
# ---------------------------------------------------------------------------

def _place_sp(sp, volume, coords, radii, zpos=None, max_attempts=100):
    """Try to find a non-overlapping random position for one molecule."""
    radius = sp.estimate_specie_radius()
    for _ in range(max_attempts):
        new_coord = radius + 1 + np.random.rand(3) * (volume - 2 * (radius + 1))
        if zpos is not None:
            new_coord[2] = zpos
        if not coords or np.all(
            np.linalg.norm(np.array(coords) - new_coord, axis=1)
            >= np.array(radii) + radius + 1
        ):
            return new_coord
    logger.warning("Failed to place %s after %d attempts", sp, max_attempts)
    return None


def _place_conmodel(solute, conmodel, volume):
    """Fixed-coordinate placement from a spatial concentration profile."""
    coords = []
    radii = []
    instructions = []
    for cc, sp in enumerate(solute):
        z_coords, conc_profile = conmodel[cc]
        z_positions = discretize_concentration(sp, conc_profile, z_coords, volume)
        for zpos in z_positions:
            new_coord = _place_sp(sp, volume, coords, radii, zpos=zpos)
            if new_coord is not None:
                coords.append(new_coord)
                radii.append(sp.estimate_specie_radius())
                instructions.append((sp.to_universe(), new_coord, "fixed"))
    return instructions


# ---------------------------------------------------------------------------
# Solute placement — public API
# ---------------------------------------------------------------------------

def populate_solutes(solute, nsolute, volume, solute_pos=None, conmodel=None):
    """
    Build PACKMOL/fixed-coord placement instructions for solute molecules.

    Parameters
    ----------
    solute : list of Specie
        Solute species.
    nsolute : int or list of int
        Number of molecules per species.
    volume : list of float
        Box dimensions [x, y, z] in Å.
    solute_pos : str or None
        Placement mode:

        - ``None`` / ``"packmol"`` — PACKMOL random placement in the full box
          (default).
        - ``"left"``   — PACKMOL random placement in the left half (z ≤ z/2).
        - ``"right"``  — PACKMOL random placement in the right half (z ≥ z/2).
        - ``"center"`` — each molecule fixed at the box centre.
        - ``"box"``    — deprecated alias for ``"packmol"``.

    conmodel : dict or None
        Spatially varying concentration model; takes precedence over
        *solute_pos* when provided.
    """
    import warnings

    volume = np.array(volume)
    # Padded PACKMOL bounding box [xmin, ymin, zmin, xmax, ymax, zmax]
    box = [1.0, 1.0, 1.0, volume[0] - 1.0, volume[1] - 1.0, volume[2] - 1.0]

    # conmodel takes precedence — fixed z-coordinate per molecule
    if conmodel is not None:
        assert len(conmodel) == len(solute), "Need one profile per solute species"
        return _place_conmodel(solute, conmodel, volume)

    # center — fixed at box centre (one fixed instruction per molecule)
    if solute_pos == "center":
        coord = list(volume / 2)
        return [
            (sp.to_universe(), coord, "fixed")
            for cc, sp in enumerate(solute)
            for _ in range(nsolute if isinstance(nsolute, int) else nsolute[cc])
        ]

    # deprecated alias
    if solute_pos == "box":
        warnings.warn(
            "solute_pos='box' is deprecated; use solute_pos='packmol'.",
            DeprecationWarning,
            stacklevel=3,
        )
        solute_pos = "packmol"

    # PACKMOL-based placement: full box, left half, or right half
    if solute_pos in (None, "packmol"):
        bounds = box
    elif solute_pos == "left":
        bounds = [box[0], box[1], box[2], box[3], box[4], volume[2] / 2]
    elif solute_pos == "right":
        bounds = [box[0], box[1], volume[2] / 2, box[3], box[4], box[5]]
    else:
        raise ValueError(
            f"Unknown solute_pos {solute_pos!r}. "
            "Choose from: None, 'packmol', 'left', 'right', 'center'."
        )

    logger.debug("  >> solute placement: %s", solute_pos or "packmol (full box)")
    return [
        (
            sp.to_universe(),
            nsolute if isinstance(nsolute, int) else nsolute[cc],
            "box",
            bounds,
        )
        for cc, sp in enumerate(solute)
    ]


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def make_solvent_box(
    species: List[Any],
    solvent: Optional[Any],
    solute: Optional[List[Any]],
    volume: List[float],
    density: Optional[float],
    nsolute: Optional[Union[int, List[int]]],
    concentration: Optional[float],
    conmodel: Optional[Dict[int, Tuple[List[float], List[float]]]],
    solute_pos: Optional[str],
    nsolvent=None,
    tolerance: float = 2.0,
    ratio: Optional[List[float]] = None,
) -> Optional[mda.Universe]:
    """
    Build a solvent box with optional dissolved species.

    Parameters
    ----------
    species : list
        All registered species (used for topology look-up after PACKMOL).
    solvent : Specie, list of Specie, or None
        Solvent molecule(s). Pass a list for mixed-solvent boxes.
    solute : list of Specie or None
        Species to dissolve (ions, neutral molecules, …).
    volume : list of float
        Box dimensions [x, y, z] in Å.
    density : float or None
        Solvent density in g/cm³.  Ignored when *nsolvent* is given.
    nsolute : int, list, or None
        Number of each solute species.
    concentration : float or None
        Solute concentration in Molar (alternative to *nsolute*).
    conmodel : dict or None
        Spatially varying concentration model.
    solute_pos : str or None
        Solute placement strategy: ``"box"``, ``"center"``, ``"left"``, or
        ``None`` for random fixed placement.
    nsolvent : int, list, or None
        Number of solvent molecules.  List length must match *solvent* list.
    tolerance : float
        PACKMOL minimum inter-molecule distance (Å).  Default 2.0.
    ratio : list of float or None
        Molar mixing ratio for multi-solvent boxes.  Must be used with
        either *density* or *nsolvent* (int).

    Returns
    -------
    MDAnalysis.Universe or None
    """

    # Normalise solvent to a list (backward-compatible).
    if solvent is None:
        solvents = []
    elif isinstance(solvent, (list, tuple)):
        solvents = list(solvent)
    else:
        solvents = [solvent]

    _validate_solvent_box_parameters(nsolute, concentration, conmodel, solute,
                                     solvents, density, nsolvent, ratio)

    # convert concentration to number of solute molecules
    if concentration is not None:
        nsolute = int(concentration * np.prod(volume) * units.mol / ((units.m / 10) ** 3))

    instructions = []

    # solute
    if (conmodel is not None) or (nsolute is not None and solute is not None):
        solute_instr = populate_solutes(solute, nsolute, volume, solute_pos=solute_pos,
                                        conmodel=conmodel)
        instructions.extend(solute_instr)

    # solvents
    if solvents:
        solvent_volume = 1e-24 * np.prod(volume)

        if ratio is not None:
            ratio_arr = np.array(ratio, dtype=float)
            if nsolvent is not None:
                counts = [max(1, int(r / ratio_arr.sum() * nsolvent)) for r in ratio_arr]
            else:
                masses = np.array([s.atoms.get_masses().sum() for s in solvents])
                mass_per_unit = float(np.dot(ratio_arr, masses))
                n_units = units.mol * density * solvent_volume / mass_per_unit
                counts = [max(1, int(r * n_units)) for r in ratio_arr]
        elif isinstance(nsolvent, (list, tuple)):
            counts = list(nsolvent)
        elif nsolvent is not None:
            counts = [int(nsolvent)]
        else:
            mass = solvents[0].atoms.get_masses().sum()
            counts = [int(units.mol * density * (1.0 / mass) * solvent_volume)]

        for sp, n in zip(solvents, counts):
            if n > 0:
                instructions.append([sp.to_universe(), n, "box"])

    universe = populate_box(volume, instructions, tolerance=tolerance)

    if universe is None:
        return None

    species_dict = {specie.residues.resnames[0]: specie for specie in species}

    alist = []
    for res in universe.residues:
        resname = res.resname
        if resname in species_dict:
            nmol = species_dict[resname].copy()
            nmol.atoms.positions = res.atoms.positions
            alist.append(nmol.atoms)

    solution = mda.Merge(*alist)
    solution.dimensions = volume + [90, 90, 90]

    return solution
