#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 19:43:33 2025

@author: roncofaber
"""

# repo
from mdinterface.io.read import read_lammps_data_file

# not repo
import logging
import math
import copy
import numpy as np
import os
import ase
import ase.io
import tempfile
import shutil
import subprocess
from ase.data import atomic_numbers, covalent_radii

logger = logging.getLogger(__name__)

#%%

# ---------------------------------------------------------------------------
# Internal helpers for refine_large_specie_topology
# ---------------------------------------------------------------------------

_HETEROATOMS = {"N", "O", "S", "P", "F", "Cl", "Br", "I"}
# Elements that are never acceptable as cut-bond endpoints (true functional
# group centres / halogens -- NOT Si, which is a backbone element in silicones)
_FORBIDDEN = {"N", "P", "F", "Cl", "Br", "I"}


def _candidate_cut_bonds(specie, n_needed=1):
    """Return (i, j) edges suitable for splitting the molecule.

    Uses a tiered strategy: strictest criteria first, progressively relaxed
    until at least *n_needed* candidates are found.

    Tier 1 -- preferred, organic backbones
        Both atoms are C, non-ring, non-terminal, neither adjacent to a
        heteroatom or ring atom.
    Tier 2 -- inorganic / mixed backbones (e.g. silicones)
        Both atoms are non-forbidden, non-ring, non-terminal, non-H, neither
        adjacent to a forbidden atom or ring atom.  Allows Si–C, Si–O, Si–Si.
    Tier 3 -- fallback
        Same as Tier 2 but drops the adjacency restriction entirely.
    """
    rings = specie._find_rings()
    ring_atoms = set()
    for ring in rings:
        ring_atoms.update(ring)

    g = specie.graph
    elements = {n: g.nodes[n]["element"] for n in g.nodes}

    def _nbrs_match(node, avoid, exclude):
        return any(elements[n] in avoid or n in avoid
                   for n in g.neighbors(node) if n != exclude)

    def _base_ok(i, j):
        """Conditions common to all tiers."""
        return (
            g.degree(i) > 1 and g.degree(j) > 1
            and elements[i] != "H" and elements[j] != "H"
            and i not in ring_atoms and j not in ring_atoms
        )

    tiers = [
        # Tier 1: C-C only, no adjacent heteroatoms or ring atoms
        lambda i, j: (
            _base_ok(i, j)
            and elements[i] == "C" and elements[j] == "C"
            and not _nbrs_match(i, _HETEROATOMS | ring_atoms, j)
            and not _nbrs_match(j, _HETEROATOMS | ring_atoms, i)
        ),
        # Tier 2: any non-forbidden bond, no adjacent forbidden/ring atoms
        lambda i, j: (
            _base_ok(i, j)
            and elements[i] not in _FORBIDDEN
            and elements[j] not in _FORBIDDEN
            and not _nbrs_match(i, _FORBIDDEN | ring_atoms, j)
            and not _nbrs_match(j, _FORBIDDEN | ring_atoms, i)
        ),
        # Tier 3: any non-forbidden, non-ring, non-terminal bond
        lambda i, j: (
            _base_ok(i, j)
            and elements[i] not in _FORBIDDEN
            and elements[j] not in _FORBIDDEN
        ),
    ]

    import networkx as nx

    # Imbalance threshold: best achievable split must put at least 25% of
    # atoms on the minority side.  Tiers that can't meet this fall through.
    n_atoms    = len(specie.atoms)
    max_diff   = n_atoms * 0.5   # 25% minority → diff ≤ 50% of total

    for tier_num, check in enumerate(tiers):
        candidates = [(i, j) for i, j in g.edges() if check(i, j)]
        if len(candidates) < n_needed:
            continue

        # Check whether any candidate gives a reasonably balanced split
        best_diff = float("inf")
        for ci, cj in candidates:
            sg = g.copy()
            sg.remove_edge(ci, cj)
            comps = list(nx.connected_components(sg))
            if len(comps) == 2:
                best_diff = min(best_diff, abs(len(comps[0]) - len(comps[1])))

        if best_diff <= max_diff:
            if tier_num > 0:
                logger.info(
                    "  >> using tier-%d cut criteria "
                    "(tier-1 bonds cannot give a balanced split)",
                    tier_num + 1)
            return candidates

    return []


def _balanced_cuts(specie, n_cuts, candidates):
    """Select *n_cuts* bonds from *candidates* to partition the molecule into
    *n_cuts + 1* roughly equal segments using recursive halving.

    Returns
    -------
    cut_edges : list of (i, j)
    segments  : list of sets of atom indices
    """
    import networkx as nx

    available = set(map(frozenset, candidates))

    def _best_cut(node_set):
        subg = specie.graph.subgraph(node_set)
        best_edge, best_diff = None, float("inf")
        for i, j in subg.edges():
            if frozenset([i, j]) not in available:
                continue
            sg = subg.copy()
            sg.remove_edge(i, j)
            comps = list(nx.connected_components(sg))
            if len(comps) != 2:
                continue
            diff = abs(len(comps[0]) - len(comps[1]))
            if diff < best_diff:
                best_diff, best_edge = diff, (i, j)
        return best_edge

    def _recurse(node_set, remaining):
        if remaining == 0:
            return [], [node_set]
        edge = _best_cut(node_set)
        if edge is None:
            logger.warning(
                "No valid cut bond found in a segment of %d atoms -- "
                "that segment may exceed 200 atoms.", len(node_set))
            return [], [node_set]
        i, j = edge
        available.discard(frozenset([i, j]))
        subg = specie.graph.subgraph(node_set).copy()
        subg.remove_edge(i, j)
        comp_a, comp_b = sorted(nx.connected_components(subg), key=len)
        cuts_a = round(remaining * len(comp_a) / len(node_set))
        cuts_b = remaining - 1 - cuts_a
        edges_a, segs_a = _recurse(comp_a, cuts_a)
        edges_b, segs_b = _recurse(comp_b, cuts_b)
        return [(i, j)] + edges_a + edges_b, segs_a + segs_b

    return _recurse(set(specie.graph.nodes()), n_cuts)


def _make_capped_segment(specie, seg_indices, cut_edges, ending="H"):
    """Build a capped ASE Atoms for *seg_indices*, adding one *ending* atom
    per bond that crosses the segment boundary.

    Returns
    -------
    capped : ase.Atoms  (real atoms first, then caps)
    n_real : int        (number of non-cap atoms)
    """
    seg_set = set(seg_indices)
    atoms   = specie.atoms

    cap_positions, cap_nominal = [], []
    for ci, cj in cut_edges:
        if ci in seg_set and cj not in seg_set:
            inner, outer = ci, cj
        elif cj in seg_set and ci not in seg_set:
            inner, outer = cj, ci
        else:
            continue

        pos_i = atoms.positions[inner]
        pos_o = atoms.positions[outer]
        d = pos_o - pos_i
        d /= np.linalg.norm(d)
        r_inner = covalent_radii[atomic_numbers[atoms[inner].symbol]]
        r_cap   = covalent_radii[atomic_numbers[ending]]
        cap_positions.append(pos_i + d * (r_inner + r_cap))
        cap_nominal.append(0)

    seg_atoms = atoms[list(seg_indices)].copy()
    n_real = len(seg_atoms)

    if cap_positions:
        caps   = ase.Atoms(symbols=[ending] * len(cap_positions),
                           positions=cap_positions)
        capped = seg_atoms + caps
    else:
        capped = seg_atoms

    if "nominal_charge" in atoms.arrays:
        nc = atoms.arrays["nominal_charge"][list(seg_indices)].tolist() + cap_nominal
    else:
        nc = [0] * len(capped)
    capped.set_array("nominal_charge", np.array(nc, dtype=int))

    return capped, n_real

def run_ligpargen(system, charge=None, is_snippet=False):
    """
    Runs the ligpargen command for the given xyz file.

    Parameters:
    system (ase.Atoms): The atoms system to be processed.

    Returns:
    tuple: Containing system, atoms, bonds, angles, dihedrals, impropers.
    """

    if "BOSSdir" not in os.environ:
        mdint = os.environ["MDINT_CONFIG_DIR"]
        logger.warning(
            "BOSSdir not set. Please either:\n"
            "  os.environ['BOSSdir'] = '/path/to/your/boss'\n"
            "  or add 'BOSSdir = /path/to/boss' in [settings] in %s/config.ini",
            mdint,
        )

    # all ligpargen files go in a temp dir; kept on failure for inspection
    tmpdir   = tempfile.mkdtemp(prefix="ligpargen_")
    mol_name = os.path.basename(tmpdir)
    xyz_file = os.path.join(tmpdir, f"{mol_name}.xyz")

    ase.io.write(xyz_file, system)

    # use relative filenames and cwd=tmpdir -- ligpargen does not accept
    # absolute paths for -i
    ligpargen_command = ["ligpargen", "-i", f"{mol_name}.xyz", "-p", tmpdir,
                         "-debug", "-o", "0", "-cgen", "CM1A"]
    if charge is not None:
        ligpargen_command.extend(["-c", str(charge)])

    try:
        result = subprocess.run(ligpargen_command, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                cwd=tmpdir)
        logger.debug("ligpargen completed successfully")
        logger.debug("ligpargen stdout:\n%s", result.stdout.decode())

    except subprocess.CalledProcessError as e:
        # write stdout/stderr into the temp dir for inspection, then keep it
        error_log = os.path.join(tmpdir, "error_log.txt")
        with open(error_log, "w") as fh:
            fh.write("STDOUT:\n" + e.stdout.decode() + "\n")
            fh.write("STDERR:\n" + e.stderr.decode() + "\n")
        logger.error("ligpargen failed; temp files kept at: %s", tmpdir)
        logger.debug("ligpargen stderr:\n%s", e.stderr.decode())
        raise

    # read result
    system, atoms, bonds, angles, dihedrals, impropers = read_lammps_data_file(
        os.path.join(tmpdir, f"{mol_name}.lammps.lmp"), is_snippet=is_snippet)

    # success -- clean up
    shutil.rmtree(tmpdir, ignore_errors=True)

    return system, atoms, bonds, angles, dihedrals, impropers


def refine_large_specie_topology(specie, Nmax=12, ending="H", offset=True,
                                 segment_size=200):
    """
    Assign OPLS-AA force-field parameters to a ``Specie`` via LigParGen using
    a segment-and-junction strategy.

    The molecule is split into segments of at most *segment_size* atoms along
    clean backbone bonds (avoiding rings, heteroatoms, and their immediate
    neighbours).  Each segment is capped with *ending* atoms and passed to
    LigParGen independently.  A local snippet centred on each cut bond is then
    refined to correct parameters for atoms adjacent to a cap.

    Parameters
    ----------
    specie : Specie
        The molecule to parametrise.  Modified **in place**.
    Nmax : int
        Neighbourhood radius (in bonds) used for junction snippet creation.
        Default 12.
    ending : str
        Element used to cap dangling bonds at segment boundaries.  Default
        ``"H"``.
    offset : bool
        If ``True``, redistribute any charge rounding error uniformly so the
        total partial charge matches ``nominal_charge.sum()``.
    segment_size : int
        Maximum number of atoms per segment.  Default 200.  Lower values
        force more splits and are useful for testing.

    Raises
    ------
    ValueError
        If ``nominal_charge`` is not set on ``specie.atoms``.
    RuntimeError
        If not enough valid cut bonds can be found.
    """
    from mdinterface.build.snippets import make_snippet, remap_snippet_topology

    natoms = len(specie.atoms)

    if "nominal_charge" not in specie.atoms.arrays:
        raise ValueError(
            "nominal_charge array not found on specie.atoms.  "
            "Set it before calling refine_large_specie_topology() -- "
            "use 0 for uncharged atoms and the formal integer charge for "
            "charged ones (same convention as the Polymer class)."
        )

    n_segments = math.ceil(natoms / segment_size)
    n_cuts     = n_segments - 1
    logger.info("Refining topology for %d-atom molecule -- %d segment(s) via LigParGen",
                natoms, n_segments)

    # ------------------------------------------------------------------
    # 1. Find cut bonds
    # ------------------------------------------------------------------
    candidates = _candidate_cut_bonds(specie, n_needed=n_cuts)
    if len(candidates) < n_cuts:
        raise RuntimeError(
            f"Only {len(candidates)} valid cut bond(s) found but "
            f"{n_cuts} cut(s) are needed.  Consider relaxing the molecule "
            f"structure or cutting manually."
        )

    cut_edges, segments = _balanced_cuts(specie, n_cuts, candidates)
    for ii, seg in enumerate(segments):
        logger.info("  >> segment %d: %d atoms", ii, len(seg))
    logger.info("  >> cut bond(s): %s",
                ", ".join(f"{int(a)} -- {int(b)}" for a, b in cut_edges) or "none")

    # ------------------------------------------------------------------
    # 2. Run LigParGen on each segment, accumulate topology
    # ------------------------------------------------------------------
    charges            = specie.charges.copy()
    atom_types_ordered = [None] * natoms  # Atom-type objects in original order
    all_bonds, all_angles, all_dihs, all_imps = [], [], [], []

    for seg_idx, seg_set in enumerate(segments):
        seg_indices = sorted(seg_set)
        logger.info("  >> ligpargen on segment %d (%d atoms)...",
                    seg_idx, len(seg_indices))

        capped, n_real = _make_capped_segment(
            specie, seg_indices, cut_edges, ending)
        sn_charge = int(capped.arrays["nominal_charge"].sum())

        sn_sys, sn_atypes, sn_bonds, sn_angles, sn_dihs, sn_imps = \
            run_ligpargen(capped, charge=sn_charge, is_snippet=True)

        # Remap: real atoms -> original sids; caps -> placeholders
        real_sids     = list(specie._sids[seg_indices])
        cap_sids      = [f"__cap_{seg_idx}_{c}__"
                         for c in range(len(sn_atypes) - n_real)]
        original_idxs = np.array(real_sids + cap_sids)
        local_idxs    = np.array(real_sids)

        b, a, d, i = remap_snippet_topology(
            original_idxs, sn_sys, sn_atypes,
            sn_bonds, sn_angles, sn_dihs, sn_imps,
            local_idxs)
        all_bonds  += b
        all_angles += a
        all_dihs   += d
        all_imps   += i

        seg_charges = sn_sys.get_initial_charges()
        for pos, orig_idx in enumerate(seg_indices):
            charges[orig_idx]            = seg_charges[pos]
            atom_types_ordered[orig_idx] = sn_atypes[pos]

    # ------------------------------------------------------------------
    # 3. Rebuild topology from assembled segment results
    # ------------------------------------------------------------------
    specie._setup_topology(atom_types_ordered,
                           all_bonds, all_angles, all_dihs, all_imps)
    specie.atoms.set_initial_charges(charges)

    # ------------------------------------------------------------------
    # 4. Refine each junction with a local snippet run
    # ------------------------------------------------------------------
    for (ci, cj) in cut_edges:
        # One snippet per cut, centred on ci (it is bonded to cj so the
        # snippet naturally spans both sides of the cut)
        center = ci
        logger.info("  >> refining junction at atoms %d -- %d...", ci, cj)

        snippet, snippet_idxs = make_snippet(specie, center, Nmax,
                                             ending=ending)

        ldxs    = list(set(np.concatenate(
            specie.find_relevant_distances(4, centers=center))))
        mapping = [int(np.argwhere(snippet_idxs == ll)[0][0]) for ll in ldxs]

        if "nominal_charge" not in snippet.arrays:
            raise ValueError(
                f"nominal_charge missing in junction snippet at atom {center}.")
        sn_charge = int(snippet.arrays["nominal_charge"].sum())

        sn_sys, sn_atypes, sn_bonds, sn_angles, sn_dihs, sn_imps = \
            run_ligpargen(snippet, charge=sn_charge, is_snippet=True)

        original_idxs = specie._sids[snippet_idxs]
        local_idxs    = specie._sids[ldxs]

        new_bonds, new_angles, new_dihs, new_imps = remap_snippet_topology(
            original_idxs, sn_sys, sn_atypes,
            sn_bonds, sn_angles, sn_dihs, sn_imps,
            local_idxs)

        new_charges    = sn_sys.get_initial_charges()
        charges[ldxs]  = new_charges[mapping]

        # Correct OPLS types for the two cut-bond atoms (they had H caps in
        # their respective segments so their types may be wrong).
        specie._update_junction_lj_types([ci, cj], snippet_idxs, sn_atypes)

        specie._add_to_topology(bonds=new_bonds, angles=new_angles,
                                dihedrals=new_dihs, impropers=new_imps)

    # ------------------------------------------------------------------
    # 5. Final charge assignment (+ optional offset correction)
    # ------------------------------------------------------------------
    if offset:
        target   = int(specie.atoms.arrays["nominal_charge"].sum())
        charges -= (charges.sum() - target) / natoms

    specie.atoms.set_initial_charges(charges)
    logger.info("  >> done.  Total charge: %.4f  (target %d)",
                charges.sum(),
                int(specie.atoms.arrays["nominal_charge"].sum()))
