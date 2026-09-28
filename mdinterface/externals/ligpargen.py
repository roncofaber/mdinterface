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


class LigParGenError(RuntimeError):
    """LigParGen configuration, execution, or output-processing failure."""

    def __init__(self, message, tempdir=None, log_path=None, returncode=None):
        self.tempdir = tempdir
        self.log_path = log_path
        self.returncode = returncode
        details = [message]
        if returncode is not None:
            details.append(f"return code: {returncode}")
        if tempdir is not None:
            details.append(f"temporary files: {tempdir}")
        if log_path is not None:
            details.append(f"log: {log_path}")
        super().__init__("; ".join(details))


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
    import networkx as nx
    bridges = {frozenset(edge) for edge in nx.bridges(specie.graph)}
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
            frozenset((i, j)) in bridges
            and g.edges[i, j].get("bond_order", 1) == 1
            and not g.nodes[i].get("formal_charge", 0)
            and not g.nodes[j].get("formal_charge", 0)
            and g.degree(i) > 1 and g.degree(j) > 1
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
        cuts_a = min(remaining - 1, round((remaining - 1) * len(comp_a) / len(node_set)))
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
    from mdinterface.core.chemistry import capped_molecule

    capped, _ = capped_molecule(specie.atoms, specie.to_rdkit(), seg_indices, ending=ending)
    return capped, len(seg_indices)


def run_ligpargen(system, charge=None, is_snippet=False):
    """Generate OPLS-AA parameters by running LigParGen.

    Parameters
    ----------
    system : ase.Atoms
        Atomic system to parameterize. Stored RDKit chemistry is transferred
        through MOL files; coordinate-only inputs use XYZ and Open Babel.
    charge : int or None, default None
        Total molecular charge. LigParGen detects it when omitted.
    is_snippet : bool, default False
        Whether the system is a capped molecular snippet.

    Returns
    -------
    tuple
        Parameterized system, atom types, bonds, angles, dihedrals, and impropers.

    Raises
    ------
    LigParGenError
        If LigParGen or BOSS is not configured, execution fails, or the output
        is missing or unreadable.
    """

    if len(system) > 200:
        raise ValueError(f"LigParGen accepts at most 200 atoms including caps; received {len(system)}.")

    ligpargen_executable = shutil.which("ligpargen")
    if ligpargen_executable is None:
        raise LigParGenError(
            "LigParGen executable was not found on PATH. Install it in the active "
            "environment with `python -m pip install "
            "\"git+https://github.com/roncofaber/ligpargen.git@ad78036842318f166531be41cfcbc3563d7c5476\"` and verify "
            "the installation with `ligpargen -h`"
        )

    from rdkit import Chem
    from mdinterface.core.chemistry import stored_molecule

    mol = stored_molecule(system)
    if mol is not None:
        formal_charge = Chem.GetFormalCharge(mol)
        if charge is not None and charge != formal_charge:
            raise ValueError("LigParGen charge conflicts with the molecular formal charge.")
        charge = formal_charge
    if mol is None and shutil.which("obabel") is None:
        raise LigParGenError(
            "Open Babel executable `obabel` was not found on PATH. LigParGen "
            "requires it to read mdinterface's XYZ input. Install it with "
            "`conda install -c conda-forge openbabel` and verify the installation "
            "with `obabel -V`"
        )

    if not os.environ.get("BOSSdir"):
        from mdinterface.config import load_config

        load_config()

    if not os.environ.get("BOSSdir"):
        mdint = os.environ.get("MDINT_CONFIG_DIR", "~/.config/mdinterface")
        config_file = os.path.join(os.path.expanduser(mdint), "config.ini")
        raise LigParGenError(
            "BOSSdir is not configured. Export BOSSdir or set it under [settings] "
            f"in {config_file}"
        )

    # all ligpargen files go in a temp dir; kept on failure for inspection
    tmpdir   = tempfile.mkdtemp(prefix="ligpargen_")
    mol_name = os.path.basename(tmpdir)
    extension = "mol" if mol is not None else "xyz"
    input_file = os.path.join(tmpdir, f"{mol_name}.{extension}")
    log_file = os.path.join(tmpdir, "ligpargen.log")

    if mol is None:
        ase.io.write(input_file, system)
    else:
        Chem.MolToMolFile(mol, input_file)

    # use relative filenames and cwd=tmpdir -- ligpargen does not accept
    # absolute paths for -i
    ligpargen_command = [ligpargen_executable, "-i", f"{mol_name}.{extension}", "-p", tmpdir,
                         "-debug", "-o", "0", "-cgen", "CM1A"]
    if charge is not None:
        ligpargen_command.extend(["-c", str(charge)])

    try:
        result = subprocess.run(ligpargen_command, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                cwd=tmpdir, text=True, encoding="utf-8",
                                errors="replace")
        with open(log_file, "w") as fh:
            fh.write("STDOUT:\n" + result.stdout + "\n")
            fh.write("STDERR:\n" + result.stderr + "\n")
        logger.debug("ligpargen completed successfully")
        logger.debug("ligpargen stdout:\n%s", result.stdout)

    except subprocess.CalledProcessError as e:
        with open(log_file, "w") as fh:
            fh.write("STDOUT:\n" + (e.stdout or "") + "\n")
            fh.write("STDERR:\n" + (e.stderr or "") + "\n")
        logger.error("ligpargen failed; temp files kept at: %s", tmpdir)
        logger.debug("ligpargen stderr:\n%s", e.stderr)
        raise LigParGenError(
            "LigParGen exited unsuccessfully",
            tmpdir,
            log_file,
            returncode=e.returncode,
        ) from e
    except OSError as e:
        raise LigParGenError(
            f"LigParGen could not be started: {e}",
            tmpdir,
            log_file,
        ) from e

    output_file = os.path.join(tmpdir, f"{mol_name}.lammps.lmp")
    if not os.path.isfile(output_file):
        diagnostic = (result.stderr or result.stdout).strip()
        message = f"LigParGen did not create the expected output file {output_file}"
        if diagnostic:
            message += f". LigParGen output: {diagnostic[-500:]}"
        raise LigParGenError(
            message,
            tmpdir,
            log_file,
            returncode=result.returncode,
        )

    try:
        original_numbers = system.numbers.copy()
        system, atoms, bonds, angles, dihedrals, impropers = read_lammps_data_file(
            output_file, is_snippet=is_snippet)
        if not np.array_equal(original_numbers, system.numbers):
            raise ValueError("LigParGen output atom ordering differs from the input.")
    except Exception as e:
        raise LigParGenError(
            f"LigParGen output could not be read from {output_file}",
            tmpdir,
            log_file,
            returncode=result.returncode,
        ) from e

    # success -- clean up
    shutil.rmtree(tmpdir, ignore_errors=True)

    return system, atoms, bonds, angles, dihedrals, impropers


def refine_large_specie_topology(specie, snippet_radius=12, cap_element="H",
                                 charge_correction="none", segment_size=200):
    """Assign LigParGen parameters atomically; see ``Specie.parameterize``.

    Parameters
    ----------
    specie : Specie
        Species to parameterize without changing its coordinates.
    snippet_radius : int, default 12
        Graph radius for junction snippets.
    cap_element : str, default "H"
        Neutral monovalent capping element.
    charge_correction : {"none", "uniform"}, default "none"
        Optional uniform correction to the molecular charge.
    segment_size : int, default 200
        Maximum segment size including caps, at most 200.

    Returns
    -------
    dict
        Charge audit. No changes are applied if any calculation fails.
    """
    if not isinstance(segment_size, (int, np.integer)) or not 4 <= segment_size <= 200:
        raise ValueError("segment_size must be an integer between 4 and 200, including caps.")
    if not isinstance(snippet_radius, (int, np.integer)) or snippet_radius < 4:
        raise ValueError("snippet_radius must be an integer of at least 4.")
    if charge_correction not in {"none", "uniform"}:
        raise ValueError("charge_correction must be 'none' or 'uniform'.")
    if cap_element not in {"H", "F", "Cl", "Br", "I"}:
        raise ValueError("cap_element must be a neutral monovalent element.")
    staged, attributes = specie._parameterization_copy()
    initial = float(staged.charges.sum())
    target = staged._resolve_charge(None)
    if len(staged.atoms) <= segment_size:
        result, atom_types, bonds, angles, dihedrals, impropers = run_ligpargen(staged.atoms, charge=target)
        if not np.array_equal(result.numbers, staged.atoms.numbers) or len(atom_types) != len(result):
            raise ValueError("Parameterization changed atom ordering or omitted atom types.")
        staged._setup_topology(atom_types, bonds, angles, dihedrals, impropers)
        staged.atoms.set_initial_charges(result.get_initial_charges())
        junctions = 0
    else:
        junctions = _refine_large_specie_topology(staged, snippet_radius, cap_element, segment_size)
    charges = staged.charges
    if not np.isfinite(charges).all():
        raise ValueError("Parameterization returned nonfinite partial charges.")
    refined = float(charges.sum())
    residual = refined - target
    correction = -residual / len(charges) if charge_correction == "uniform" else 0.0
    staged.atoms.set_initial_charges(charges + correction)
    report = dict(formal_charge=target, initial_charge=initial, refined_charge=refined,
                  residual=residual, correction_per_atom=correction,
                  final_charge=float(staged.charges.sum()), junctions=junctions)
    staged.validate_force_field()
    specie._apply_parameterization(staged, attributes)
    logger.info("Parameterization charge audit: %s", report)
    return report


def _refine_large_specie_topology(specie, snippet_radius, cap_element, segment_size):
    from mdinterface.build.snippets import make_snippet, remap_snippet_topology

    natoms = len(specie.atoms)
    mol = specie.to_rdkit()
    from mdinterface.core.chemistry import graph_from_molecule, store_molecule
    store_molecule(specie.atoms, mol)
    specie._graph = graph_from_molecule(mol)
    capped_segments = []
    for n_cuts in range(max(1, math.ceil(natoms / segment_size) - 1), natoms):
        candidates = _candidate_cut_bonds(specie, n_needed=n_cuts)
        if len(candidates) < n_cuts:
            raise ValueError("Cannot split this molecule into chemically valid capped segments within segment_size; use a larger limit (at most 200) or another parameterization method.")
        cut_edges, segments = _balanced_cuts(specie, n_cuts, candidates)
        capped_segments = [_make_capped_segment(specie, sorted(seg), cut_edges, cap_element) for seg in segments]
        if all(len(capped) <= segment_size for capped, _ in capped_segments):
            break
    snippets = [make_snippet(specie, int(ci), snippet_radius, ending=cap_element) for ci, _ in cut_edges]
    for pair, (snippet, _) in zip(cut_edges, snippets):
        if len(snippet) > 200:
            raise ValueError(f"Junction {pair} expands to {len(snippet)} atoms including caps; LigParGen accepts at most 200. Reduce snippet_radius or use another parameterization method.")

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

        capped, n_real = capped_segments[seg_idx]
        sn_charge = int(capped.arrays["nominal_charge"].sum())

        sn_sys, sn_atypes, sn_bonds, sn_angles, sn_dihs, sn_imps = \
            run_ligpargen(capped, charge=sn_charge, is_snippet=True)

        # Remap: real atoms -> globally unique labels (original atom index
        # avoids collisions across segments); caps -> placeholders
        real_sids     = [f"{sn_atypes[pos].symbol}_{seg_indices[pos]:03d}"
                         for pos in range(n_real)]
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
            atom_types_ordered[orig_idx].set_label(real_sids[pos])

    # ------------------------------------------------------------------
    # 3. Rebuild topology from assembled segment results
    # ------------------------------------------------------------------
    specie._setup_topology(atom_types_ordered,
                           all_bonds, all_angles, all_dihs, all_imps)
    specie.atoms.set_initial_charges(charges)

    # ------------------------------------------------------------------
    # 4. Refine each junction with a local snippet run
    # ------------------------------------------------------------------
    for (ci, cj), snippet_data in zip(cut_edges, snippets):
        # One snippet per cut, centred on ci (it is bonded to cj so the
        # snippet naturally spans both sides of the cut)
        center = ci
        logger.info("  >> refining junction at atoms %d -- %d...", ci, cj)

        snippet, snippet_idxs = snippet_data

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

    specie.atoms.set_initial_charges(charges)
    return len(cut_edges)
