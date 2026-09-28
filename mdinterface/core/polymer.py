#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Polymer class: a Specie built from one or more repeating monomer units.

Handles chain assembly, LigParGen-based topology refinement at junction
points while preserving molecular connectivity and formal-charge sites.
"""

# repo stuff
from .specie import Specie
from mdinterface.externals import run_ligpargen
from mdinterface.build.polymerize import build_polymer
from mdinterface.build.snippets import make_snippet, remap_snippet_topology

# import random
import numpy as np
import logging

logger = logging.getLogger(__name__)

#%%

class Polymer(Specie):
    """
    A polymer chain built from one or more repeating monomer Specie objects.

    Inherits from :class:`~mdinterface.core.specie.Specie`.  The monomers are
    assembled into a linear chain with :func:`build_polymer`, then the parent
    class handles topology and force-field setup.  Optionally, LigParGen is
    called at each junction point to refine atom types and charges.

    Parameters
    ----------
    monomers : Specie or list of Specie
        Monomer unit(s) to polymerize. A list defines a co-polymer sequence.
        Parameters are inherited from Specie inputs when no explicit topology
        overrides are supplied, with distinct labels for each repeat.
        Each monomer's ASE ``Atoms`` must carry a ``polymerize`` array marking
        the leaving atom (e.g. H or F) at each chain end: value ``1`` = head,
        ``2`` = tail.  Those atoms are deleted during assembly and the bond
        forms between the heavy atoms they were attached to.
    nrep : int, optional
        Number of monomer repetitions in the chain.
    sequence : list of int, optional
        Explicit monomer sequence indices when *monomers* is a list (e.g.
        ``[0, 1, 0, 1]`` alternates two monomers).
    refine_polymer : bool, default False
        If True, run LigParGen at every junction point to obtain accurate
        OPLS-AA parameters for the chain interior.  Requires LigParGen.
    charge_correction : {"none", "uniform"}, default "none"
        Optional uniform correction to the formal charge after refinement.
    cap_element : str, default "H"
        Element symbol used to cap dangling bonds at chain termini during
        LigParGen snippet calculations.

    Notes
    -----
    Parameters not listed above (``charges``, ``atom_types``, ``lj``,
    ``cutoff``, ``name``, ``lammps_data``, ``fix_missing``, ``chg_scaling``,
    ``pbc``, ``ligpargen``, ``tot_charge``) are forwarded to
    :class:`~mdinterface.core.specie.Specie`.

    Examples
    --------
    ::

        from mdinterface import Polymer
        chain = Polymer(monomer, nrep=10)
        chain = Polymer([monomer_A, monomer_B], sequence=[0,1,0,1,0,1])
    """

    def __init__(self, monomers=None, charges=None, atom_types=None, bonds=None,
                 angles=None, dihedrals=None, impropers=None, lj={}, cutoff=1.0,
                 name=None, lammps_data=None, fix_missing=False, chg_scaling=1.0,
                 pbc=False, ligpargen=False, tot_charge=None, nrep=None,
                 sequence=None, refine_polymer=False, charge_correction="none", cap_element="H"):

        # initialize polymer stuff
        self._sequence = sequence
        
        templates = monomers if isinstance(monomers, list) else [monomers]
        inherit = all(value is None for value in (atom_types, bonds, angles, dihedrals, impropers))
        if inherit:
            order = sequence if sequence is not None else [0] * nrep
            templates = [templates[index] for index in order]
            sequence = list(range(len(templates)))
        prepared, inherited = self._prepare_monomers(templates, inherit)
        if inherit:
            bonds, angles, dihedrals, impropers = [inherited[key] for key in ("bonds", "angles", "dihedrals", "impropers")]
        polymer = build_polymer(prepared, sequence=sequence, nrep=nrep)
        
        # Initialize the parent class with polymerized monomers
        super().__init__(atoms=polymer, charges=charges, atom_types=atom_types,
                         bonds=bonds, angles=angles, dihedrals=dihedrals,
                         impropers=impropers, lj=lj, cutoff=cutoff, name=name,
                         lammps_data=lammps_data, fix_missing=fix_missing,
                         chg_scaling=chg_scaling, pbc=pbc, ligpargen=ligpargen,
                         tot_charge=tot_charge)
        
        if refine_polymer:
            self.refine_junctions(charge_correction=charge_correction, cap_element=cap_element)
        
        return
    
    @staticmethod
    def _prepare_monomers(monomers, inherit):
        prepared = []
        inherited = {key: [] for key in ("bonds", "angles", "dihedrals", "impropers")}
        for index, monomer in enumerate(monomers):
            if not isinstance(monomer, Specie):
                prepared.append(monomer)
                continue
            atoms = monomer.atoms.copy()
            types = []
            for atom_index, sid in enumerate(monomer._sids):
                atom_type = monomer._stype[monomer._smap[sid]].copy()
                atom_type.set_label(f"{atom_type.symbol}_M{index}_{atom_index}" if inherit else str(sid))
                types.append(atom_type)
            atoms.set_array("stype", np.array(types, dtype=object))
            if inherit:
                for key, type_key in (("bonds", "_btype"), ("angles", "_atype"),
                                      ("dihedrals", "_dtype"), ("impropers", "_itype")):
                    interactions, type_indices = getattr(monomer, key)
                    for indices, type_index in zip(interactions, type_indices):
                        parameter = getattr(monomer, type_key)[type_index].copy()
                        parameter.update(**{f"a{i + 1}": types[atom_index].label for i, atom_index in enumerate(indices)})
                        inherited[key].append(parameter)
            prepared.append(atoms)
        return prepared, inherited

    @property
    def junction_bonds(self):
        """list of tuple of int: Inter-monomer bonds in current ASE atom indices."""
        monomers = self.atoms.arrays["mon_id"]
        return [(a, b) for a, b in self.graph.edges if monomers[a] != monomers[b]]

    def _get_start_end(self):
        
        str_idx = np.argwhere(self.atoms.arrays["is_connected"] == -1).flatten()
        end_idx = np.argwhere(self.atoms.arrays["is_connected"] == -2).flatten()
        
        return np.concatenate([str_idx, end_idx])
    

    def _update_connection(self, center, partner, Nmax, charges, ending="H", snippet_data=None):

        # make a lil snippet
        snippet, snippet_idxs = snippet_data if snippet_data is not None else make_snippet(self, center, Nmax, ending=ending)
        
        # get local indexes (within dihedral from center)
        ldxs = list(set(np.concatenate(self.find_relevant_distances(4, centers=center))))
        
        # find mapping between indexes
        mapping = [np.argwhere(snippet_idxs == ll)[0][0] for ll in ldxs]
    
        sn_charge = int(snippet.arrays["nominal_charge"].sum())
        sn_atoms, sn_atypes, sn_bonds, sn_angles, sn_dihedrals, sn_impropers = run_ligpargen(
            snippet, charge=sn_charge, is_snippet=True,
        )
        new_charges = sn_atoms.get_initial_charges()
        if len(sn_atoms) != len(snippet) or not np.array_equal(sn_atoms.numbers, snippet.numbers):
            raise ValueError("Junction parameterization changed the snippet atom ordering.")
        if not np.isfinite(new_charges).all():
            raise ValueError("Junction parameterization returned nonfinite partial charges.")

        # update topology of section
        original_idxs = self._sids[snippet_idxs]
        local_idxs = self._sids[ldxs]
        
        # remap topology of local indexes
        new_bonds, new_angles, new_dihedrals, new_impropers = remap_snippet_topology(
            original_idxs, sn_atoms, sn_atypes, sn_bonds, sn_angles, sn_dihedrals,
            sn_impropers, local_idxs)
        
        charges[ldxs] = new_charges[mapping]

        # Correct OPLS types for the two junction atoms (they had H caps
        # substituting their bonded partner during monomer LigParGen runs).
        self._update_junction_lj_types([center, partner], snippet_idxs, sn_atypes)

        self._add_to_topology(bonds=new_bonds, angles=new_angles,
                              dihedrals=new_dihedrals, impropers=new_impropers)
        
        return
    
    def refine_junctions(self, snippet_radius=12, charge_correction="none", cap_element="H"):
        """Refine junction parameters with LigParGen and optionally correct total charge.

        Parameters
        ----------
        snippet_radius : int, default 12
            Graph distance in bonds used to select each junction snippet.
        charge_correction : {"none", "uniform"}, default "none"
            Whether to distribute the residual uniformly over all atoms.
        cap_element : str, default "H"
            Neutral capping element for cut single bonds.

        Returns
        -------
        dict
            Charge audit in elementary-charge units: ``formal_charge``,
            ``initial_charge``, ``refined_charge`` before correction,
            ``residual`` (refined minus formal), ``correction_per_atom`` and
            ``final_charge``. ``junctions`` is the number parameterized.

        Raises
        ------
        ValueError
            If the chain is disconnected or a snippet cannot preserve its
            chemical structure, or an option is invalid.

        Notes
        -----
        All changes are staged on an independent topology. If a junction
        calculation or validation fails, the original chain is unchanged.
        """
        import networkx as nx

        if not nx.is_connected(self.graph):
            raise ValueError("Polymer refinement requires a connected chemical graph.")
        if not isinstance(snippet_radius, (int, np.integer)) or snippet_radius < 4:
            raise ValueError("snippet_radius must be an integer of at least 4.")
        if charge_correction not in {"none", "uniform"}:
            raise ValueError("charge_correction must be 'none' or 'uniform'.")
        if cap_element not in {"H", "F", "Cl", "Br", "I"}:
            raise ValueError("cap_element must be a neutral monovalent element.")
        self.to_rdkit()
        staged, attributes = self._parameterization_copy()
        report = staged._refine_junctions(snippet_radius, charge_correction, cap_element)
        self._apply_parameterization(staged, attributes)
        return report

    def _refine_junctions(self, snippet_radius, charge_correction, cap_element):

        # Clean topology first to remove any invalid interactions from polymerization
        self._cleanup_topology()

        # get charges and connection elements
        charges = self.charges
        initial_charge = float(charges.sum())
        pairs = self.junction_bonds


        snippets = [make_snippet(self, int(pair[1]), snippet_radius, ending=cap_element) for pair in pairs]
        for pair, (snippet, _) in zip(pairs, snippets):
            if len(snippet) > 200:
                raise ValueError(f"Junction {pair} has {len(snippet)} atoms after capping and chemical expansion; LigParGen accepts at most 200. Reduce snippet_radius or use another parameterization method.")
        for pair, snippet_data in zip(pairs, snippets):
            ci, cj = int(pair[0]), int(pair[1])
            self._update_connection(cj, ci, snippet_radius, charges, ending=cap_element, snippet_data=snippet_data)

        if not np.isfinite(charges).all():
            raise ValueError("Junction refinement returned nonfinite partial charges.")
        target = int(self.atoms.arrays["nominal_charge"].sum())
        refined_charge = float(charges.sum())
        residual = refined_charge - target
        correction = 0.0
        logger.info("Polymer charge after junction refinement: %.8f e; formal: %d e; residual: %+.8f e",
                    charges.sum(), target, residual)
        if charge_correction == "uniform":
            correction = -residual / len(charges)
            logger.info("Applying uniform partial-charge correction: %+.8f e per atom", correction)
            charges += correction

        self.atoms.set_initial_charges(charges)
        return {
            "formal_charge": target,
            "initial_charge": initial_charge,
            "refined_charge": refined_charge,
            "residual": residual,
            "correction_per_atom": correction,
            "final_charge": float(charges.sum()),
            "junctions": len(pairs),
        }
