#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Optional integrations with third-party tools.

Exposes wrappers for LigParGen (OPLS-AA parameters), Open Babel (charge
models), PySCF/gpu4pyscf (RESP charges), ASE optimisers (structure
relaxation), and FAIRChem (AIMD).  Missing optional dependencies are handled
gracefully at import time.
"""

from .ligpargen import run_ligpargen, refine_large_specie_topology
from .obabel import run_OBChargeModel
from .pyscf import calculate_RESP_charges
from .optimization import relax_structure


def run_aimd(*args, **kwargs):
    """Run AIMD while importing the optional FAIRChem integration on demand."""
    from .aimd import run_aimd as _run_aimd

    return _run_aimd(*args, **kwargs)
