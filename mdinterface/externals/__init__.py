#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
External library integrations for mdinterface.

Interfaces to external libraries including LigParGen, OpenBabel,
and PySCF for force field generation and charge calculations.

Author: Fabrice Roncoroni
Created: 2025-01-15
"""

from .ligpargen import run_ligpargen
from .obabel import run_OBChargeModel
from .pyscf import calculate_RESP_charges
