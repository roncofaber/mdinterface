#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Packmol integration utilities.

Template strings and utilities for generating Packmol input files
to pack molecules into simulation boxes.

Author: Fabrice Roncoroni
Created: 2023-11-22
"""


header = """\
tolerance 2.0
output    {}
nloop     100
seed      {}

filetype  pdb
movebadrandom
        """

box_place = """
structure mol_{}.pdb
    number      {}
    inside box  {}
    resnumbers  2
end structure
            """

fix_place = """
structure mol_{}.pdb
    number      1
    center
    fixed {} {} {} 0. 0. 0. 
    resnumbers  2
end structure
"""
