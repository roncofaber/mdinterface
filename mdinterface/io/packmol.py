#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PACKMOL input-file string templates.

Provides the ``header``, ``box_place``, and ``fix_place`` format strings used
by :func:`~mdinterface.build.box.populate_box` to generate PACKMOL ``.in``
files.
"""


header = """\
tolerance {}
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
