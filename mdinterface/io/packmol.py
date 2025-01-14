#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 13:57:02 2023

@author: roncoroni
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
