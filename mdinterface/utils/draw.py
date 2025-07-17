#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 18:55:16 2025

@author: roncofaber
"""

import matplotlib.pyplot as plt

def draw_bond_markers(ax, mol, node_pos, jmol_colors):
    # Assign marker and color for each bond type
    bond_markers = ['o', 's', 'D', '^', 'v', '*', 'P', 'X']
    bond_colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'cyan', 'magenta']
    bondtype_to_marker = {}
    bondtype_to_color = {}
    for idx, bond in enumerate(mol._btype):
        bondtype_to_marker[idx] = bond_markers[idx % len(bond_markers)]
        bondtype_to_color[idx] = bond_colors[idx % len(bond_colors)]
    
    # Draw markers at bond midpoints
    for (a, b), btype_idx in zip(mol.bonds[0], mol.bonds[1]):
        marker = bondtype_to_marker[btype_idx]
        color = bondtype_to_color[btype_idx]
        x0, y0 = node_pos[a]
        x1, y1 = node_pos[b]
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        ax.scatter(xm, ym, marker=marker, color=color, s=100, zorder=5)
    
    # Legend
    legend_handles = []
    for idx, bond in enumerate(mol._btype):
        marker = bondtype_to_marker[idx]
        color = bondtype_to_color[idx]
        label = f"$k_r={bond.kr}$, $r_0={bond.r0}$"
        h = plt.Line2D([0], [0], marker=marker, color='w', markerfacecolor=color, markersize=10, label=label, linestyle='None')
        legend_handles.append(h)
    ax.legend(handles=legend_handles, loc='best', fontsize=10, framealpha=1.0, edgecolor='k')
    
    return