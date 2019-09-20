# -*- coding: utf-8 -*-

"""
Distance Matrix.

Computation of the distance matrix from a gene tree.

Methods in this module:
    - distance_matrix
"""

import itertools
import numpy as np


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


def distance_matrix(tree, leaves=None, leaf_index=None):
    """Compute the distance matrix from a gene tree.
    
    Additionally a list of leaves and a dictionary containing
    the indices as values is computed (if not supplied)
    and returned.
    
    Keyword arguments:
        leaves -- set of leaves.
        leaf_index -- dictionary having the leaves as keys
                      and matrix indices as values.
    """
    distance_dict, _ = tree.distances_from_root()
    
    if not leaves:
        tree.supply_leaves()
        color_dict = {}
        for leaf in tree.root.leaves:
            if leaf.color not in color_dict:
                color_dict[leaf.color] = []
            color_dict[leaf.color].append(leaf)
        leaves = []                                                 # color-sorted leaf list
        for color, leaf_list in color_dict.items():
            for leaf in leaf_list:
                leaves.append(leaf)
    
    if not leaf_index:                                              # maps leaf to index in matrix
        leaf_index = {leaf: i for i, leaf in enumerate(leaves)}
    
    D = np.zeros((len(leaves),len(leaves)), dtype=np.float)
    
    for v in tree.preorder():
        if v.children:
            for c1, c2 in itertools.combinations(v.children, 2):
                for x in c1.leaves:
                    x_index = leaf_index[x]
                    x_dist = distance_dict[x] - distance_dict[v]
                    for y in c2.leaves:
                        y_index = leaf_index[y]
                        y_dist = distance_dict[y] - distance_dict[v]
                        D[x_index, y_index] = x_dist + y_dist
                        D[y_index, x_index] = x_dist + y_dist
    
    return leaves, leaf_index, D