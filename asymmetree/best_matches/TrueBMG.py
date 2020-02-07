# -*- coding: utf-8 -*-

"""
Orthology graph, cBMG, cRBMG.

This module provides classes concerning colored graphs and colored
(phylogentic) trees, including a cBMG(PhyloTree)-function.

Methods in this module:
    - true_orthology_graph(tree)
    - best_match_graphs(tree)
    - build_cBMG_alternative(tree)
    - RBMG_from_BMG(cBMG)
"""

import itertools

import numpy as np
import networkx as nx


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


# --------------------------------------------------------------------------
#                 True Orthology Graph (TOG), cBMG & cRBMG
#
#                         (given a known tree)
# --------------------------------------------------------------------------

def true_orthology_graph(tree):
    """Constructs the true orthology graph from an event-labeled tree."""
    
    tree.supply_leaves()                                # assign list of leaves to each node
    G = nx.Graph()
    
    for v in tree.root.leaves:
        G.add_node(v.ID, label=v.label, color=v.color)
    
    for node in tree.preorder():
        if node.label == "S":
            for child1, child2 in itertools.combinations(node.children, 2):
                for u in child1.leaves:
                    for v in child2.leaves:
                        G.add_edge(u.ID, v.ID)
    return G


def best_match_graphs(tree):
    """Create an n-colored BMG and RBMG from a given tree.
    
    Own algorithm."""
    
    tree.supply_leaves()                                # assign list of leaves to each node
    BMG, RBMG = nx.DiGraph(), nx.Graph()
    colors = set()
    
    for v in tree.root.leaves:
        colors.add(v.color)
        BMG.add_node(v.ID, label=v.label, color=v.color)
        RBMG.add_node(v.ID, label=v.label, color=v.color)
    
    # ---- build BMG ----
    for u in tree.root.leaves:
        remaining = colors - set([u.color])             # colors to which no best match has yet been found
        parent = u.parent                               # start with direct parent of each node
        while remaining and parent:
            colors_here = set()
            for v in parent.leaves:
                if v.color in remaining:                # best match found
                    colors_here.add(v.color)
                    BMG.add_edge(u.ID, v.ID)            # insert edge (u,v)
            remaining -= colors_here                    # update remaining colors
            parent = parent.parent
    
    # ---- build RBMG as symmetric part of the BMG ----
    for x, neighbors in BMG.adjacency():
        for y in neighbors:
            if BMG.has_edge(y,x):
                RBMG.add_edge(x,y)
                
    return BMG, RBMG


def build_cBMG_alternative(tree):
    """Create an n-colored BMG from a given tree.
    
    Algorithm presented in the cobmg paper."""
    
    V = []
    L = []
    colors = set()
    for v in tree.postorder():
        V.append(v)
        if not v.children:
            L.append(v)
            colors.add(v.color)
    
    l_matrix = np.zeros((len(V), len(colors)), dtype=np.int8)
    v_index = {V[i]: i for i in range(len(V))}
    col_index = {color: i for i, color in enumerate(colors)}
    
    G = nx.DiGraph()
    
    for v in L:
        v.leaves = [v]
        l_matrix[v_index[v], col_index[v.color]] = 1
        G.add_node(v.ID, label=v.label, color=v.color)
    
    for v in tree.postorder():
        if not v.children:
            continue                                        # only traverse inner vertices
        
        for u1, u2 in itertools.combinations(v.children, 2):
            for x in u1.leaves:
                for y in u2.leaves:
                    if l_matrix[v_index[u1], col_index[y.color]] == 0:
                        G.add_edge(x.ID, y.ID)                  # edge (x,y)
                    if l_matrix[v_index[u2], col_index[x.color]] == 0:
                        G.add_edge(y.ID, x.ID)                  # edge (y,x)
        v.leaves = []
        v_nr = v_index[v]
        for u in v.children:
            v.leaves.extend(u.leaves)                       # compute leaves L(v) on the fly
            u_nr = v_index[u]
            for r_nr in range(len(colors)):
                if l_matrix[u_nr, r_nr] == 1:
                    l_matrix[v_nr, r_nr] = 1
    return G


def RBMG_from_BMG(cBMG):
    
    G = nx.Graph()
    G.add_nodes_from(cBMG.nodes(data=True))
    
    for x, neighbors in cBMG.adjacency():
        for y in neighbors:
            if cBMG.has_edge(y,x):
                G.add_edge(x,y)
    return G