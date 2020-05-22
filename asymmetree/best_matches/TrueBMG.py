# -*- coding: utf-8 -*-

"""
Orthology graph, cBMG, cRBMG.

This module provides classes concerning colored graphs and colored
(phylogentic) trees, including a cBMG(PhyloTree)-function.
"""

import itertools

import networkx as nx


__author__ = "David Schaller"


# --------------------------------------------------------------------------
#                 True Orthology Graph (TOG), cBMG & cRBMG
#
#                         (given a known tree)
# --------------------------------------------------------------------------

def orthology_from_tree(tree):
    """Constructs the true orthology graph from an event-labeled tree."""
    
    tree.supply_leaves()                                # assign list of leaves to each node
    G = nx.Graph()
    
    for v in tree.root.leaves:
        G.add_node(v.ID, label=v.label, color=v.color)
    
    for node in tree.preorder():
        if node.label == 'S':
            for child1, child2 in itertools.combinations(node.children, 2):
                for u in child1.leaves:
                    for v in child2.leaves:
                        G.add_edge(u.ID, v.ID)
    return G


def BMG_from_tree(tree, supply_RBMG=False):
    """Create an n-colored BMG (and optionally RBMG) from a given tree."""
    
    tree.supply_leaves()                                # assign list of leaves to each node
    BMG = nx.DiGraph()
    colors = set()
    
    for v in tree.root.leaves:
        colors.add(v.color)
        BMG.add_node(v.ID, label=v.label, color=v.color)
    
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
    
    if not supply_RBMG:
        return BMG
    else:
        return BMG, RBMG_from_BMG(BMG)


def RBMG_from_BMG(BMG):
    
    RBMG = nx.Graph()
    RBMG.add_nodes_from(BMG.nodes(data=True))
    
    for x, neighbors in BMG.adjacency():
        for y in neighbors:
            if BMG.has_edge(y,x):
                RBMG.add_edge(x,y)
    return RBMG