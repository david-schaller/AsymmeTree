# -*- coding: utf-8 -*-

"""
Orthology graph, BMG, RBMG.

This module provides classes concerning colored graphs and colored
(phylogentic) trees, including a BMG(PhyloTree)-function.
"""

import itertools

import networkx as nx

from asymmetree.tools.GraphTools import sort_by_colors


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                 True Orthology Graph (TOG), BMG & RBMG
#
#                         (given a known tree)
# --------------------------------------------------------------------------

def orthology_from_tree(tree):
    """Constructs the true orthology graph from an event-labeled tree."""
    
    tree.supply_leaves()                                # assign list of leaves to each node
    TOG = nx.Graph()
    
    for v in tree.root.leaves:
        TOG.add_node(v.ID, label=v.label, color=v.color)
    
    for node in tree.preorder():
        if node.label == 'S':
            for child1, child2 in itertools.combinations(node.children, 2):
                for u in child1.leaves:
                    for v in child2.leaves:
                        TOG.add_edge(u.ID, v.ID)
                        
    return TOG


def bmg_from_tree(tree, supply_rbmg=False):
    """Construct a BMG (and optionally RBMG) from a given tree."""
    
    tree.supply_leaves()                                # assign list of leaves to each node
    bmg = nx.DiGraph()
    colors = set()
    
    for v in tree.root.leaves:
        colors.add(v.color)
        bmg.add_node(v.ID, label=v.label, color=v.color)
    
    for u in tree.root.leaves:
        remaining = colors - set([u.color])             # colors to which no best match has yet been found
        parent = u.parent                               # start with direct parent of each node
        while remaining and parent:
            colors_here = set()
            for v in parent.leaves:
                if v.color in remaining:                # best match found
                    colors_here.add(v.color)
                    bmg.add_edge(u.ID, v.ID)            # insert edge (u,v)
            remaining -= colors_here                    # update remaining colors
            parent = parent.parent
    
    if not supply_rbmg:
        return bmg
    else:
        return bmg, bmg.to_undirected(reciprocal=True)
    
    
def bmg_from_tree_quadratic(tree, supply_rbmg=False):
    """Construct a BMG (and optionally RBMG) from a given tree in O(|L|^2).
    
    Implementation of the quadratic algorithm in Geiss et al. 2020. Proven to
    run in O(|L|^2).
    """
    
    tree.supply_leaves()
    bmg = nx.DiGraph()
    colors = set()
    
    # maps (v, color) --> 0 / 1
    l = {}
    
    for v in tree.root.leaves:
        colors.add(v.color)
        bmg.add_node(v.ID, label=v.label, color=v.color)
        l[v, v.color] = 1
        
    for v in tree.postorder():
        
        if not v.children:
            continue
        
        for u1, u2 in itertools.combinations(v.children, 2):
            for x, y in itertools.product(u1.leaves, u2.leaves):
                
                if not l.get((u1, y.color)):
                    bmg.add_edge(x.ID, y.ID)
                if not l.get((u2, x.color)):
                    bmg.add_edge(y.ID, x.ID)
                    
        for u, r in itertools.product(v.children, colors):
            if l.get((u, r)):
                l[v, r] = 1
    
    if not supply_rbmg:
        return bmg
    else:
        return bmg, bmg.to_undirected(reciprocal=True)
    

# --------------------------------------------------------------------------
#                 Analysis of the structure of (R)BMGs
# --------------------------------------------------------------------------


def classify_good_ugly(bmg, rbmg, fp):
    
    for x, y in fp.edges():
        fp.edges[x, y]['middle_in_good'] = False
        fp.edges[x, y]['first_in_ugly'] = False
        fp.edges[x, y]['first_in_bad'] = False
        
    color_dict = sort_by_colors(bmg)
    
    for x, y in fp.edges():
        
        # middle edges of good quartets
        for c in color_dict:
            
            if c == bmg.nodes[x]['color'] or c == bmg.nodes[y]['color']:
                continue
            
            for z1, z2 in itertools.permutations(color_dict[c], 2):
                
                if (rbmg.has_edge(z1, x) and rbmg.has_edge(z2, y) and
                    bmg.has_edge(z1, y) and not bmg.has_edge(y, z1) and
                    bmg.has_edge(z2, x) and not bmg.has_edge(x, z2)):
                    
                    fp.edges[x, y]['middle_in_good'] = True
        
        # first edges of ugly/bad quartets          
        for _ in range(2):
            
            for x2 in color_dict[bmg.nodes[x]['color']]:
                if x == x2:
                    continue
                
                for z in rbmg.neighbors(x2):
                    
                    if (bmg.nodes[z]['color'] == bmg.nodes[x]['color'] or
                        bmg.nodes[z]['color'] == bmg.nodes[y]['color']):
                        continue
                    # first edges of ugly quartets
                    if (rbmg.has_edge(y, x2) and not rbmg.has_edge(z, y) and
                        not rbmg.has_edge(z, x)):
                        fp.edges[x, y]['first_in_ugly'] = True
                    
                    # first edges of bad quartets
                    elif (not rbmg.has_edge(y, x2) and rbmg.has_edge(z, y) and
                          not bmg.has_edge(x, z)):
                        fp.edges[x, y]['first_in_bad'] = True
            
            # swap for second iteration
            x, y = y, x
    
    return fp   


def count_good_ugly(bmg, rbmg, fp):
    
    classify_good_ugly(bmg, rbmg, fp)
    
    good, ugly, both = 0, 0, 0
    
    for x, y in fp.edges():
        
        if fp.edges[x, y]['middle_in_good'] and fp.edges[x, y]['first_in_ugly']:
            good += 1
            ugly += 1
            both += 1
        elif fp.edges[x, y]['middle_in_good']:
            good += 1
        elif fp.edges[x, y]['first_in_ugly']:
            ugly += 1
    
    return good, ugly, both


def count_good_ugly_bad(bmg, rbmg, fp):
    
    classify_good_ugly(bmg, rbmg, fp)
    
    gub, gu, gb, ub, g, u, b = 0, 0, 0, 0, 0, 0, 0
    
    for x, y in fp.edges():
        
        if (fp.edges[x, y]['middle_in_good'] and
            fp.edges[x, y]['first_in_ugly'] and
            fp.edges[x, y]['first_in_bad']):
            gub += 1
            gu += 1
            gb += 1
            ub += 1
            g += 1
            u += 1
            b += 1
            
        elif fp.edges[x, y]['middle_in_good'] and fp.edges[x, y]['first_in_ugly']:
            gu += 1
            g += 1
            u += 1
        elif fp.edges[x, y]['middle_in_good'] and fp.edges[x, y]['first_in_bad']:
            gb += 1
            g += 1
            b += 1
        elif fp.edges[x, y]['first_in_ugly'] and fp.edges[x, y]['first_in_bad']:
            ub += 1
            u += 1
            b += 1
        elif fp.edges[x, y]['middle_in_good']:
            g += 1
        elif fp.edges[x, y]['first_in_ugly']:
            u += 1
        elif fp.edges[x, y]['first_in_bad']:
            b += 1
    
    return gub, gu, gb, ub, g, u, b
    

# --------------------------------------------------------------------------
#                             P4 EDITING
# --------------------------------------------------------------------------
    
def list_path(t, k, path, G, P4_list):
    for u in G[path[-1]]:
        G.nodes[u]['vis'] += 1
    for u in G[path[-1]]:
        if G.nodes[u]['vis'] == 1:
            if t < k:
                list_path(t+1, k, path + [u], G, P4_list)
            elif path[0] < u:
                P4_list.append(tuple(path) + (u,))
    for u in G[path[-1]]:
        G.nodes[u]['vis'] -= 1


def find_all_P4(G):
    P4_list = []
    for v in G:
            G.nodes[v]['vis'] = 0
    for v in G:
        G.nodes[v]['vis'] = 1
        list_path(2, 4, [v], G, P4_list)
        G.nodes[v]['vis'] = 0
    return P4_list


def is_good_quartet(path, bmg):
    """Determine whether an induced P4 is a good quartet in the bmg."""
    if (bmg.has_edge(path[0], path[2]) and 
        bmg.has_edge(path[3], path[1]) and
        bmg.nodes[path[0]]['color'] == bmg.nodes[path[3]]['color']):
        return True
    else:
        return False
    

def remove_all_P4(rbmg, bmg, P4_list=None):
    """Find all good quartets and removes the middle edges."""
    GP4 = rbmg.copy()
    
    if P4_list is None:
        P4_list = find_all_P4(GP4)

    for path in P4_list:
        if is_good_quartet(path, bmg) and GP4.has_edge(path[1], path[2]):
            GP4.remove_edge(path[1], path[2])
            print(path[1], "-|-", path[2])
    
    return GP4

# --------------------------------------------------------------------------
#                            C6 DETECTION
# --------------------------------------------------------------------------

def validate_C6(G, C6_path):
    
    C6 = G.subgraph(C6_path)
    if C6.size() != 6:
        print(C6_path)
        print('Nodes', [n for n in C6.nodes()], C6.order())
        print('Edges', [edge for edge in C6.edges()], C6.size())
        print('Size not equal to 6!')
        return False
    for n in C6.nodes():
        if not C6.degree(n) == 2:
            print("Degree not equal to 2 for", n)
            return False
        for neighbor in G.neighbors(n):
            if G.nodes[n]['color'] == G.nodes[neighbor]['color']:
                print('Same color', n, neighbor)
                return False
    return True


def find_all_C6(G, P4_list=None):
    """Find induces C6 of the form <x1 y1 z1 x2 y2 z2>."""
    
    if P4_list is None:
        P4_list = find_all_P4(G)
        
    endpoints = {}
    for path in P4_list:
        a, b, c, d = path if path[0] < path[3] else path[::-1]
        a_col, b_col, c_col, d_col = (G.nodes[a]['color'], G.nodes[b]['color'],
                                      G.nodes[c]['color'], G.nodes[d]['color'])
        # avoid multiple recording for the 3 different colors:
        if a_col == d_col and a_col < b_col and a_col < c_col:
            if (a, d) not in endpoints:
                endpoints[(a, d)] = {}
            if (b_col, c_col) not in endpoints[(a, d)]:
                endpoints[(a, d)][(b_col, c_col)] = []
            endpoints[(a, d)][(b_col, c_col)].append(path)
    
    C6_list = []
    for endpoint, colors in endpoints.items():
        x1, x2 = endpoint
        for color_pair, paths in colors.items():
            y_color, z_color = color_pair
            # avoid multiple recording for the 2 directions:
            if y_color < z_color and (z_color, y_color) in colors:
                for path1 in paths:
                    for path2 in colors[(z_color, y_color)]:
                        if ((not G.has_edge(path1[1], path2[1])) and
                            (not G.has_edge(path1[2], path2[2]))):
                            C6_list.append(path1 + (path2[2], path2[1]))
                            if not validate_C6(G, C6_list[-1]):
                                raise RuntimeError('problem in C6 detection')
                            print(C6_list[-1])
    return C6_list


def graph_type(G):
    
    P4_list = find_all_P4(G)
    C6_list = find_all_C6(G, P4_list=P4_list)
    
    if C6_list:
        return 'C', P4_list, C6_list
    elif P4_list:
        return 'B', P4_list, C6_list
    else:
        return 'A', P4_list, C6_list