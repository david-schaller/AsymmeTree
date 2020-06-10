# -*- coding: utf-8 -*-

"""
Graph Tools.
"""

import itertools, random
import numpy as np
import networkx as nx


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                            GENERAL GRAPH TOOLS
# --------------------------------------------------------------------------

def graphs_equal(G1, G2):
    """Returns whether two NetworkX graphs (directed or undirected) are equal."""
    
    if (not G1.order() == G2.order()) or (not G1.size() == G2.size()):
        return False
    
    if set(G1.nodes()) != set(G2.nodes()):
        return False
    
    for u, v in G1.edges():
        if not G2.has_edge(u,v):
            return False
    
    return True


def symmetric_part(G):
    """Returns the (undirected) symmetric part of a directed graph."""
    
    if not isinstance(G, nx.DiGraph):
        raise TypeError('directed graph required')
    
    G_sym = nx.Graph()
    G_sym.add_nodes_from(G.nodes(data=True))
    
    for x, neighbors in G.adjacency():
        for y in neighbors:
            if G.has_edge(y,x):
                G_sym.add_edge(x,y)
                
    return G_sym


def symmetric_diff(G1, G2):
    """Returns the number of edges in the symmetric difference."""
    
    set1 = {x for x in G1.nodes()}
    set2 = {x for x in G2.nodes()}
    
    if set1 != set2:
        raise RuntimeError('graphs do not have the same vertex set')
        return
    
    sym_diff_number = 0
    
    for x, y in itertools.combinations(set1, 2):
        if G1.has_edge(x, y) and not G2.has_edge(x, y):
            sym_diff_number += 1
        elif not G1.has_edge(x, y) and G2.has_edge(x, y):
            sym_diff_number += 1
    
    return sym_diff_number


def random_graph(N, p=0.5):
    """Construct a random graph on N nodes.
    
    Keyword arguments:
        p - probability that an edge xy is inserted
    """
        
    G = nx.Graph()
    G.add_nodes_from(range(1, N+1))
    
    for x, y in itertools.combinations(range(1, N+1), 2):
        if random.random() < p:
            G.add_edge(x, y)
    
    return G


def performance(true_graph, graph):
    """Returns order, size, tp, tn, fp, fn, accuracy, precision and recall
    of a (directed or undirected) graph w.r.t. true graph."""
    
    order, size, tp, fp, fn = graph.order(), graph.size(), 0, 0, 0
    
    for u, v in graph.edges():
        if true_graph.has_edge(u,v):
            tp += 1
        else:
            fp += 1
    
    for u, v in true_graph.edges():
        if not graph.has_edge(u, v):
            fn += 1
    
    if isinstance(graph, nx.DiGraph):
        tn = (order * (order-1)) - (tp + fp + fn)
    else:
        tn = (order * (order-1) / 2) - (tp + fp + fn)
    
    accuracy = (tp + tn) / (tp + tn + fp + fn) if tp + tn + fp + fn > 0 else float('nan')
    precision = tp / (tp + fp) if tp + fp > 0 else float('nan')
    recall = tp / (tp + fn) if tp + fn > 0 else float('nan')
    
    return order, size, tp, tn, fp, fn, accuracy, precision, recall


def false_edges(true_graph, graph):
    """Returns a graph containing fn and a graph containg fp edges."""
    
    if isinstance(true_graph, nx.DiGraph):
        fn_graph = nx.DiGraph()
        fp_graph = nx.DiGraph()
    else:
        fn_graph = nx.Graph()
        fp_graph = nx.Graph()
        
    fn_graph.add_nodes_from(true_graph.nodes(data=True))
    fp_graph.add_nodes_from(true_graph.nodes(data=True))
        
    for u, v in graph.edges():
        if not true_graph.has_edge(u,v):
            fp_graph.add_edge(u, v)
    
    for u, v in true_graph.edges():
        if not graph.has_edge(u, v):
            fn_graph.add_edge(u, v)
        
    return fn_graph, fp_graph


def build_adjacency_matrix(G):
    """Return an adjacency matrix."""
    
    # maps node --> row/column index
    index = {node: i for i, node in enumerate(G.nodes())}
    M = np.zeros( (len(index), len(index)), dtype=np.int8)
    
    for x, neighbors in G.adjacency():
        for y in neighbors:
            M[index[x], index[y]] = 1
            
    return M


def sort_by_colors(graph):
    
    color_dict = {}
    
    for v in graph.nodes():
        color = graph.nodes[v]['color']
        
        if color not in color_dict:
            color_dict[color] = [v]
        else:
            color_dict[color].append(v)
    
    return color_dict
            

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
        print("Size not equal to 6!")
        return False
    for n in C6.nodes():
        if not C6.degree(n) == 2:
            print("Degree not equal to 2 for", n)
            return False
        for neighbor in G.neighbors(n):
            if G.nodes[n]['color'] == G.nodes[neighbor]['color']:
                print("Same color", n, neighbor)
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
