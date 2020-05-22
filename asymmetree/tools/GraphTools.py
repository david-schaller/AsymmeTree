# -*- coding: utf-8 -*-

"""
Graph Tools.
"""

import itertools
import numpy as np
import networkx as nx


__author__ = "David Schaller"


# --------------------------------------------------------------------------
#                            GENERAL GRAPH TOOLS
# --------------------------------------------------------------------------

def graphs_equal(graph1, graph2):
    """Returns whether two NetworkX graphs (directed or undirected) are equal."""
    
    if (not graph1.order() == graph2.order()) or (not graph1.size() == graph2.size()):
        return False
    
    for u, v in graph1.edges():
        if not graph2.has_edge(u,v):
            return False
    
    return True


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
    
    index = {node: i for i, node in enumerate(G.nodes())}   # maps node --> row/column index
    M = np.zeros( (len(index), len(index)), dtype=np.int8)
    
    for x, neighbors in G.adjacency():
        for y in neighbors:
            M[index[x], index[y]] = 1                       # add out-neighbor key_j to key_i
            
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
            

def classify_good_ugly(BMG, RBMG, fp):
    
    for x, y in fp.edges():
        fp.edges[x, y]['middle_in_good'] = False
        fp.edges[x, y]['first_in_ugly'] = False
        
    color_dict = sort_by_colors(BMG)
    
    
    for x, y in fp.edges():
        
        # middle edges of good quartets
        for c in color_dict:
            
            if c == BMG.nodes[x]['color'] or c == BMG.nodes[y]['color']:
                continue
            
            for z1, z2 in itertools.permutations(color_dict[c], 2):
                
                if (RBMG.has_edge(z1, x) and RBMG.has_edge(z2, y) and
                    BMG.has_edge(z1, y) and not BMG.has_edge(y, z1) and
                    BMG.has_edge(z2, x) and not BMG.has_edge(x, z2)):
                    
                    fp.edges[x, y]['middle_in_good'] = True
                    
        # first edges of ugly quartets
        for _ in range(2):
            
            for x2 in color_dict[BMG.nodes[x]['color']]:
                if x == x2 or not RBMG.has_edge(y, x2):
                    continue
                
                for z in RBMG.neighbors(x2):
                    
                    if (BMG.nodes[z]['color'] == BMG.nodes[x]['color'] or
                        BMG.nodes[z]['color'] == BMG.nodes[y]['color']):
                        continue
                    
                    if not RBMG.has_edge(z, y) and not RBMG.has_edge(z, x):
                        fp.edges[x, y]['first_in_ugly'] = True
            
            # swap for second iteration
            x, y = y, x
    
    return fp   


def count_good_ugly(BMG, RBMG, fp):
    
    classify_good_ugly(BMG, RBMG, fp)
    
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


def is_good_quartet_OLD(path, BMG):
    """Determine whether an induced P4 is a good quartet in the BMG."""
    if (BMG.has_edge(path[0], path[2]) and 
        BMG.has_edge(path[3], path[1]) and
        BMG.nodes[path[0]]['color'] == BMG.nodes[path[3]]['color']):
        return True
    else:
        return False


def is_good_quartet(path, BMG):
    """Determine whether an induced P4 is a good quartet in the BMG."""
    induced = BMG.subgraph(path)
    degrees = [(induced.out_degree(v), induced.in_degree(v)) for v in induced]
    if (degrees[0] == (2,1) and 
        degrees[1] == (2,3) and
        degrees[2] == (2,3) and 
        degrees[3] == (2,1) and
        BMG.nodes[path[0]]['color'] == BMG.nodes[path[3]]['color']):
        return True
    else:
        return False
    

def remove_all_P4(RBMG, BMG, P4_list=None):
    """Find all good quartets and removes the middle edges."""
    GP4 = RBMG.copy()
    
    if P4_list is None:
        P4_list = find_all_P4(GP4)

    for path in P4_list:
        if is_good_quartet(path, BMG) and GP4.has_edge(path[1], path[2]):
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
                                raise Exception("Problem in C6 detection")
                            print(C6_list[-1])
    return C6_list


def graph_type(G):
    
    P4_list = find_all_P4(G)
    C6_list = find_all_C6(G, P4_list=P4_list)
    
    if C6_list:
        return "C", P4_list, C6_list
    elif P4_list:
        return "B", P4_list, C6_list
    else:
        return "A", P4_list, C6_list
