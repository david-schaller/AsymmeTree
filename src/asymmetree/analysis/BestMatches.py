# -*- coding: utf-8 -*-

"""
Orthology graph, BMG, RBMG.

This module provides classes concerning colored graphs and colored
(phylogenetic) trees, including a BMG(Tree)-function, and computation of an
LRT from a (not necessarily valid) BMG or from an observable gene tree.
"""

import itertools

import networkx as nx

from tralda.datastructures.Tree import Tree, TreeNode

from tralda.tools.GraphTools import (sort_by_colors,
                                     is_properly_colored,
                                     graphs_equal)
from tralda.supertree import Build

from asymmetree.tools.PhyloTreeTools import (topology_only,
                                             reconstruct_color_from_graph)


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                 True Orthology Graph (TOG), BMG & RBMG
#
#                         (given a known tree)
# --------------------------------------------------------------------------

def orthology_from_tree(tree):
    """Constructs the true orthology graph from an event-labeled tree.
    
    Parameters
    ----------
    tree: tree
        A tree with leaves that have the `label`, and inner vertices that have
        the `event` attribute set.
    
    Returns
    -------
    networx.Graph
        The ortholy graph explained by the tree, i.e., a graph on the set of
        leaf labels and edges between two of them if and only if their last
        common ancestor was a speciation event.
    """
    
    if not isinstance(tree, Tree):
        raise TypeError("not of type 'Tree'")
    
    leaves = tree.leaf_dict()
    TOG = nx.Graph()
    
    for v in leaves[tree.root]:
        TOG.add_node(v.label, color=v.color)
    
    for node in tree.preorder():
        if node.event == 'S':
            for c1, c2 in itertools.combinations(node.children, 2):
                for u in leaves[c1]:
                    for v in leaves[c2]:
                        TOG.add_edge(u.label, v.label)
                        
    return TOG


def bmg_from_tree(tree, supply_rbmg=False):
    """Construct a BMG (and optionally RBMG) from a given tree.
    
    Parameters
    ----------
    tree: tree
        A tree with leaves that have the `label` and `color` attribute set.
    supply_rbmg : bool, optional
        If True, also return the symmetric part of the constructed BMG.
        The default is False.
    
    Returns
    -------
    networx.DiGraph or pair of networx.DiGraph and networx.Graph
        The constructed BMG and optionally its symmetric part.
    """
    
    if not isinstance(tree, Tree):
        raise TypeError("not of type 'Tree'")
    
    leaves = tree.leaf_dict()
    bmg = nx.DiGraph()
    colors = set()
    
    for v in leaves[tree.root]:
        colors.add(v.color)
        bmg.add_node(v.label, color=v.color)
    
    for u in leaves[tree.root]:
        remaining = colors - set([u.color])             # colors to which no best match has yet been found
        parent = u.parent                               # start with direct parent of each node
        while remaining and parent:
            colors_here = set()
            for v in leaves[parent]:
                if v.color in remaining:                # best match found
                    colors_here.add(v.color)
                    bmg.add_edge(u.label, v.label)      # insert edge (u,v)
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
    
    Parameters
    ----------
    tree: tree
        A tree with leaves that have the `label` and `color` attribute set.
    supply_rbmg : bool, optional
        If True, also return the symmetric part of the constructed BMG.
        The default is False.
    
    Returns
    -------
    networx.DiGraph or pair of networx.DiGraph and networx.Graph
        The constructed BMG and optionally its symmetric part.
    """
    
    if not isinstance(tree, Tree):
        raise TypeError("not of type 'Tree'")
    
    leaves = tree.leaf_dict()
    bmg = nx.DiGraph()
    colors = set()
    
    # maps (v, color) --> 0 / 1
    l = {}
    
    for v in leaves[tree.root]:
        colors.add(v.color)
        bmg.add_node(v.label, color=v.color)
        l[v, v.color] = 1
        
    for v in tree.postorder():
        
        if not v.children:
            continue
        
        for u1, u2 in itertools.combinations(v.children, 2):
            for x, y in itertools.product(leaves[u1], leaves[u2]):
                
                if not l.get((u1, y.color)):
                    bmg.add_edge(x.label, y.label)
                if not l.get((u2, x.color)):
                    bmg.add_edge(y.label, x.label)
                    
        for u, r in itertools.product(v.children, colors):
            if l.get((u, r)):
                l[v, r] = 1
    
    if not supply_rbmg:
        return bmg
    else:
        return bmg, bmg.to_undirected(reciprocal=True)
    
def extended_best_hits(leaves, D, epsilon=1e-8, supply_rbmg=False):
    """Compute BMG and RBMG from a distances matrix D.
    
    Parameters
    ----------
    leaves : list of TreeNode instances with `label` attribute
    D : 2-dimensional numpy array
        The distances of the leaves where the order of the rows/columns
        corresponds to that in the list.
    epsilon : float, optional
        relative threshold for the inclusion of suboptimal hits, i.e.,
        (x,y) in BMG if D(x,y) <= (1+epsilon) * min d(x,y'). The default
        is 1e-8 (for limited float precision).
    supply_rbmg : bool, optional
        If True, also return the symmetric part of the constructed BMG.
        The default is False.
    
    Returns
    -------
    networx.DiGraph or pair of networx.DiGraph and networx.Graph
        The constructed BMG and optionally its symmetric part.
    """
    
    bmg = nx.DiGraph()
    colors = set()
    relative_threshold = 1 + epsilon
    
    for v in leaves:
        bmg.add_node(v.label, color=v.color)
        colors.add(v.color)
    
    # ---- build BMG ----
    for u in range(len(leaves)):
        minima = {color: float('inf') for color in colors}
        for v in range(len(leaves)):
            if D[u,v] < minima[leaves[v].color]:
                minima[leaves[v].color] = D[u,v]
        for v in range(len(leaves)):
            if (leaves[u].color != leaves[v].color and
                D[u,v] <= relative_threshold * minima[leaves[v].color]):
                bmg.add_edge(leaves[u].label, leaves[v].label,
                             distance = D[u,v])
    
    if not supply_rbmg:
        return bmg
    else:
        return bmg, bmg.to_undirected(reciprocal=True)


# --------------------------------------------------------------------------
#                         Least resolved trees
# --------------------------------------------------------------------------

def informative_triples(graph, color_dict=None):
    """Compute the informative triples of a colored digraph."""
    
    if not isinstance(graph, nx.DiGraph):
        raise TypeError("must be a NetworkX 'Digraph'")
    
    if not color_dict:
        color_dict = sort_by_colors(graph)
    
    R = []
    
    for c1, c2 in itertools.permutations(color_dict.keys(), 2):
        for a in color_dict[c1]:
            for b1, b2 in itertools.permutations(color_dict[c2], 2):
                if graph.has_edge(a, b1) and (not graph.has_edge(a, b2)):
                    R.append( (a, b1, b2) )
    
    return R


def forbidden_triples(graph, color_dict=None):
    """Compute the forbidden triples of a colored digraph."""
    
    if not isinstance(graph, nx.DiGraph):
        raise TypeError("must be a NetworkX 'Digraph'")
    
    if not color_dict:
        color_dict = sort_by_colors(graph)
    
    F = []
    
    for c1, c2 in itertools.permutations(color_dict.keys(), 2):
        for a in color_dict[c1]:
            for b1, b2 in itertools.combinations(color_dict[c2], 2):
                if graph.has_edge(a, b1) and graph.has_edge(a, b2):
                    F.append( (a, b1, b2) )
                    F.append( (a, b2, b1) )
    
    return F


def informative_forbidden_triples(graph, color_dict=None):
    """Compute the informative and forbidden triples of a colored digraph."""
    
    if not isinstance(graph, nx.DiGraph):
        raise TypeError("must be a NetworkX 'Digraph'")
    
    if not color_dict:
        color_dict = sort_by_colors(graph)
    
    R, F = [], []
    
    for c1, c2 in itertools.permutations(color_dict.keys(), 2):
        for a in color_dict[c1]:
            for b1, b2 in itertools.combinations(color_dict[c2], 2):
                
                if graph.has_edge(a, b1) and graph.has_edge(a, b2):
                    F.append( (a, b1, b2) )
                    F.append( (a, b2, b1) )
                elif graph.has_edge(a, b1):
                    R.append( (a, b1, b2) )
                elif graph.has_edge(a, b2):
                    R.append( (a, b2, b1) )
    
    return R, F


def binary_explainable_triples(graph, color_dict=None):
    """Extended informative triple set for binary-explainable graphs."""
    
    if not isinstance(graph, nx.DiGraph):
        raise TypeError("must be a NetworkX 'Digraph'")
    
    if not color_dict:
        color_dict = sort_by_colors(graph)
    
    R_binary = []
    
    for c1, c2 in itertools.permutations(color_dict.keys(), 2):
        for a in color_dict[c1]:
            for b1, b2 in itertools.combinations(color_dict[c2], 2):
                
                if graph.has_edge(a, b1) and graph.has_edge(a, b2):
                    try:
                        if b2 < b1:
                            b1, b2 = b2, b1
                    except NotImplemented: pass
                    R_binary.append( (b1, b2, a) )
                elif graph.has_edge(a, b1):
                    R_binary.append( (a, b1, b2) )
                elif graph.has_edge(a, b2):
                    R_binary.append( (a, b2, b1) )
    
    return R_binary


def _finalize(tree, G):
    
    if not tree:
        return None
    
    if isinstance(tree, TreeNode):
        tree = Tree(tree)
    
    # assign colors to the leaves
    reconstruct_color_from_graph(tree, G)
    
    return tree


# --------------------------------------------------------------------------
#                      LRT FROM OBSERVABLE GENE TREE
# --------------------------------------------------------------------------

def lrt_from_observable_tree(T):
    """Computes the Least Resolved Tree from a tree.
    
    The unique Least Resolved Tree from a leaf-colored (observable)
    gene tree is computed by contraction of all redundant edges.
    """
    
    lrt = T.copy()
    if not lrt.root:
        return lrt
    
    # remove planted root if existent
    if len(lrt.root.children) == 1:
        new_root = lrt.root.children[0]
        new_root.detach()
        lrt.root = new_root
        new_root.dist = 0.0
    
    # assign list of leaves to each node
    leaves = lrt.leaf_dict()
    
    subtree_colors = {}
    for v in lrt.preorder():
        subtree_colors[v] = {leaf.color for leaf in leaves[v]}
        
    arc_colors = _arc_colors(lrt, leaves, subtree_colors)
    red_edges = redundant_edges(lrt, subtree_colors, arc_colors)
    lrt.contract(red_edges)
    lrt = topology_only(lrt)
    
    return lrt


def _arc_colors(T, leaves, subtree_colors):
    """Color sets relevant for redundant edge computation.
    
    Computes for all inner vertices v the color set of y such that y with (x,y)
    is an arc in the BMG and lca(x,y) = v.
    """
    
    all_colors = subtree_colors[T.root]
    
    # color sets for all v
    arc_colors = {v: set() for v in T.preorder()}
    
    for u in leaves[T.root]:
        
        # colors to which no best match has yet been found
        remaining = all_colors - {u.color}
        
        # start with direct parent of each node
        current = u.parent
        
        while remaining and current:
            colors_here = set()
            for v in leaves[current]:
                
                # best match found
                if v.color in remaining:
                    colors_here.add(v.color)
            
            arc_colors[current].update(colors_here)
            remaining -= colors_here
            current = current.parent
    
    return arc_colors


def redundant_edges(T, subtree_colors, arc_colors):
    
    red_edges = []
    
    for u, v in T.inner_edges():
        
        # colors s in sigma( L(T(u) \ T(v)) )
        aux_set = set()
        
        for v2 in u.children:
            if v2 is not v:
                aux_set.update(subtree_colors[v2])
        
        if not arc_colors[v].intersection(aux_set):
            red_edges.append((u, v))
    
    return red_edges


# --------------------------------------------------------------------------
#                               LRT FROM BMG
# --------------------------------------------------------------------------

def is_bmg(G):
    """Determine whether a colored digraph is a BMG.
    
    If the graph is a BMG, then its LRT is returned.
    """
    
    if not isinstance(G, nx.DiGraph):
            raise TypeError('not a digraph')
    if not is_properly_colored(G):
        raise RuntimeError('not a properly colored digraph')
    
    subtrees = []
    colors = None
    
    for wcc in nx.weakly_connected_components(G):
        
        sg = G.subgraph(wcc).copy()
        color_dict = sort_by_colors(sg)
        
        # the color sets of all components must be equal
        if colors is None:
            colors = set(color_dict)
        elif colors != set(color_dict):
            return False
        
        L = {v for v in sg.nodes()}
        R = informative_triples(sg, color_dict=color_dict)
        build = Build(R, L, mincut=False)
        subtree = build.build_tree()
        
        if not subtree: 
            return False
        else:
            # a digraph is a BMG iff its equal to the BMG of BUILD(R)
            reconstruct_color_from_graph(subtree, sg)
            if not graphs_equal(sg, bmg_from_tree(subtree)):
                return False
            subtrees.append(subtree)
    
    if len(subtrees) == 1:
        root = subtrees[0].root
    else:
        root = TreeNode()
        for subtree in subtrees:
            root.add_child(subtree.root)
    
    return _finalize(root, G)
    

def lrt_from_colored_graph(G, mincut=False, weighted_mincut=False):
    
    L = {v for v in G.nodes()}
    R = informative_triples(G)
    
    build = Build(R, L, mincut=mincut, weighted_mincut=weighted_mincut)
    tree = build.build_tree()
    
    return _finalize(tree, G)


def correct_bmg(bmg_original):
    """Build the LRT (using min cut in BUILD algorithm) and return its BMG."""
    
    subtrees = []
    for sg in (bmg_original.subgraph(c)
               for c in nx.weakly_connected_components(bmg_original)):
        tree = lrt_from_colored_graph(sg, mincut=True)
        if tree:
            subtrees.append(tree)
            
    if len(subtrees) == 0:
        return None
    elif len(subtrees) == 1:
        tree = subtrees[0]
    else:
        tree = Tree(TreeNode())
        for subtree in subtrees:
            tree.root.add_child(subtree.root)
    
    return bmg_from_tree(tree)


# --------------------------------------------------------------------------
#                          LRT FROM 2-col. BMG
#                         (new characterization)
# --------------------------------------------------------------------------
            
class TwoColoredLRT:
    
    def __init__(self, digraph):
        
        if not isinstance(digraph, nx.DiGraph):
            raise TypeError('not a digraph')
            
        self.digraph = digraph
        self.color_dict = sort_by_colors(digraph)
                
    
    def build_tree(self):
        
        if not is_properly_colored(self.digraph):
            raise RuntimeError('not a properly colored digraph')
        if len(self.color_dict) > 2:
            raise RuntimeError('more than 2 colors')
        
        # star tree if there is only one color
        if len(self.color_dict) == 1:
            root = TreeNode()
            for v in self.digraph.nodes():
                root.add_child(TreeNode(label=v))
        # 2 colors
        else:
            subtrees = []
            for wcc in nx.weakly_connected_components(self.digraph):
                if len(wcc) == 1:
                    return False
                
                subroot = self._build_tree(self.digraph.subgraph(wcc).copy())
                
                if not subroot:
                    return False
                else:
                    subtrees.append(subroot)
            
            if len(subtrees) == 1:
                root = subtrees[0]
            else:
                root = TreeNode()
                for subroot in subtrees:
                    root.add_child(subroot)
        
        return _finalize(root, self.digraph)
    
        
    def _build_tree(self, G):
            
        color_count = {color: 0 for color in self.color_dict.keys()}
        for v in G.nodes():
            color_count[G.nodes[v]['color']] += 1
        other_color = {color: G.order() - count
                       for color, count in color_count.items()}
        
        umbrella = {v for v in G.nodes() if 
                    other_color[G.nodes[v]['color']] == G.out_degree(v)}
        S_1 = {v for v in umbrella if umbrella.issuperset(G.predecessors(v))}
        S_2 = {v for v in S_1 if S_1.issuperset(G.predecessors(v))}
        
        if not S_2 or len(S_1) != len(S_2):
            return False
            
        node = TreeNode()
        for v in S_2:
            node.add_child(TreeNode(label=v))
            G.remove_node(v)
        
        for wcc in nx.weakly_connected_components(G):
            
            if len(wcc) == 1:
                return False
            
            child = self._build_tree(G.subgraph(wcc).copy())
            
            if not child:
                return False
            else:
                node.add_child(child)
                
        return node
    

def lrt_from_2bmg(G):
    """Effieciently constructs the LRT for a 2-BMG via the support vertices.
    
    Returns false if the graph is not a 2-colored BMG.
    """
    
    builder = TwoColoredLRT(G)
    return builder.build_tree()


# --------------------------------------------------------------------------
#                      Binary-refinable tree (BRT)
# --------------------------------------------------------------------------

def binary_refinable_tree(G, mincut=False, weighted_mincut=False):
    """Compute the binary-refinable tree (BRT) for a graph.
    
    The BRT of a BMG can, if it exists, be refined arbitrarily and such that
    it still explains the BMG. In particular, all binary explainations are 
    refinements of the BRT."""
    
    L = {v for v in G.nodes()}
    R = binary_explainable_triples(G)
    
    build = Build(R, L, mincut=mincut, weighted_mincut=weighted_mincut)
    brt = build.build_tree()
    
    return _finalize(brt, G)


# --------------------------------------------------------------------------
#                           Augmented tree
# --------------------------------------------------------------------------

def augment_and_label(tree, inplace=False):
    """Augment tree and add event labeling based on color intersections."""
    
    if not inplace:
        tree = tree.copy()
    
    # precompute since tree will be modified
    inner_nodes = [u for u in tree.inner_nodes()]
    leaves = tree.leaf_dict()
    
    for u in inner_nodes:
        
        # should not appear except if tree is planted, otherwise ignore
        if len(u.children) == 1:
            continue
        
        C = _color_intersection_components(u, leaves)
        
        if len(C) == 1:
            u.event = 'D'
            
        else:
            u.event = 'S'
            
            for cc in C:
                if len(cc) > 1:
                    w = TreeNode(event='D')
                    u.add_child(w)
                    for v in cc:
                        u.remove_child(v)
                        w.add_child(v)
        
    return tree


def _color_intersection_components(u, leaves):
    
    aux_graph = nx.Graph()
    
    for v in u.children:
        aux_graph.add_node(v)
        
        for x in leaves[v]:
            aux_graph.add_edge(v, x.color)
            
    result = []
    
    for cc in nx.connected_components(aux_graph):
        
        # discard the colors in the components
        result.append( [v for v in cc if isinstance(v, TreeNode)] )
        
    return result


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