# -*- coding: utf-8 -*-

"""
Construction of Least Resolved Trees (LRT).

Computes an LRT from a (not necessarily valid) cBMG or
from an observable gene tree.
"""

import itertools

import networkx as nx

from asymmetree.best_matches.TrueBMG import bmg_from_tree
from asymmetree.datastructures.PhyloTree import PhyloTree, PhyloTreeNode
from asymmetree.tools.GraphTools import (sort_by_colors,
                                         is_properly_colored,
                                         graphs_equal)
from asymmetree.tools.Build import Build


__author__ = 'David Schaller'


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
    
    if isinstance(tree, PhyloTreeNode):
        tree = PhyloTree(tree)
    
    # assign IDs to inner nodes
    tree.reconstruct_IDs()
    # assign label and colors to the leaves
    tree.reconstruct_info_from_graph(G)
    
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
    lrt.remove_planted_root()
    
    # assign list of leaves to each node
    lrt.supply_leaves()
    
    subtree_colors = {}
    for v in lrt.preorder():
        subtree_colors[v] = {leaf.color for leaf in v.leaves}
        
    arc_colors = _arc_colors(lrt, subtree_colors)
    red_edges = redundant_edges(lrt, subtree_colors, arc_colors)
    lrt.contract(red_edges)
    lrt = lrt.topology_only()
    
    return lrt


def _arc_colors(T, subtree_colors):
    """Color sets relevant for redundant edge computation.
    
    Computes for all inner vertices v the color set of y such that y with (x,y)
    is an arc in the BMG and lca(x,y) = v.
    """
    
    all_colors = subtree_colors[T.root]
    
    # color sets for all v
    arc_colors = {v: set() for v in T.preorder()}
    
    for u in T.root.leaves:
        
        # colors to which no best match has yet been found
        remaining = all_colors - {u.color}
        
        # start with direct parent of each node
        current = u.parent
        
        while remaining and current:
            colors_here = set()
            for v in current.leaves:
                
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
            subtree.reconstruct_info_from_graph(sg)
            if not graphs_equal(sg, bmg_from_tree(subtree)):
                return False
            subtrees.append(subtree)
    
    if len(subtrees) == 1:
        root = subtrees[0].root
    else:
        root = PhyloTreeNode(-1)
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
        tree = PhyloTree(PhyloTreeNode(-1))
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
            root = PhyloTreeNode(-1)
            for v in self.digraph.nodes():
                root.add_child(PhyloTreeNode(v))
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
                root = PhyloTreeNode(-1)
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
            
        node = PhyloTreeNode(-1)
        for v in S_2:
            node.add_child(PhyloTreeNode(v))
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

        
if __name__ == '__main__':
    
    import time
    
    N = 100
    
    T = PhyloTree.random_colored_tree(N, 2)
    print('--- T ---\n', T.to_newick())
    
    bmg = bmg_from_tree(T)
    
    start_time1 = time.time()
    lrt1 = lrt_from_colored_graph(bmg, mincut=False)
    # lrt1 = is_bmg(bmg)
    end_time1 = time.time()
    
    lrt2 = lrt_from_observable_tree(T)
    
    print('--- LRT1 ---\n', lrt1.to_newick())
    print('--- LRT2 ---\n', lrt2.to_newick())
    print('LRTs equal: {}'.format( lrt1.compare_topology(lrt2) ))
    
    bmg = bmg_from_tree(T)
    
    start_time2 = time.time()
    tc = TwoColoredLRT(bmg)
    lrt3 = tc.build_tree()
    end_time2 = time.time()
    
    lrt4 = lrt_from_observable_tree(T)
    print('--- LRT3 ---\n', lrt3.to_newick())
    print('--- LRT4 ---\n', lrt4.to_newick())
    print('LRTs equal: {}'.format( lrt3.compare_topology(lrt4) ))
    print(end_time1 - start_time1, end_time2 - start_time2)