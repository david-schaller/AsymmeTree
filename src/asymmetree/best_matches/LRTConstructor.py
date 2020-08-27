# -*- coding: utf-8 -*-

"""
Construction of Least Resolved Trees (LRT).

Computes an LRT from a (not necessarily valid) cBMG or
from an observable gene tree.
"""

import itertools

import networkx as nx

from asymmetree.best_matches import TrueBMG
from asymmetree import PhyloTree, PhyloTreeNode
from asymmetree.tools.GraphTools import sort_by_colors


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
                elif graph.has_edge(a, b1):
                    R.append( (a, b1, b2) )
                elif graph.has_edge(a, b2):
                    R.append( (a, b2, b1) )
    
    return R, F


def aho_graph(R, L, weighted=False):
    """Construct the auxiliary graph (Aho graph) for BUILD.
        
    Edges (a,b) are optionally weighted by the number of occurrences in
    triples of the form ab|x or ba|x.
    """
    
    G = nx.Graph()
    G.add_nodes_from(L)

    for t1, t2, t3 in R:
        if not G.has_edge(t1, t2):
            G.add_edge(t1, t2, weight=1)
        elif weighted:
            G[t1][t2]['weight'] += 1
    
    return G
    


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
    redundant_edges = _redundant_edges(lrt, subtree_colors, arc_colors)
    lrt.contract(redundant_edges)
    lrt = lrt.topology_only()
    
    return lrt


def _arc_colors(T, subtree_colors):
    """Color sets relevant for redundant edge computation.
    
    Computes for all inner vertices v the color set of y such that y with (x,y)
    is an arc in the BMG and lca(x,y) = v.
    """
    
    all_colors = subtree_colors[T.root]
    
    arc_colors = {v: set() for v in T.preorder()}     # color sets for all v
    
    for u in T.root.leaves:
        remaining = all_colors - {u.color}            # colors to which no best match has yet been found
        current = u.parent                            # start with direct parent of each node
        while remaining and current:
            colors_here = set()
            for v in current.leaves:
                if v.color in remaining:              # best match found
                    colors_here.add(v.color)
            
            arc_colors[current].update(colors_here)
            remaining -= colors_here
            current = current.parent
    
    return arc_colors


def _redundant_edges(T, subtree_colors, arc_colors):
    
    redundant_edges = []
    
    for u, v in T.inner_edges():
        
        aux_set = set()                                 # colors s in sigma( L(T(u) \ T(v)) )
        for v2 in u.children:
            if v2 is not v:
                aux_set.update(subtree_colors[v2])
        
        if not arc_colors[v].intersection(aux_set):
            redundant_edges.append((u, v))
    
    return redundant_edges


# --------------------------------------------------------------------------
#                               LRT FROM BMG
# --------------------------------------------------------------------------

class LRTConstructor:
    
    def __init__(self, G, mincut=True, weighted_mincut=False):
        
        self.G = G
        self.mincut = mincut
        self.weighted_mincut = weighted_mincut
    
    
    def build_tree(self):
        """Build the least resolved tree from informative triples."""
        
        self.L = {v for v in self.G.nodes()}
        self.color_dict = sort_by_colors(self.G)
        self.R = informative_triples(self.G, color_dict=self.color_dict)
        
        self.cut_value = 0
        self.cut_list = []
        
        root = self.aho(self.L, self.R)
        if not root:
            print('no such tree exists')
            return None
        else:
            return PhyloTree(root)
    
    
    def aho(self, L, R):
        """Recursive Aho-algorithm."""
        
        # trivial case: only one leaf left in L
        if len(L) == 1:
            leaf = L.pop()
            return PhyloTreeNode(leaf, label=self.G.nodes[leaf]['label'],
                                 color=self.G.nodes[leaf]['color'])
            
        help_graph = aho_graph(R, L, weighted=self.weighted_mincut)
        conn_comps = self.connected_components(help_graph)
        
        # return False if less than 2 connected components
        if len(conn_comps) <= 1:
            print('Connected component:\n', conn_comps)
            return False
        
        child_nodes = []
        for cc in conn_comps:
            Li = set(cc)                            # list/dictionary --> set
            Ri = []
            for t in R:                             # determine which triples are in the subtree
                if Li.issuperset(t):
                    Ri.append(t)
            Ti = self.aho(Li, Ri)                   # recursive call
            if not Ti:
                return False                        # raise False to previous call of aho()
            else:
                child_nodes.append(Ti)
                
        subtree_root = PhyloTreeNode(0, label='')   # place new inner node
        for Ti in child_nodes:
            subtree_root.add_child(Ti)              # add roots of the subtrees to the new node
   
        return subtree_root                         # return the new node
    
    
    def connected_components(self, aho_graph):
        """Determines the connected components of the graph."""
        
        conn_comps = list(nx.connected_components(aho_graph))
        if (not self.mincut) or len(conn_comps) > 1:
            return conn_comps
        else:
            # Stoer–Wagner algorithm
            cut_value, partition = nx.stoer_wagner(aho_graph)
            self.cut_value += cut_value
            if len(partition[0]) < len(partition[1]):
                smaller_comp = partition[0]
            else:
                smaller_comp = partition[1]
            for edge in self.G.edges:
                if ((edge[0] in smaller_comp and edge[1] not in smaller_comp)
                    or
                    (edge[1] in smaller_comp and edge[0] not in smaller_comp)):
                    self.cut_list.append(edge)
            return partition
    
    
    @staticmethod
    def correct_bmg(bmg_original):
        """Build the LRT (using Mincut in Aho algorithm) and
        return its BMG and RBMG."""
        
        subtrees = []
        for sg in (bmg_original.subgraph(c)
                   for c in nx.weakly_connected_components(bmg_original)):
            lrt_constructor = LRTConstructor(sg, mincut=True)
            tree = lrt_constructor.build_tree()
            if tree:
                subtrees.append(tree)
                
        if len(subtrees) == 0:
            return None
        elif len(subtrees) == 1:
            tree = subtrees[0]
        else:
            tree = PhyloTree(PhyloTreeNode(0))
            for subtree in subtrees:
                tree.root.add_child(subtree.root)
        
        return TrueBMG.bmg_from_tree(tree)
        
        
# --------------------------------------------------------------------------
#                          LRT FROM 2-col. BMG
#                         (new characterization)
# --------------------------------------------------------------------------
            
class TwoColoredLRT:
    
    
    def __init__(self, digraph):
        
        if not isinstance(digraph, nx.DiGraph):
            raise TypeError('not a digraph')
            
        self.digraph = digraph
        self._color_list()
        
        
    def _color_list(self):
        
        self.colors = []
        
        for v in self.digraph.nodes():
            col = self.digraph.nodes[v]['color']
            if col not in self.colors and len(self.colors) < 2:
                self.colors.append(col)
            elif col not in self.colors:
                raise RuntimeError('more than 2 colors in digraph')
        
    
    def _color_count(self, G):
        
        color_count = [0 for _ in range(2)]
        
        for v in G.nodes():
            color_count[self.colors.index(G.nodes[v]['color'])] += 1
        
        return color_count
                
    
    def build(self):
        
        def _build_tree(G):
            
            root_leaves = []
            color_count = self._color_count(G)
            
            for v in G.nodes():
                col_idx = self.colors.index(G.nodes[v]['color'])
                
                if color_count[(col_idx + 1) % 2] == G.out_degree(v):
                    valid = True
                    for v2 in G.predecessors(v):
                        if color_count[col_idx % 2] != G.out_degree(v2):
                            valid = False
                    if valid:
                        root_leaves.append(v)
                    
            if not root_leaves and nx.is_weakly_connected(G):
                return False
                
            node = PhyloTreeNode(0, label='')
            for v in root_leaves:
                node.add_child(PhyloTreeNode(v, label=str(v),
                                             color=G.nodes[v]['color']))
                G.remove_node(v)
            
            for wcc in nx.weakly_connected_components(G):
                child = _build_tree(G.subgraph(wcc).copy())
                
                if not child:
                    return False
                else:
                    node.add_child(child)
                    
            return node
        
        root = _build_tree(self.digraph.copy())
        if not root:
            return False
        else:
            return PhyloTree(root) 

        
        
if __name__ == '__main__':
    
    import time
    
    N = 100
    
    T = PhyloTree.random_colored_tree(N, 2)
    print('--- T ---\n', T.to_newick())
    
    bmg = TrueBMG.bmg_from_tree(T)
    
    start_time1 = time.time()
    lrt_constr = LRTConstructor(bmg, mincut=False)
    lrt1 = lrt_constr.build_tree()
    end_time1 = time.time()
    
    lrt2 = lrt_from_observable_tree(T)
    
    print('--- LRT1 ---\n', lrt1.to_newick())
    print('--- LRT2 ---\n', lrt2.to_newick())
    print('LRTs equal: {}'.format( lrt1.compare_topology(lrt2) ))
    
    bmg = TrueBMG.bmg_from_tree(T)
    
    start_time2 = time.time()
    tc = TwoColoredLRT(bmg)
    lrt3 = tc.build()
    end_time2 = time.time()
    
    lrt4 = lrt_from_observable_tree(T)
    print('--- LRT3 ---\n', lrt3.to_newick())
    print('--- LRT4 ---\n', lrt4.to_newick())
    print('LRTs equal: {}'.format( lrt3.compare_topology(lrt4) ))
    print(end_time1 - start_time1, end_time2 - start_time2)