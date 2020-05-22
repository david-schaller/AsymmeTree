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


__author__ = "David Schaller"


# --------------------------------------------------------------------------
#                      LRT FROM OBSERVABLE GENE TREE
# --------------------------------------------------------------------------


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


def LRT_from_observable_tree(T):
    """Computes the Least Resolved Tree from a tree.
    
    The unique Least Resolved Tree from a leaf-colored (observable)
    gene tree is computed by contraction of all redundant edges.
    """
    
    LRT = T.copy()
    if not LRT.root:
        return LRT
    
    # remove planted root if existent
    LRT.remove_planted_root()
    
    # assign list of leaves to each node
    LRT.supply_leaves()
    
    subtree_colors = {}
    for v in LRT.preorder():
        subtree_colors[v] = {leaf.color for leaf in v.leaves}
        
    arc_colors = _arc_colors(LRT, subtree_colors)
    redundant_edges = _redundant_edges(LRT, subtree_colors, arc_colors)
    LRT.contract(redundant_edges)
    LRT = LRT.topology_only()
    
    return LRT


# --------------------------------------------------------------------------
#                               LRT FROM BMG
# --------------------------------------------------------------------------

class LRTConstructor:
    
    def __init__(self, G, mincut=True):
        
        self.G = G
        self.R = []
        self.L = set()
        self.color_dict = {}
        self.mincut = mincut
        self.cut_value = 0
        self.cut_list = []
        
    
    def initialize(self):
        """Initialize color dictionary and leaf set."""
        
        for v in self.G.nodes:
            color = self.G.nodes[v]["color"]
            if color not in self.color_dict:
                self.color_dict[color] = [v]
            else:
                self.color_dict[color].append(v)
            self.L.add(v)
    
    
    def informative_triples(self):
        """Compute the informative triples."""
        
        for col1, col2 in itertools.permutations(self.color_dict.keys(), 2):
            list1, list2 = self.color_dict[col1], self.color_dict[col2]
            for a in list1:
                for b, c in itertools.permutations(list2, 2):
                    if self.G.has_edge(a,b) and (not self.G.has_edge(a,c)):   # X1 / X2 / X3 / X4
                        self.R.append( (a,b,c) )
    
    
    def build_tree(self):
        """Build the least resolved tree from informative triples."""
        
        self.initialize()
        self.informative_triples()
        
        root = self.aho(self.L, self.R)
        if not root:
            print("No such tree exists!")
            return None
        else:
            return PhyloTree(root)
    
    
    def aho(self, L, R):
        """Recursive Aho-algorithm."""
        
        if len(L) == 1:                                 # trivial case: only one leaf left in L
            leaf = L.pop()
            return PhyloTreeNode(leaf, label=self.G.nodes[leaf]["label"],
                            color=self.G.nodes[leaf]["color"])
            
        help_graph_A = HelpGraph(L, R, self)            # construct the help graph A(R)
                                                        # determine connected components A1, ..., Ak
        conn_comps = help_graph_A.connected_comp(mincut=self.mincut)
        
        if len(conn_comps) <= 1:                        # return False if less than 2 connected components
            print("Connected component:\n", conn_comps)
            return False
        child_nodes = []
        for cc in conn_comps:
            Li = set(cc)                                # list/dictionary --> set
            Ri = []
            for t in R:                                 # determine which triples are in the subtree
                if Li.issuperset(t):
                    Ri.append(t)
            Ti = self.aho(Li, Ri)                       # recursive call
            if not Ti:
                return False                            # raise False to previous call of aho()
            else:
                child_nodes.append(Ti)
                
        subtree_root = PhyloTreeNode(0, label="")       # place new inner node
        for Ti in child_nodes:
            subtree_root.add_child(Ti)                  # add roots of the subtrees to the new node
   
        return subtree_root                             # return the new node
    
    
    @staticmethod
    def correct_BMG(BMG_original):
        """Build the LRT (using Mincut in Aho algorithm and return its BMG and RBMG."""
        
        subtrees = []
        for sg in (BMG_original.subgraph(c) for c in nx.weakly_connected_components(BMG_original)):
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
        
        return TrueBMG.BMG_from_tree(tree)


class HelpGraph:
    """Help graph for Aho-algorithm.
    
    Takes a list of leaves L and a list of tuples R of the form (x,y|z) and
    builds a graph with V = L and E = { (x,y) | (x,y|z) in R }.
    Provides optional Minimal Edge Cut for function connected_comp().
    """
    
    def __init__(self, L, R, aho):
        """Constructs the help graph as NetworkX graph.
        
        Edges (a,b) are weighted by the number of occurrences in triples
        of the form ab|x or ba|x.
        """
        
        self.G = nx.Graph()
        self.G.add_nodes_from(L)
        self.aho = aho

        for t1, t2, t3 in R:
            if not self.G.has_edge(t1, t2):
                self.G.add_edge(t1, t2, weight=1)
#            if self.G.has_edge(t1, t2):
#                self.G[t1][t2]["weight"] += 1
#            else:
#                self.G.add_edge(t1, t2, weight=1)
    
    
    def connected_comp(self, mincut=False):
        """Determines the connected components of the graph.
        
        Keyword arguments:
            mincut -- if True and graph has only one component, a 
                      minimal edge cut is applied first, function then
                      always return > 1 components (default=False).
        """
        
        if (not mincut) or nx.number_connected_components(self.G) > 1:
            return list(nx.connected_components(self.G))
        else:
            cut = nx.stoer_wagner(self.G)           # Stoer–Wagner algorithm
            self.aho.cut_value += cut[0]
            if len(cut[1][0]) < len(cut[1][1]):
                smaller_comp = cut[1][0]
            else:
                smaller_comp = cut[1][1]
            for edge in self.G.edges:
                if ((edge[0] in smaller_comp and edge[1] not in smaller_comp) or
                    (edge[1] in smaller_comp and edge[0] not in smaller_comp)):
                    self.aho.cut_list.append(edge)
            return cut[1]
        
        
if __name__ == '__main__':
    
    T = PhyloTree.random_colored_tree(30, 5)
    print('--- T ---\n', T.to_newick())
    
    BMG = TrueBMG.BMG_from_tree(T)
    
    lrt_constr = LRTConstructor(BMG, mincut=False)
    LRT1 = lrt_constr.build_tree()
    
    LRT2 = LRT_from_observable_tree(T)
    
    print('--- LRT1 ---\n', LRT1.to_newick())
    print('--- LRT2 ---\n', LRT2.to_newick())
    
    print('LRTs equal: {}'.format( LRT1.compare_topology(LRT2) ))