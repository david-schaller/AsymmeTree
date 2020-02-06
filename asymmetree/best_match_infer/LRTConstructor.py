# -*- coding: utf-8 -*-

"""
Construction of Least Resolved Trees (LRT).

Computes an LRT from a not necessarily valid BMG.
"""

import itertools

import networkx as nx

from best_match_infer import TrueBMG
from tools.PhyloTree import PhyloTree, PhyloTreeNode


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


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
                
        subtree_root = PhyloTreeNode(0, label="inner")  # place new inner node
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
        
        return TrueBMG.best_match_graphs(tree)


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