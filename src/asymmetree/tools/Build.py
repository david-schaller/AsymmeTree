# -*- coding: utf-8 -*-

"""BUILD algorithm (Aho et al. 1981)."""

import networkx as nx

from asymmetree.tools.PhyloTree import PhyloTree, PhyloTreeNode


__author__ = 'David Schaller'


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


class Build:
    
    def __init__(self, R, L, mincut=False, weighted_mincut=False):
        
        self.R = R
        self.L = L
        self.mincut = mincut
        self.weighted_mincut = weighted_mincut
    
    
    def build_tree(self):
        """Build a tree displaying all triples in R if possible."""
        
        self.cut_value = 0
        self.cut_list = []
        
        root = self.aho(self.L, self.R)
        if not root:
            print('no such tree exists')
            return None
        else:
            tree = PhyloTree(root)
            
            # restore leaf IDs and assign unique IDs to all other nodes
            tree.reconstruct_IDs()
            
            return tree
    
    
    def aho(self, L, R):
        """Recursive Aho-algorithm."""
        
        # trivial case: only one leaf left in L
        if len(L) == 1:
            leaf = L.pop()
            return PhyloTreeNode(-1, label=str(leaf))
            
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
                
        subtree_root = PhyloTreeNode(-1, label='')  # place new inner node
        for Ti in child_nodes:
            subtree_root.add_child(Ti)              # add roots of the subtrees to the new node
   
        return subtree_root                         # return the new node
    
    
    def connected_components(self, aho_graph):
        """Determines the connected components of the graph.
        
        And optionally executes a min cut if there is only one component."""
        
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