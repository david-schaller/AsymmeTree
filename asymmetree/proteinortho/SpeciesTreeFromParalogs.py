# -*- coding: utf-8 -*-

import itertools

import networkx as nx

from cograph.Cograph import SimpleGraph
from cograph.CographEditor import CographEditor

from simulator.Tree import Tree, TreeNode


class TreeReconstructor:
    
    
    def __init__(self):
        self.R = {}             # (weighted) species triples
        self.L = set()          # set of species leaves
        
        
    def add_ortho_graph(self, ortho_graph):
        
        for cc in nx.connected_components(ortho_graph):
            
            G = SimpleGraph()
            color_dict = {}
            
            for u in cc:
                G.add_node(u)
                for v in ortho_graph.neighbors(u):
                    G.add_edge(u, v)
                u_color = ortho_graph.nodes[u]['color']
                color_dict[u] = u_color
                self.L.add(u_color)
                
            # apply minimal cograph editing (heuristic)
            ce = CographEditor(G)
            cotree = ce.cograph_edit()
            
            # add the informative triples (root is a speciation)
            self._informative_triples(cotree, color_dict)
            
            
    def _informative_triples(self, cotree, color_dict):
        """Add informative species triples from a (co)tree to R."""
        
        cotree.supply_leaves()
        
        for u in cotree.preorder():
            if u.label == "parallel" or not u.children:
                continue
            
            for v1, v2 in itertools.permutations(u.children, 2):
                for x, y in itertools.combinations(v1.leaves, 2):
                    for z in v2.leaves:
                        X, Y, Z = (color_dict[x.ID],
                                   color_dict[y.ID],
                                   color_dict[z.ID])
                        
                        if X != Y and X != Z and Y != Z:
                            self._add_triple(X, Y, Z)
                            
        
    def _add_triple(self, a, b, c, weight=1):
        """Add a triple ab|c (= ba|c)."""
        
        if a <= b:
            triple = (a, b, c)
        else:
            triple = (b, a, c)
            
        if triple in self.R:
            self.R[triple] += weight
        else:
            self.R[triple] = weight
            
            
    def build_species_tree(self):
        """Build a species tree from the informative triples."""
        
        root = self._BUILD(self.L, self.R)
        if not root:
            raise Exception("Could not build a species tree!")
        
        return Tree(root)
    
    
    def _BUILD(self, L, R):
        """Aho's BUILD-algorithm with minimal edge cut."""
        
        if len(L) == 1:                                 # trivial case: only one leaf left in L
            leaf = L.pop()
            return TreeNode(leaf, label=leaf)
        
        aho_graph = self._aho_graph(L, R)               # construct the Aho-graph
                                                        # determine connected components A1, ..., Ak
        conn_comps = self._connected_comp(aho_graph)
        
        if len(conn_comps) <= 1:                        # return False if less than 2 connected components
            print("Connected component:\n", conn_comps)
            print(R)
            return False
        
        child_nodes = []
        for cc in conn_comps:
            Li = set(cc)                                # list/dictionary --> set
            Ri = {}
            for t, weight in R.items():                 # determine which triples are in the subtree
                if Li.issuperset(t):
                    Ri[t] = weight
            Ti = self._BUILD(Li, Ri)                    # recursive call
            if not Ti:
                return False                            # raise False to previous call of aho()
            else:
                child_nodes.append(Ti)
                
        subtree_root = TreeNode(0, label="")            # place new inner node
        for Ti in child_nodes:
            Ti.parent = subtree_root                    # set the parent
            subtree_root.children.append(Ti)            # add roots of the subtrees to the new node
        return subtree_root                             # return the new node
    
    
    def _aho_graph(self, L, R):
        """Auxiliary graph for Aho-algorithm.
    
        Takes a list of leaves L and a list of tuples R of the form (xy|z) and
        builds a graph with V = L and E = { {x,y} | (xy|z) in R }.
        Applies Minimal Edge Cut if necessary.
        """
        G = nx.Graph()
        G.add_nodes_from(L)

        for a, b, c in R:
            if G.has_edge(a, b):
                G[a][b]["weight"] += R[(a,b,c)]
            else:
                G.add_edge(a, b, weight=R[(a,b,c)])
                
        return G
                
                
    def _connected_comp(self, G):
#        return list(nx.connected_components(G))
        if nx.number_connected_components(G) > 1:
            return list(nx.connected_components(G))
        else:
            cut = nx.stoer_wagner(G)                    # Stoer–Wagner algorithm
                                                        # for minimal weighted edge cut
            return cut[1]