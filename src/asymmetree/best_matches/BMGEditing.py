# -*- coding: utf-8 -*-

"""Heuristics for BMG editing."""

import itertools

import networkx as nx

from asymmetree.datastructures.PhyloTree import PhyloTree, PhyloTreeNode

import asymmetree.best_matches.LeastResolvedTree as LRT
from asymmetree.best_matches.TrueBMG import bmg_from_tree
from asymmetree.tools.GraphTools import sort_by_colors
from asymmetree.tools.TreeTools import LCA
from asymmetree.tools.Build import (aho_graph,
                                    Build,
                                    best_pair_merge_first)
from asymmetree.tools.Partitioning import Karger


__author__ = 'David Schaller'



class BMGEditor:
    
    def __init__(self, G):
        
        self.G = G
        self.color_dict = sort_by_colors(G)
        
        self.L = [v for v in G.nodes()]
        
        # informative and forbidden triples
        self.R, self.F = LRT.informative_forbidden_triples(G,
                             color_dict=self.color_dict)
        
        # current tree built by one of the heuristics
        self._tree = None
        
    
    def extract_consistent_triples(self):
        
        if not self._tree:
            raise RuntimeError('no tree has been built yet')
        
        lca = LCA(self._tree)
        
        return lca.consistent_triples(self.R)
    
    
    def build_mincut(self, weighted=True):
        
        build = Build(self.R, self.L,
                      mincut=True, weighted_mincut=weighted)
        self._tree = build.build_tree()
        
    
    def build_bpmf(self):
        
        self._tree = best_pair_merge_first(self.R, self.L,
                                           triple_weights=None)
        
    def build_karger(self):
        
        build = BuildKarger(self.R, self.L, self.G)
        self._tree = build.build_tree()
    
    
    def get_bmg(self, extract_triples_first=False):
        
        if not extract_triples_first:
            self._tree.reconstruct_info_from_graph(self.G)
            return bmg_from_tree(self._tree)
        
        else:
            R_consistent = self.extract_consistent_triples()
            build = Build(R_consistent, self.L, mincut=False)
            tree = build.build_tree()
            tree.reconstruct_info_from_graph(self.G)
            return bmg_from_tree(tree)
        
        
def unsatisfiability_cost(partition, G):
    """Return the unsatisfiability cost of a partition.
    
    This cost is the number of unsatisfiable (non-)arcs induced by the
    bipartition.
    """
    
    color_sets = [{} for V in partition]
    
    for V, colors in zip(partition, color_sets):
        for v in V:
            c = G.nodes[v]['color']
            if c not in colors:
                colors[c] = 1
            else:
                colors[c] += 1
    
    cost = 0
    
    for V, colors in zip(partition, color_sets):
        
        if not isinstance(V, set):
            V = set(V)
        
        for x, y in itertools.product(V, G.nodes()):
            
            y_color = G.nodes[y]['color']
            
            # skip pairs with the same color
            if G.nodes[x]['color'] == y_color:
                continue
            
            if y not in V:
                if G.has_edge(x, y):
                    if colors.get(y_color):
                        cost += 1
                elif not colors.get(y_color):
                    cost += 1
            elif not G.has_edge(x, y) and colors.get(y_color) == 1:
                cost += 1
                    
    return cost


class BuildKarger:
    """BUILD algorithm with minimal cost bipartition."""
    
    def __init__(self, R, L, digraph):
        
        self.R = R
        self.L = L
        
        # colored digraph
        self.G = digraph
    
    
    def build_tree(self, return_root=False):
        """Build a tree displaying all triples in R if possible.
        
        Keyword arguments:
            return_root - if True, return 'PhyloTreeNode' instead of
                'PhyloTree' instance
        """
        
        root = self._aho(self.L, self.R)
        return root if return_root else PhyloTree(root)
    
    
    def _aho(self, L, R):
        """Recursive Aho-algorithm."""
        
        # trivial case: only one leaf left in L
        if len(L) == 1:
            leaf = L.pop()
            return PhyloTreeNode(leaf, label=str(leaf))
            
        help_graph = aho_graph(R, L, weighted=False)
        conn_comps = self._connected_components(help_graph)
        
        # there should always be >= 2 components
        if len(conn_comps) < 2:
            raise RuntimeError('only one component')
        
        child_nodes = []
        for cc in conn_comps:
            Li = set(cc)                    # list/dictionary --> set
            Ri = []
            for t in R:                     # construct triple subset
                if Li.issuperset(t):
                    Ri.append(t)
            Ti = self._aho(Li, Ri)          # recursive call
            if not Ti:
                return False                # raise False to previous call
            else:
                child_nodes.append(Ti)
                
        node = PhyloTreeNode(-1)            # place new inner node
        for Ti in child_nodes:
            node.add_child(Ti)              # add roots of the subtrees
   
        return node
    
    
    def _connected_components(self, aho_graph):
        """Determines the connected components of the graph.
        
        And optionally executes a min cut if there is only one component."""
        
        conn_comps = list(nx.connected_components(aho_graph))
        if len(conn_comps) > 1:
            return conn_comps
        else:
            karger = Karger(aho_graph)
            
            best_cost, best_V1, best_V2 = None, None, None
            
            for V1, V2, cutvalue in karger.generate():
                cost = unsatisfiability_cost([V1, V2], self.G)
                
                if best_cost is None or cost < best_cost:
                    best_cost, best_V1, best_V2 = cost, V1, V2
            
            return best_V1, best_V2


if __name__ == '__main__':
    
    from asymmetree.tools.GraphTools import disturb_graph, symmetric_diff
    
    random_tree = PhyloTree.random_colored_tree(20, 4,
                                                force_all_colors=True)
    bmg = bmg_from_tree(random_tree)
    disturbed = disturb_graph(bmg, 0.1, 0.1, preserve_properly_colored=True)
    
    editor = BMGEditor(disturbed)
    
    editor.build_mincut()
    bmg1 = editor.get_bmg(extract_triples_first=False)
    bmg2 = editor.get_bmg(extract_triples_first=True)
    
    editor.build_bpmf()
    bmg3 = editor.get_bmg(extract_triples_first=False)
    bmg4 = editor.get_bmg(extract_triples_first=True)
    
    editor.build_karger()
    bmg5 = editor.get_bmg(extract_triples_first=False)
    bmg6 = editor.get_bmg(extract_triples_first=True)
    
    print('-------------')
    print([(v, bmg.nodes[v]['color']) for v in bmg.nodes()])
    print('orginal', bmg.edges(), bmg.size())
    print('disturbed', disturbed.edges(), disturbed.size())
    print('BMG1', bmg1.edges(), bmg1.size())
    print('BMG2', bmg2.edges(), bmg2.size())
    print('BMG3', bmg3.edges(), bmg3.size())
    print('BMG4', bmg4.edges(), bmg4.size())
    print('BMG5', bmg5.edges(), bmg5.size())
    print('BMG6', bmg6.edges(), bmg6.size())
    print('-------------')
    print('original BMG is BMG:   ', bool(LRT.is_bmg(bmg)  ) )
    print('disturbed graph is BMG:', bool(LRT.is_bmg(disturbed)) )
    print('original vs disturbed:', symmetric_diff(bmg, disturbed))
    print('-------------')
    
    print('orig./dist. vs BMG1 (Mincut, no extr):',
          symmetric_diff(bmg, bmg1), symmetric_diff(disturbed, bmg1))
    print('orig./dist. vs BMG2 (Mincut, extr):',
          symmetric_diff(bmg, bmg2), symmetric_diff(disturbed, bmg2))
    print('orig./dist.vs BMG3 (BPMF, no extr):',
          symmetric_diff(bmg, bmg3), symmetric_diff(disturbed, bmg3))
    print('orig./dist. vs BMG4 (BPMF, extr):',
          symmetric_diff(bmg, bmg4), symmetric_diff(disturbed, bmg4))
    print('orig./dist.vs BMG5 (Karger min cost, no extr):',
          symmetric_diff(bmg, bmg5), symmetric_diff(disturbed, bmg5))
    print('orig./dist. vs BMG6 (Karger min cost, extr):',
          symmetric_diff(bmg, bmg6), symmetric_diff(disturbed, bmg6))
    
