# -*- coding: utf-8 -*-

"""Heuristics for BMG editing."""

import itertools

# import networkx as nx

from asymmetree.datastructures.PhyloTree import PhyloTree#, PhyloTreeNode

import asymmetree.best_matches.LeastResolvedTree as LRT
from asymmetree.best_matches.TrueBMG import bmg_from_tree
from asymmetree.tools.GraphTools import sort_by_colors
from asymmetree.datastructures.Tree import LCA
from asymmetree.tools.Build import Build2, best_pair_merge_first


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
        
        build = Build2(self.R, self.L,
                       allow_inconsistency=True,
                       bipart_method='mincut',
                       weighted_mincut=weighted)
        self._tree = build.build_tree()
        
    
    def build_bpmf(self):
        
        self._tree = best_pair_merge_first(self.R, self.L,
                                           triple_weights=None)
        
    def build_karger(self):
        
        build = Build2(self.R, self.L,
                       allow_inconsistency=True,
                       bipart_method='karger',
                       cost_function=unsatisfiability_cost,
                       cost_function_args=(self.G,),)
        self._tree = build.build_tree()
    
    
    def build_greedy(self):
        
        build = Build2(self.R, self.L,
                       allow_inconsistency=True,
                       bipart_method='greedy',
                       cost_function=unsatisfiability_cost,
                       cost_function_args=(self.G,),)
        self._tree = build.build_tree()
    
    
    def build_gradient_walk(self):
        
        build = Build2(self.R, self.L,
                       allow_inconsistency=True,
                       bipart_method='gradient_walk',
                       cost_function=unsatisfiability_cost,
                       cost_function_args=(self.G,),)
        self._tree = build.build_tree()
    
    
    def get_bmg(self, extract_triples_first=False):
        
        if not extract_triples_first:
            self._tree.reconstruct_info_from_graph(self.G)
            return bmg_from_tree(self._tree)
        
        else:
            R_consistent = self.extract_consistent_triples()
            build = Build2(R_consistent, self.L, allow_inconsistency=False)
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


if __name__ == '__main__':
    
    from time import time
    
    from asymmetree.tools.GraphTools import disturb_graph, symmetric_diff
    
    random_tree = PhyloTree.random_colored_tree(20, 4,
                                                force_all_colors=True)
    bmg = bmg_from_tree(random_tree)
    disturbed = disturb_graph(bmg, 0.1, 0.1, preserve_properly_colored=True)
    
    editor = BMGEditor(disturbed)
    
    start_time = time()
    editor.build_mincut()
    time_mincut = time() - start_time
    bmg1 = editor.get_bmg(extract_triples_first=False)
    bmg2 = editor.get_bmg(extract_triples_first=True)
    
    start_time = time()
    editor.build_bpmf()
    time_bpmf = time() - start_time
    bmg3 = editor.get_bmg(extract_triples_first=False)
    bmg4 = editor.get_bmg(extract_triples_first=True)
    
    start_time = time()
    editor.build_karger()
    time_karger = time() - start_time
    bmg5 = editor.get_bmg(extract_triples_first=False)
    bmg6 = editor.get_bmg(extract_triples_first=True)
    
    start_time = time()
    editor.build_greedy()
    time_greedy = time() - start_time
    bmg7 = editor.get_bmg(extract_triples_first=False)
    bmg8 = editor.get_bmg(extract_triples_first=True)
    
    start_time = time()
    editor.build_gradient_walk()
    time_gradw = time() - start_time
    bmg9 = editor.get_bmg(extract_triples_first=False)
    bmg10 = editor.get_bmg(extract_triples_first=True)
    
    print('-------------')
    print([(v, bmg.nodes[v]['color']) for v in bmg.nodes()])
    print('orginal', bmg.edges(), bmg.size())
    print('disturbed', disturbed.edges(), disturbed.size())
    # print('BMG1', bmg1.edges(), bmg1.size())
    # print('BMG2', bmg2.edges(), bmg2.size())
    # print('BMG3', bmg3.edges(), bmg3.size())
    # print('BMG4', bmg4.edges(), bmg4.size())
    # print('BMG5', bmg5.edges(), bmg5.size())
    # print('BMG6', bmg6.edges(), bmg6.size())
    # print('BMG7', bmg7.edges(), bmg7.size())
    # print('BMG8', bmg8.edges(), bmg8.size())
    # print('BMG9', bmg9.edges(), bmg9.size())
    # print('BMG8', bmg10.edges(), bmg10.size())
    print('-------------')
    print('original BMG is BMG:   ', bool(LRT.is_bmg(bmg)  ) )
    print('disturbed graph is BMG:', bool(LRT.is_bmg(disturbed)) )
    print('original vs disturbed:', symmetric_diff(bmg, disturbed))
    print('-------------')
    
    print('Mincut time:', time_mincut)
    print('orig./dist. vs BMG1 (Mincut, no extr):',
          symmetric_diff(bmg, bmg1), symmetric_diff(disturbed, bmg1))
    print('orig./dist. vs BMG2 (Mincut, extr):',
          symmetric_diff(bmg, bmg2), symmetric_diff(disturbed, bmg2))
    
    print('BPMF time:', time_bpmf)
    print('orig./dist.vs BMG3 (BPMF, no extr):',
          symmetric_diff(bmg, bmg3), symmetric_diff(disturbed, bmg3))
    print('orig./dist. vs BMG4 (BPMF, extr):',
          symmetric_diff(bmg, bmg4), symmetric_diff(disturbed, bmg4))
    
    print('Karger time:', time_karger)
    print('orig./dist.vs BMG5 (Karger min cost, no extr):',
          symmetric_diff(bmg, bmg5), symmetric_diff(disturbed, bmg5))
    print('orig./dist. vs BMG6 (Karger min cost, extr):',
          symmetric_diff(bmg, bmg6), symmetric_diff(disturbed, bmg6))
    
    print('Greedy time:', time_greedy)
    print('orig./dist.vs BMG7 (Greedy min cost, no extr):',
          symmetric_diff(bmg, bmg7), symmetric_diff(disturbed, bmg7))
    print('orig./dist. vs BMG8 (Greedy min cost, extr):',
          symmetric_diff(bmg, bmg8), symmetric_diff(disturbed, bmg8))
    
    print('Grad. W. time:', time_gradw)
    print('orig./dist.vs BMG9 (Grad. W. min cost, no extr):',
          symmetric_diff(bmg, bmg9), symmetric_diff(disturbed, bmg9))
    print('orig./dist. vs BMG10 (Grad. W. min cost, extr):',
          symmetric_diff(bmg, bmg10), symmetric_diff(disturbed, bmg10))
    
