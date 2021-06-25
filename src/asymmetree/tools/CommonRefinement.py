# -*- coding: utf-8 -*-

"""Fast computation of a common refinement of trees on the same leaf set."""

from collections import deque

import networkx as nx

from asymmetree.datastructures.Tree import Tree, TreeNode


__author__ = 'David Schaller'


class _RefinementConstructor:
    
    def __init__(self, trees):
        
        self.L = None
        
        for T_i in trees:
            
            if not isinstance(Tree, T_i):
                raise TypeError("not a 'Tree' instance")
                
            if not self.L:
                self.L = {v.label for v in T_i.leaves()}
                if len(self.L) == 0:
                    raise RuntimeError('empty tree in tree list')
            else:
                L2 = {v.label for v in T_i.leaves()}
                if self.L != L2:
                    raise RuntimeError('trees must have the same leaf labels')
        
        self.trees = trees
        
        self.V = set()          # vertices in resulting tree
        self.vi_to_v = {}       # inner vertex in input tree --> 
                                # corresponding vertex in resulting tree 
                                
        self.p = [{} for _ in range(len(self.trees))]
                                # v --> lca_Ti (L(Ti(v)))
        
        self.J = {}             # vertex in resulting tree --> indices of
                                # the trees with corresponding vertex
        
        self.Q = deque()        # queue
        
        self.label_dict = {}
        for label in self.L:
            v = TreeNode(len(self.V), label=label)
            self.V.add(v)
            self.Q.append(v)
            self.label_dict[label] = v
            self.J[v] = {i: None for i in range(len(self.trees))}
        
        for i, T_i in enumerate(trees):
            for v_i in T_i.leaves():
                v = self.label_dict[v_i.label]
                self.p[i][v] = v_i
                self.J[v][i] = v_i
            
        
    
    def _leaf_set_cardinalities(self):
        """Compute the number of leaves below each vertex."""
        
        self.l = {}
        
        for T_i in self.trees:
            for v in T_i.postorder():
                if v.is_leaf():
                    self.l[v] = 1
                else:
                    self.l[v] = sum(self.l[w] for w in v.children)
    
    
        