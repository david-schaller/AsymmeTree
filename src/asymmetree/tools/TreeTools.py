# -*- coding: utf-8 -*-

"""
Tree tools.
"""

import itertools

from asymmetree.datastructures.Tree import Tree, TreeNode


__author__ = 'David Schaller'


# based on
# - Bender, M. A., M. Farach-Colton, G. Pemmasani, S. Skiena, and P. Sumazin:
#   Lowest common ancestors in trees and directed acyclic graphs.
#   In: Journal of Algorithms. 57, Nr. 2, November 2005, S. 75–94.
#   ISSN 0196-6774. doi:10.1016/j.jalgor.2005.08.001.
# - https://cp-algorithms.com/data_structures/sparse-table.html
class LCA:
    """Compute last common ancestors efficiently.
    
    Uses a reduction to the Range minimum query (RMQ) problem and a sparse
    table implementation.
    Preprocessing complexity: O(n * log n)
    Query complexity: O(1)
    where n is the number of vertices in the tree.
    """
    
    def __init__(self, tree):
        
        if not isinstance(tree, Tree):
            raise TypeError("tree must be of type 'Tree'")
            
        self.tree = tree
        
        self.V = [v for v in self.tree.preorder()]
        self.index = {v: i for i, v in enumerate(self.V)}
        
        # store leaf IDs for rooted triple queries
        self.leaf_ids = {v.ID: v for v in self.V if v.is_leaf()}
        
        self.euler_tour = []
        # levels of the vertices in the Euler tour
        self.L = []
        
        for v, level in self.tree.euler_and_level():
            self.euler_tour.append(self.index[v])
            self.L.append(level)
        
        # repres. of vertices in the Euler tour (index of first occurence)
        self.R = [None for _ in range(len(self.V))]
        for j, i in enumerate(self.euler_tour):
            if self.R[i] is None:
                self.R[i] = j
                
        # build sparse table for range minimum query (RMQ)
        self._precompute_logs()
        self._RMQ_sparse_table()
        
    
    def get(self, v1, v2):
        """Return the last common ancestor of two nodes."""
        
        if v1 is v2:
            return v1
        
        r1 = self.R[self.index[v1]]
        r2 = self.R[self.index[v2]]
        if r1 > r2:
            r1, r2 = r2, r1
        return self.V[ self.euler_tour[self._RMQ_query(r1, r2)] ]
    
    
    def displays_triple(self, a, b, c):
        """Return whether the tree displays the rooted triple ab|c (= ba|c)."""
        
        abc = [a, b, c]
        
        for i in range(3):
            if not isinstance(abc[i], TreeNode):
                if abc[i] not in self.leaf_ids:
                    return False
                else:
                    abc[i] = self.leaf_ids[abc[i]]
            elif abc[i] not in self.index:
                return False
        
        return self._has_triple(*abc)
    
    
    def consistent_triples(self, triples):
        """Return a list with the subset of 'triples' that are displayed."""
        
        return [t for t in triples if self.displays_triple(*t)]
    
    
    def consistent_triple_generator(self, triples):
        """Generator for the items in 'triples' that are displayed."""
        
        for t in triples:
            if self.displays_triple(*t):
                yield t
    
    
    def _precompute_logs(self):
        
        n = len(self.L)
        self.log = [0 for _ in range(n + 1)]
        for i in range(2, n + 1):
            self.log[i] = int(self.log[i//2]) + 1
        
        
    def _RMQ_sparse_table(self):
        
        n = len(self.L)
        K = self.log[n]
        
        # lookup table M
        self.M = [[0 for j in range(K + 1)] for i in range(n)]
        
        # initialize the intervals with length 1
        for i in range(n):
            self.M[i][0] = i
         
        # dynamic programming: compute values from smaller to bigger intervals  
        for j in range(1, K + 1):
            
            # compute minimum value for all intervals with size 2^j
            for i in range(n - (1 << j) + 1):
                
                if (self.L[ self.M[i][j - 1] ] <
                    self.L[ self.M[i + (1 << (j - 1))][j - 1] ]):
                    self.M[i][j] = self.M[i][j - 1]
                else:
                    self.M[i][j] = self.M[i + (1 << (j - 1))][j - 1]
    
    
    def _RMQ_query(self, i, j):  
        
        k = self.log[j - i + 1]
        if self.L[self.M[i][k]] < self.L[self.M[j - (1 << k) + 1][k]]:
            return self.M[i][k]
        else:
            return self.M[j - (1 << k) + 1][k]
        
    
    def _has_triple(self, a, b, c):
        
        lca_ab = self.get(a, b)
        return lca_ab is not self.get(lca_ab, c)
        

def lcas_naive(tree):
    """Naive computation of the last common ancestor for all leaf pairs."""
    
    lcas = {}
    tree.supply_leaves()
    
    for u in tree.preorder():
        for v1, v2 in itertools.permutations(u.children, 2):
            for x, y in itertools.product(v1.leaves, v2.leaves):
                lcas[x, y] = u
    
    return lcas


if __name__ == '__main__':
    
    import random
    
    # import time
    # tree = Tree.random_tree(4000)
    # leaves = tree.supply_leaves()
    
    # start_time1 = time.time()
    # lca_fast = LCA(tree)
    # end_time1 = time.time()
    
    # start_time2 = time.time()
    # lca_naive = lcas_naive(tree)
    # end_time2 = time.time()
    
    # print(end_time1 - start_time1, end_time2 - start_time2)
    
    # for i in range(10):
    #     v1, v2 = random.sample(leaves, 2)
        
    #     print(v1, v2, lca_fast.get(v1, v2), lca_naive[v1, v2])
    
    tree = Tree.random_tree(10)
    print(tree.to_newick())
    leaves = tree.supply_leaves()
    lca = LCA(tree)
    
    for i in range(10):
        a, b, c = random.sample(leaves, 3)
        
        print(a, b, c, lca.displays_triple(a, b, c),
                       lca.displays_triple(a.ID, b.ID, c.ID))
    