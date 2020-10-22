# -*- coding: utf-8 -*-

"""
Algorithms for graph (bi)partitioning.
"""

import random, math, itertools
import numpy as np
import networkx as nx

from asymmetree.datastructures.LinkedList import LinkedList
from asymmetree.datastructures.AVLTree import TreeSet


__author__ = 'David Schaller'

# ----------------------------------------------------------------------------
#                           Karger's algorithm
#
# Implementation of Karger's algorithm for minimal edge cut. Randomized
# algorithm that returns a minimal edge cut with high probability.
# This implementation is based on: 
# http://web.stanford.edu/class/archive/cs/cs161/cs161.1166/lectures/
# lecture15.pdf
# ----------------------------------------------------------------------------


class _Supernode:
    
    def __init__(self, initial_node):
        
        self.nodes = LinkedList()
        self.nodes.append(initial_node)
        
    
    def __len__(self):
        
        return len(self.nodes)
    
    
    def __lt__(self, other):
        
        return id(self) < id(other)
    
    
    def merge_from(self, other):
        
        self.nodes.concatenate(other.nodes)
        

def _sort2(a, b):
    
    return (a, b) if a < b else (b, a)
        
        
class Karger():
    
    def __init__(self, graph):
        
        if not isinstance(graph, nx.Graph):
            raise TypeError('must be a NetworkX graph')
            
        if not graph.order() >= 2:
            raise ValueError('graph must have >= 2 nodes')
            
        if not nx.is_connected(graph):
            raise ValueError('graph not connected')
            
        self.graph = graph
            
    
    def run(self):
        """Execute a single run of the contraction procedure."""
        
        return self._karger_main()
    
    
    def best_out_of(self, runs=None):
        """Execute multiple runs and return the one with the smallest cutvalue.
        
        Keyword arguments:
            runs -- number of runs; default is None, in which case
                O(n^2 * log n) runs are executed
        """
        
        if not runs:
            runs = self._default_run_number()
        
        best_V1, best_V2, best_cutvalue = None, None, None
        
        for _ in range(runs):
            V1, V2, cutvalue = self._karger_main()
            
            if not best_cutvalue or cutvalue < best_cutvalue:
                best_V1, best_V2, best_cutvalue = V1, V2, cutvalue
                
        return best_V1, best_V2, best_cutvalue
    
    
    def generate(self, runs=None):
        """Generator for runs of the contraction procedure.
        
        Keyword arguments:
            runs -- number of runs; default is None, in which case
                O(n^2 * log n) runs are generated
        """
        
        if not runs:
            runs = self._default_run_number()
        
        for _ in range(runs):
            yield self._karger_main()
            
            
    def _default_run_number(self):
        
        n = self.graph.order()
        return round(n**2 * math.log(n))
    
        
    def _full_edge_set(self):
        
        # construct the balanced search tree for the full edge set once
        if not hasattr(self, 'full_egde_set'):
            self.full_egde_set = TreeSet()
            for u, v in self.graph.edges():
                self.full_egde_set.add( (u, v) )
        
        # otherwise copy which is O(|E|)
        return self.full_egde_set.copy()
        
        
    def _initialize(self):
        
        self.supernodes = set()
        self.supernode_dict = {}            # maps node --> supernode
        self.superedges = {}                # maps superedge --> edges
        self.edges = self._full_edge_set()  # set of all (non-loop) edges
        
        for v in self.graph.nodes():
            sn = _Supernode(v)
            self.supernode_dict[v] = sn
            self.supernodes.add(sn)
        
        for u, v in self.graph.edges():
            u_sn = self.supernode_dict[u]
            v_sn = self.supernode_dict[v]
            e = _sort2(u_sn, v_sn)
            self.superedges[e] = LinkedList( [(u, v)] )
            
            
    def _merge(self, a, b):
        """Merge two supernodes a and b."""
        
        # bigger supernode becomes the new supernode x
        if len(a) >= len(b):
            x = a
            other = b
        else:
            x = b
            other = a
            
        for v in other.nodes:
            self.supernode_dict[v] = x
        x.merge_from(other)
        self.supernodes.remove(other)
        
        for d in self.supernodes:
            
            if d is x:
                continue
            
            # merge the edge lists of superedges ad and bd
            edges_xd = LinkedList()
            for se in (_sort2(a, d), _sort2(b, d)):
                if se in self.superedges:
                    edges_xd.concatenate(self.superedges[se])
                    del self.superedges[se]
            if edges_xd:
                self.superedges[_sort2(x, d)] = edges_xd
                
                
    def _karger_main(self):
        
        self._initialize()
        
        while len(self.supernodes) > 2:
            
            u, v = random.choice(self.edges)
            u_sn = self.supernode_dict[u]
            v_sn = self.supernode_dict[v]
            
            self._merge(u_sn, v_sn)
            
            se = _sort2(u_sn, v_sn)
            self.edges.difference_update(self.superedges[se])
            del self.superedges[se]
        
        x = self.supernodes.pop()
        y = self.supernodes.pop()

        cut_value = len(self.superedges[_sort2(x, y)])
        # print('cut_edges', list(self.superedges[_sort2(x, y)]))
            
        return cut_value, [list(x.nodes), list(y.nodes)]
    

# ----------------------------------------------------------------------------
#                           Greedy bipartition
# ----------------------------------------------------------------------------

def partition_cut_value(partition, G):
    
    cut_value = 0
    
    for Vi, Vj in itertools.combinations(partition, 2):
        for x, y in itertools.product(Vi, Vj):
            if G.has_edge(x, y):
                cut_value += 1
    
    return cut_value
    

def _move(from_set, to_set, x):
    
    from_set.remove(x)
    to_set.add(x)


def greedy_bipartition(V, f_cost, args=()):
    """Randomized greedy bipartitioning with custom cost function."""
    
    if len(V) < 2:
        raise ValueError('iterable must have >=2 elements')
    
    V1, V2 = set(), set(V)
    remaining = set(V)
    
    best_cost, best_bp = float('inf'), None
    
    for _ in range(len(V)-1):
        
        best_cost_local, best_x = float('inf'), None
        for x in remaining:
            _move(V2, V1, x)
            cost = f_cost([V1, V2], *args)
            if cost < best_cost_local:
                best_cost_local = cost
                best_x = [x]
            elif cost == best_cost_local:
                best_x.append(x) 
            _move(V1, V2, x)
            
        x = random.choice(best_x)
        _move(V2, V1, x)
        remaining.remove(x)
        
        if best_cost_local < best_cost:
            best_cost = best_cost_local
            best_bp = [(list(V1), list(V2))]
        elif best_cost_local == best_cost:
            best_bp.append( (list(V1), list(V2)) )
            
    return best_cost, random.choice(best_bp)


# ----------------------------------------------------------------------------
#                        Adaptive walk bipartition
# ----------------------------------------------------------------------------

def gradient_walk_bipartition(V, f_cost, args=(), initial_bipartition=None):
    
    if len(V) < 2:
        raise ValueError('iterable must have >=2 elements')
    
    # partition P
    P = [set(), set()]
    lookup = {}
    
    if not isinstance(V, (list, tuple)):
        V = list(V)
    
    if not initial_bipartition:
        for k, m in enumerate(np.random.permutation(len(V))):
            x = V[m]
            i = 0 if k <= len(V) / 2 else 1
            P[i].add(x)
            lookup[x] = i
    elif len(initial_bipartition) == 2:
        for i in range(2):
            for x in initial_bipartition[i]:
                P[i].add(x)
                lookup[x] = i
    else:
        raise ValueError('not a bipartition')
        
    best_cost = f_cost([*P], *args)
    best_bp = tuple(list(s) for s in P)
    
    while True:
        
        best_cost_local, best_x = float('inf'), None
        for x in V:
            V1 = P[lookup[x]]
            V2 = P[(lookup[x] + 1) % 2]
            
            # avoid that one set becomes empty
            if len(V1) == 1:
                continue
            
            _move(V1, V2, x)
            cost = f_cost([V1, V2], *args)
            if cost < best_cost_local:
                best_cost_local = cost
                best_x = [x]
            elif cost == best_cost_local:
                best_x.append(x) 
            _move(V2, V1, x)
        
        if best_cost_local < best_cost:
            x = random.choice(best_x)
            V1 = P[lookup[x]]
            V2 = P[(lookup[x] + 1) % 2]
            _move(V1, V2, x)
            lookup[x] = (lookup[x] + 1) % 2
            
            best_cost = best_cost_local
            best_bp = tuple(list(s) for s in P)
        else:
            break
    
    return best_cost, best_bp
        
    

if __name__ == '__main__':
    
    G = nx.Graph()
    G.add_edges_from([('a', 'b'), ('a', 'c'), ('a', 'd'),
                      ('b', 'd'), ('c', 'd'), ('c', 'e'), ('d', 'e')])
    
    karger = Karger(G)
    for result in karger.generate(3):
        print(*result)
        
    V = list(G.nodes())
    
    print('---------')
    for i in range(3):
        print( *greedy_bipartition(V, partition_cut_value, args=(G,)) )
    
    print('---------')
    for i in range(3):
        print( *gradient_walk_bipartition(V, partition_cut_value, args=(G,)) )
        
