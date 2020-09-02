# -*- coding: utf-8 -*-

"""
Implementation of Karger's algorithm for minimal edge cut.

Randomized algorithm that returns a minimal edge cut with high probability.
This implementation is based on:
    
http://web.stanford.edu/class/archive/cs/cs161/cs161.1166/lectures/
lecture15.pdf
"""

import random, math
import networkx as nx

from asymmetree.datastructures.LinkedList import LinkedList
from asymmetree.datastructures.AVLTree import TreeSet


__author__ = 'David Schaller'


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
            
        return list(x.nodes), list(y.nodes), cut_value
    

if __name__ == '__main__':
    
    G = nx.Graph()
    G.add_edges_from([('a', 'b'), ('a', 'c'), ('a', 'd'),
                      ('b', 'd'), ('c', 'd'), ('c', 'e'), ('d', 'e')])
    karger = Karger(G)
    
    for V1, V2, cut_value in karger.generate(3):
        print(V1)
        print(V2)
        print(cut_value)