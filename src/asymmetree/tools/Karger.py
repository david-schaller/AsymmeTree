# -*- coding: utf-8 -*-

"""
Implementation of Karger's algorithm for minimal edge cut.

Randomized algorithm that returns a minimal edge cut with high probability.
This implementation is based on:
    
http://web.stanford.edu/class/archive/cs/cs161/cs161.1166/lectures/
lecture15.pdf
"""

import random
import networkx as nx

from asymmetree.datastructures.LinkedList import LinkedList


__author__ = 'David Schaller'


class Supernode:
    
    def __init__(self, initial_node):
        
        self.nodes = LinkedList()
        self.nodes.append(initial_node)
        
    
    def __len__(self):
        
        return len(self.nodes)
    
    
    def merge_from(self, other):
        
        self.nodes.concatenate(other.nodes)
        
        
class Karger():
    
    def __init__(self, graph):
        
        if not isinstance(graph, nx.Graph):
            raise TypeError('must be a NetworkX graph')
            
        if not graph.order() >= 2:
            raise RuntimeError('graph must have >= 2 nodes')
        
        self.supernodes = set()
        self.supernode_dict = {}    # maps node --> supernode
        self.edges = set()          # set of all (non-loop) edges
        self.superedges = {}        # maps superedge --> edges
        
        for v in graph.nodes():
            sn = Supernode(v)
            self.supernode_dict[v] = sn
            self.supernodes.add(sn)
        
        for u, v in graph.edges():
            u_sn = self.supernode_dict[u]
            v_sn = self.supernode_dict[v]
            
            self.superedges[u_sn, v_sn] = LinkedList()
            self.superedges[u_sn, v_sn].append( (u, v) )
            self.edges.add( (u, v) )
            
            
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
            
            # merge the edge lists E_ad and E_bd
            edges_xd = LinkedList()
            for se in ((a, d), (d, a), (b, d), (d, b)):
                if se in self.superedges:
                    edges_xd.concatenate(self.superedges[se])
                    del self.superedges[se]
            if edges_xd:
                self.superedges[x, d] = edges_xd
                
                
    def karger(self):
        
        while len(self.supernodes) > 2:
            
            # find a better solution for choice !!!
            u, v = random.choice(tuple(self.edges))
            u_sn = self.supernode_dict[u]
            v_sn = self.supernode_dict[v]
            
            self._merge(u_sn, v_sn)
            
            for se in ((u_sn, v_sn), (v_sn, u_sn)):
                if se in self.superedges:
                    self.edges.difference_update(self.superedges[se])
                    del self.superedges[se]
        
        x, y = list(self.supernodes)
        V1 = list(x.nodes)
        V2 = list(y.nodes)
        if (x, y) in self.superedges:
            cut_value = len(self.superedges[x, y])
            print('cut_edges', list(self.superedges[x, y]))
        else:
            cut_value = len(self.superedges[y, x])
            print('cut_edges', list(self.superedges[y, x]))
            
        return V1, V2, cut_value
    

if __name__ == '__main__':
    
    G = nx.Graph()
    G.add_edges_from([('a', 'b'), ('a', 'c'), ('a', 'd'),
                      ('b', 'd'), ('c', 'd'), ('c', 'e'), ('d', 'e')])
    k = Karger(G)
    print(k.superedges)
    V1, V2, cut_value = k.karger()
    print(V1)
    print(V2)
    print(cut_value)