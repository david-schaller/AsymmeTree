# -*- coding: utf-8 -*-

import collections, itertools, random

from asymmetree.tools import Tree


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


class SimpleGraph:
    
    def __init__(self, initial=None):
        if isinstance(initial, dict):
            self.adj_list = initial
        elif isinstance(initial, collections.Iterable):
            self.adj_list = {v: set() for v in initial}
        else:
            self.adj_list = {}
    
    
    def add_node(self, v):
        if v not in self.adj_list:
            self.adj_list[v] = set()
    
    
    def has_edge(self, u, v):
        if u in self.adj_list and v in self.adj_list[u]:
            return True
        else:
            return False
        
    
    def neighbors(self, v):
        return iter(self.adj_list[v])
    
    
    def add_edge(self, u, v):
        if u == v:
            raise ValueError("No loops allowed!")
        # add nodes if not yet in the graph
        self.add_node(u)
        self.add_node(v)
        
        self.adj_list[u].add(v)
        self.adj_list[v].add(u)
    
    
    def connected_comp(self, make_undirected=False):
        """Determines connected components by breadth-first search.
        
        Return a list of dictionaries.
        """
        result = []                                             # a list of dictionaries
        visited = set()                                         # already visited nodes
        Q = collections.deque()                                 # initialize queue
        for start_node in self.adj_list.keys():
            if start_node in visited:
                continue
            conn_comp = dict()                                  # dictionary that represents a c.c.
            conn_comp[start_node] = self.adj_list[start_node]
            visited.add(start_node)
            Q.append(start_node)
            while Q:
                v = Q[0]
                for u in self.adj_list[v]:
                    if u not in visited:
                        visited.add(u)
                        Q.append(u)
                        conn_comp[u] = self.adj_list[u]
                Q.popleft()
            result.append(conn_comp)
        return result
    
    
    def graphs_equal(self, other):
        """Returns whether two undirected graphs are equal."""
        
        set1 = set(self.adj_list.keys())
        set2 = set(other.adj_list.keys())
        
        if set1 != set2:
            return False
        
        for v in set1:
            if self.adj_list[v] != other.adj_list[v]:
                return False
        
        return True
    
    
    def symmetric_diff(self, other):
        """Returns the number of edges in the symmetric difference."""
        
        set1 = set(self.adj_list.keys())
        set2 = set(other.adj_list.keys())
        
        if set1 != set2:
            raise ValueError("Graphs do not have the same vertex set")
            return
        
        V = list(set1)
        sym_diff_number = 0
        
        for i in range(0, len(V)-1):
            for j in range(i+1, len(V)):
                if self.has_edge(V[i], V[j]) and not other.has_edge(V[i], V[j]):
                    sym_diff_number += 1
                elif not self.has_edge(V[i], V[j]) and other.has_edge(V[i], V[j]):
                    sym_diff_number += 1
        
        return sym_diff_number
    
    
    def get_complement(self):
        
        compl_graph = SimpleGraph(initial=self.adj_list.keys())
        
        V = list(self.adj_list.keys())
        
        for i in range(len(V)-1):
            for j in range(i+1, len(V)):
                if not self.has_edge(V[i], V[j]):
                    compl_graph.add_edge(V[i], V[j])
        
        return compl_graph
    
    
    @staticmethod
    def random_graph(N, p=0.5):
        """Construct a random graph on N nodes.
        
        Keyword arguments:
            p - probability that an edge xy is inserted"""
        V = [i for i in range(1, N+1)]
        
        G = SimpleGraph(initial=V)
        
        for i in range(0, len(V)-1):
            for j in range(i+1, len(V)):
                if random.random() < p:
                    G.add_edge(V[i], V[j])
        
        return G
    

class CotreeNode(Tree.TreeNode):
    
    __slots__ = ['aux_counter']
    
    
    def __init__(self, ID, label=None, parent=None,
                 aux_counter=0):
        
        super().__init__(ID, label=label)
        
        self.aux_counter = aux_counter      # auxiliary counting variable
        
    
    def __str__(self):
        if not self.children:
            return str(self.ID)
        elif self.label == "series":
            return "<1>"
        elif self.label == "parallel":
            return "<0>"
        else:
            return "<>"
        

class Cotree(Tree.Tree):
    
    def __init__(self, root):
        super().__init__(root)
        
    
    def to_cograph(self):
        self.supply_leaves()
        G = SimpleGraph()
        
        for v in self.root.leaves:
            G.add_node(v.ID)
        
        for u in self.preorder():
            if u.label == "series":
                for child1, child2 in itertools.combinations(u.children, 2):
                    for v1 in child1.leaves:
                        for v2 in child2.leaves:
                            G.add_edge(v1.ID, v2.ID)
        
        return G
        
        
    @staticmethod
    def cotree(G):
        """Checks if a graph is a cograph and returns its cotree.
        
        Simple O(n^3) implementation."""
        
        def build_cotree(G, label=None):
            v = CotreeNode(None)
            v.label = label
            child_nodes = []
            
            if len(G.adj_list) == 1:
                v.label = "leaf"
                for ID in G.adj_list.keys():
                    v.ID = ID
            
            elif v.label is None:
                ccs = G.connected_comp()
                if len(ccs) > 1:
                    v.label = "parallel"
                    
                    for cc in ccs:
                        G_i = SimpleGraph(initial=cc)
                        v_i = build_cotree(G_i, label="series")
                        if v_i:
                            child_nodes.append(v_i)
                        else:
                            return False
                else:
                    G = G.get_complement()
                    ccs = G.connected_comp()
                    if len(ccs) > 1:
                        v.label = "series"
                        
                        for cc in ccs:
                            G_i = SimpleGraph(initial=cc)
                            v_i = build_cotree(G_i, label="parallel")
                            if v_i:
                                child_nodes.append(v_i)
                            else:
                                return False
                    else:
                        return False
            
            else:
                G = G.get_complement()
                ccs = G.connected_comp()
                if len(ccs) > 1:                 
                    for cc in ccs:
                        child_label = "series" if v.label=="parallel" else "parallel"
                        G_i = SimpleGraph(initial=cc)
                        v_i = build_cotree(G_i, label=child_label)
                        if v_i:
                            child_nodes.append(v_i)
                        else:
                            return False
                else:
                    return False
            
            for child in child_nodes:
                v.children.append(child)
                child.parent = v
            
            return v
        
        root = build_cotree(G)
        
        if root:
            return Cotree(root)
        else:
            print("Not a Cograph!")
            return False
       
    
    @staticmethod
    def random_cotree(N_leaves, force_series_root=False):
        """Creates a random cotree."""
        
        root = CotreeNode(ID=0)
        cotree = Cotree(root)
        node_list = [root]
        nr, leaf_count = 1, 1
        
        while leaf_count < N_leaves:
            node = random.choice(node_list)
            if not node.children:                       # avoid nodes with outdegree 1
                new_child1 = CotreeNode(nr)
                new_child2 = CotreeNode(nr+1)
                node.add_child(new_child1)
                node.add_child(new_child2)
                node_list.extend(node.children)
                nr += 2
                leaf_count += 1
            elif node.children:                         # add only one child if there are already children
                new_child = CotreeNode(nr)
                node.add_child(new_child)
                node_list.append(new_child)
                nr += 1
                leaf_count += 1
                
        for v in cotree.preorder():                     # assign labels ('series', 'parallel', 'leaf')
            if not v.children:
                v.label == "leaf"
            elif v.parent is None:
                if force_series_root:
                    v.label = "series"
                else:
                    v.label = "series" if random.random() < 0.5 else "parallel"
            else:
                v.label = "series" if v.parent.label == "parallel" else "parallel"
                
        return cotree
    
    
if __name__ == "__main__":
    
    from pprint import pprint
    
    cotree = Cotree.random_cotree(100)
    print(cotree.to_newick())
    cograph = cotree.to_cograph()
    #pprint(cograph.adj_list)
    
    cotree2 = Cotree.cotree(cograph)
    print(cotree2.to_newick())
    cograph2 = cotree2.to_cograph()
    #pprint(cograph2.adj_list)
    
    print(cograph.graphs_equal(cograph2))
    
    cograph3 = SimpleGraph()
    cograph3.add_edge("a", "b")
    cograph3.add_edge("b", "d")
    cograph3.add_edge("c", "d")
    cograph3.add_edge("b", "c")
    cotree3 = Cotree.cotree(cograph3)
    print(cotree3.to_newick())
    
    cograph4 = SimpleGraph.random_graph(10, p=0.5)
    pprint(cograph4.adj_list)
    cotree4 = Cotree.cotree(cograph4)
    if cotree4:
        print(cotree4.to_newick())