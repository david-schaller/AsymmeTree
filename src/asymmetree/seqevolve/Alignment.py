# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx


__author__ = 'David Schaller'


class AlignmentBuilder:
    
    
    def __init__(self, tree, sequence_dict, alphabet, include_inner=True):
        
        self.tree = tree
        self.sequence_dict = sequence_dict
        self.include_inner = include_inner
        
        if alphabet[-1] == '-':
            self.alphabet = alphabet
        else:
            self.alphabet = alphabet + '-'
        
        
    def build(self):
        
        self._get_preorder()
        self._sort_sites()
        self._alignment_matrix()
        
        return self._sequences()
        
    
    def _get_preorder(self):
        
        if self.include_inner:
            self.nodes = [v for v in self.tree.preorder()]
        else:
            self.nodes = [v for v in self.tree.preorder() if not v.children]
            
        self.node_index = {item: index for index, item in enumerate(self.nodes)}
        
    
    def _sort_sites(self):
        
        G = nx.DiGraph()
        
        for node in self.nodes:
            
            seq = self.sequence_dict[node]
            
            # handle case |seq| = 1 correctly
            if len(seq) == 1:
                G.add_node(seq[0].site_id)
            
            # add nodes and edges to auxiliary graph
            for first, second in seq.element_pairs():
                G.add_edge(first.site_id, second.site_id)
                
        self.sites = list(nx.topological_sort(G))
                
    
    def _alignment_matrix(self):
        
        self.alignment = np.zeros((len(self.nodes), len(self.sites)),
                             dtype=np.int8)
        positions = [self.sequence_dict[v]._first for v in self.nodes]
        
        for j in range(len(self.sites)):
            for i in range(len(positions)):
                
                if positions[i] and self.sites[j] == positions[i].site_id:
                    self.alignment[i, j] = positions[i]._value
                    positions[i] = positions[i]._next
                else:
                    self.alignment[i, j] = -1
    
    
    def _sequences(self):
        
        result = {}
        
        for i in range(len(self.nodes)):
            sequence = "".join(self.alphabet[x] for x in self.alignment[i, :])
            result[self.nodes[i]] = sequence
        
        return result