# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx


class AlignmentBuilder:
    
    
    def __init__(self, tree, sequence_dict, include_inner=True):
        
        self.tree = tree
        self.sequence_dict = sequence_dict
        self.include_inner
        
    
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
            
            for first, second in seq.element_pairs:
                G.add_edge(first.site_id, second.site_id)
                
        self.sites = list(nx.topological_sort(G))
                
    
    def _alignment_matrix(self):
        
        alignment = np.zeros((len(self.nodes), len(self.sites)),
                             dtype=np.int8)
        positions = [self.sequence_dict[v]._first for v in self.nodes]
        
        for j in range(len(self.sites)):
            for i in range(len(positions)):
                
                if self.sites[j] == positions[i].site_id:
                    alignment[i, j] = positions[i].site_id._value
                    positions[i] = positions[i]._next
                else:
                    alignment[i, j] = -1
                    
        return alignment
        
        
        
            
            
        