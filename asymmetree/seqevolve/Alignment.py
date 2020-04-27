# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx


__author__ = "David Schaller"


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
    

def write_to_file(filename, alignment, al_format='phylip'):
    
    with open(filename, 'w') as f:
        
        if al_format == 'phylip':
            _write_phylip(f, alignment)
    
        elif al_format == 'clustal':
            _write_clustal(f, alignment)
        
        elif al_format == 'pretty':
            _write_pretty(f, alignment)
            
        else:
            raise ValueError("Alignment format '{}' is not available!".format(al_format))
            

def _check_alignment(alignment):
    
    max_length = 0                          # maximal label length
    seq_length = None                       # length of the aligned sequences
    for node, seq in alignment.items():
        
        if len(node.label) > max_length:
            max_length = len(node.label)
            
        if seq_length is None:
            seq_length = len(seq)
        elif seq_length != len(seq):
            raise ValueError("Aligned sequences must have the same length!")
            
    return max_length, seq_length
            
            
def _write_phylip(f, alignment):
    
    max_length, seq_length = _check_alignment(alignment)
    max_length = max(max_length, 9)
    
    
    f.write("  {} {} i".format(len(alignment), seq_length))
    
    format_str = "\n{:" + str(max_length+1) + "}"
    current = 0
    
    while current < seq_length:
        
        end = min(seq_length, current+50)
        
        for node, seq in alignment.items():
            
            if current == 0:
                f.write(format_str.format(node.label))
            else:
                f.write(format_str.format(''))
                
            f.write(" ".join( seq[i:min(i+10,end)] for i in range(current, end, 10) ))
                
        if end != seq_length:
            f.write("\n")
        
        current += 50


def _write_clustal(f, alignment):
    
    max_length, seq_length = _check_alignment(alignment)
    
    f.write('CLUSTAL W (1.8) multiple sequence alignment\n\n')
    
    format_str = "\n{:" + str(max_length+4) + "}{}"
    current = 0
    
    while current < seq_length:
        
        end = min(seq_length, current+60)
        
        for node, seq in alignment.items():
            f.write(format_str.format(node.label, seq[current:end]))
                
        if end != seq_length:
            f.write("\n\n")
        
        current += 60
    
        
def _write_pretty(f, alignment):
    
    _, seq_length = _check_alignment(alignment)
    
    f.write('  1          11         21         31         41       50\n')
    f.write('  |          |          |          |          |        |')
    
    current = 0
    
    while current < seq_length:
        
        end = min(seq_length, current+50)
        
        for node, seq in alignment.items():
            
            count = end - current - seq[current:end].count('-')
            seq_string = " ".join( seq[i:min(i+10,end)] for i in range(current, end, 10) )
            f.write("\n  {:54}{:>6} {}".format(seq_string, count, node.label))
                
        if end != seq_length:
            f.write("\n\n\n")
        
        current += 50
