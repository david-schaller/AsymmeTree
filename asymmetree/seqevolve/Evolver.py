# -*- coding: utf-8 -*-


import numpy as np

from EvolvingSequence import EvoSeq, State
from Matrices import diagonalize


class Evolver:
    
    def __init__(self, model):
        
        self.model = model
        self.Q, self.freqs = model.Q, model.freqs
        self.dim = self.Q.shape[0]
        
        # compute the eigensystem
        self.eigvals, self.U, self.U_inv = diagonalize(self.Q, self.freqs)
        
    
    def evolve_along_tree(self, T, start_length=200, start_seq=None):
        
        self.site_counter = 0
        self.sequences = {}
        
        if start_seq is None:
            root_seq = self._random_sequence(start_length)
        
        else:
            root_seq = self._initialize_root(start_seq)
            
        for v in T.preorder():
            
            # root of the tree, assign start sequence
            if v.parent is None:
                self.sequences[v] = root_seq
                
            # evolve the sequence along an edge
            else:
                parent_seq = self.sequences[v.parent]
                distance = v.dist
                self.sequences[v] = self._evolve(parent_seq, distance)
          
            
    def _random_positions(self, n):
        
        return np.random.choice(self.dim, n, p=self.freqs)
    
    
    def _random_sequence(self, n):
        
        seq = EvoSeq()
        
        for x in self._random_positions(n):
            seq.append(x, State.ROOT, self.site_counter)
            self.site_counter += 1
        
        return seq
    
    
    def choose_index(self, p, r):
        
        index = 0
        current_sum = p[index]
        
        while r > current_sum and index < len(p)-1:
            index += 1
            current_sum += p[index]
        
        return index
    
    
    def _initialize_root(self, sequence):
        
        seq = EvoSeq()
        
        for x in self.model.to_indices(sequence):
            seq.append(x, State.ROOT, self.site_counter)
            self.site_counter += 1
        
        return seq
    
    
    def _evolve(self, parent_seq, distance):
        
        child_seq = parent_seq.clone()
        
        # indels not yet implemented
        
        # substitutions
        self._substitute(child_seq, distance)        
        
        return child_seq
    
    
    def _substitute(self, sequence, t):
        
        # transtion probabilty matrix
        P = self.U  @  np.diag( np.exp(self.eigvals * t) )  @  self.U_inv
        
        r = np.random.random(len(sequence))
        pos = 0
        
        for site in sequence:
            
            # mutate according to matrix P
            if site.status == State.INHERITED:
                site._value = self.choose_index(P[site._value, :], r[pos])
            
            # choose random character
            else:
                site._value = self.choose_index(self.freqs, r[pos])
            
            pos += 1

 
if __name__ == "__main__":
    
    from Model import Model
    import asymmetree.simulator.TreeSimulator as ts
    
    model = Model("nucleotide", "JC69")
    evolver = Evolver(model)
    print(evolver.Q)
    
    T = ts.build_species_tree(5)
    evolver.evolve_along_tree(T, start_length=30)
    
    for node, sequence in evolver.sequences.items():
        print(node.label, model.to_sequence(sequence))
    
    
