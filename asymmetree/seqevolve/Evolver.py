# -*- coding: utf-8 -*-


import numpy as np

from asymmetree.seqevolve.EvolvingSequence import EvoSeq, State
from asymmetree.seqevolve.Matrices import diagonalize
from asymmetree.seqevolve.Alignment import build_alignment


class Evolver:
    
    def __init__(self, subst_model, indel_model=None):
        
        self.subst_model = subst_model
        self.Q, self.freqs = subst_model.Q, subst_model.freqs
        self.dim = self.Q.shape[0]
        
        # compute the eigensystem
        self.eigvals, self.U, self.U_inv = diagonalize(self.Q, self.freqs)
        
        self.indel_model = indel_model
        
    
    def evolve_along_tree(self, T, start_length=200, start_seq=None):
        
        self.T = T
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
        
        while r > current_sum and index < len(p) - 1:
            index += 1
            current_sum += p[index]
        
        return index
    
    
    def _initialize_root(self, sequence):
        
        seq = EvoSeq()
        
        for x in self.subst_model.to_indices(sequence):
            seq.append(x, State.ROOT, self.site_counter)
            self.site_counter += 1
        
        return seq
    
    
    def _evolve(self, parent_seq, distance):
        
        child_seq = parent_seq.clone()
        
        if self.indel_model:
            self._generate_indels(child_seq, distance)
        
        self._substitute(child_seq, distance)        
        
        return child_seq
    
    
    def _substitute(self, sequence, t):
        
        # transtion probabilty matrix
        P = self.U  @  np.diag( np.exp(self.eigvals * t) )  @  self.U_inv
        
        r = np.random.random(len(sequence))
        pos = 0
        
        for site in sequence:
            
            if site.status == State.INHERITED:
                # mutate according to matrix P
                site._value = self.choose_index(P[site._value, :], r[pos])
            
            else:
                # choose random character
                site._value = self.choose_index(self.freqs, r[pos])
            
            pos += 1
    
    
    def _generate_indels(self, sequence, t):
        """Generates indels by a Gillespie process."""
        
        current_time = 0.0
        
        while current_time < t:
            
            ins_rate, del_rate = self.indel_model.get_rates(len(sequence))
            total_rate = ins_rate + del_rate
            
            waiting_time = np.random.exponential(1/total_rate) if total_rate > 0.0 else float('inf')
            current_time += waiting_time
            
            if current_time + waiting_time < t:
                
                r = np.random.random()
                if r < ins_rate / total_rate:
                    self._insertion(sequence)
                else:
                    self._deletion(sequence)
    
    
    def _insertion(self, sequence):
        
        d = self.indel_model.draw_length()
        pos = np.random.randint(-1, high=len(sequence))
        
        # inintialize insertion before the first item
        if pos == -1:
            current_element = sequence.append_left(None, State.INSERTION,
                                                   self.site_counter)
            self.site_counter += 1
            d -= 1
        # go to the element at position pos
        else:
            current_element = sequence.element_at(pos)
            
        for _ in range(d):
            current_element = sequence.insert_right_of(current_element,
                                                       None, State.INSERTION,
                                                       self.site_counter)
            self.site_counter += 1
    
    
    def _deletion(self, sequence):
        
        d = self.indel_model.draw_length()
        pos = np.random.randint(-d + 1, high=len(sequence))
        
        # deletion begins before or at the start of the sequence
        if pos <= 0:
            d += pos
            pos = 0
            
        sequence.remove_range(pos, d)
    
    
    def true_alignment(self):
        
        return build_alignment(self.T, self.sequences,
                               self.subst_model.get_alphabet())

 
if __name__ == "__main__":
    
    from SubstModel import SubstModel
    from IndelModel import IndelModel
    import asymmetree.simulator.TreeSimulator as ts
    
    subst_model = SubstModel("nucleotide", "JC69")
    indel_model = IndelModel(0.01, 0.01, length_model="zipf")
    
    evolver = Evolver(subst_model, indel_model=indel_model)
    print(evolver.Q)
    
    T = ts.build_species_tree(5)
    evolver.evolve_along_tree(T, start_length=30)
    
    for node, sequence in evolver.sequences.items():
        print(node.label, subst_model.to_sequence(sequence))
        
    alg_seq = evolver.true_alignment()
    for node, sequence in alg_seq.items():
        print(node.label, sequence)
    
    
    
