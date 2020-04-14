# -*- coding: utf-8 -*-


import numpy as np

from asymmetree.seqevolve.EvolvingSequence import EvoSeq, State
from asymmetree.seqevolve.Alignment import AlignmentBuilder, write_to_file


class Evolver:
    
    def __init__(self, subst_model, indel_model=None, het_model=None):
        
        self.subst_model = subst_model
        self.eigvals, self.U, self.U_inv = subst_model.eigensystem()
        
        self.indel_model = indel_model
        self.het_model = het_model
        
        self.jump_chain = True
        
    
    def evolve_along_tree(self, T, start_length=200, start_seq=None):
        
        self.T = T
        self.site_counter = 0
        self.sequences = {}
        
        if start_seq is None:
            root_seq = self._random_sequence(start_length)
        else:
            root_seq = self._initialize_root(start_seq)
            
        if self.het_model:
            self.het_model.assign(root_seq)
        
        for v in T.preorder():
            
            if v.parent is None:
                self.sequences[v] = root_seq
            else:
                self.sequences[v] = self._evolve(self.sequences[v.parent], v.dist)
          
            
    def _random_positions(self, n):
        
        return np.random.choice(len(self.subst_model.alphabet), n,
                                p=self.subst_model.freqs)
    
    
    def _random_sequence(self, n):
        
        seq = EvoSeq()
        
        for x in self._random_positions(n):
            seq.append(x, State.ROOT, self.site_counter)
            self.site_counter += 1
        
        return seq
    
    
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
            
            if self.het_model:
                self.het_model.assign(child_seq)
        
        if not self.jump_chain:
            self._substitute(child_seq, distance)               # apply matrix substitution
        else:
            self._substitute_jump_chain(child_seq, distance)    # apply jump chain substitution
        
        return child_seq
    
    
    def true_alignment(self, include_inner=True):
        
        alg_builder = AlignmentBuilder(self.T, self.sequences,
                                       self.subst_model.alphabet,
                                       include_inner=include_inner)
        
        return alg_builder.build()
    
    
    # --------------------------------------------------------------------------
    #                           Matrix substitution
    # --------------------------------------------------------------------------
    
    def _substitute(self, sequence, t):
        
        # transtion probability matrices
        P = self._P_matrices(sequence, t)
        
        r = np.random.random(len(sequence))
        pos = 0
        
        for site in sequence:
            
            if site.status == State.INHERITED:
                # mutate according to matrix P_c of the corresponding class
                P_c = P[site.rate_class]
                site._value = self._choose_index(P_c[site._value, :], r[pos])
            
            else:
                # choose random character
                site._value = self._choose_index(self.subst_model.freqs, r[pos])
            
            pos += 1
            
    
    def _choose_index(self, p, r):
        
        index = 0
        current_sum = p[index]
        
        while r > current_sum and index < len(p) - 1:
            index += 1
            current_sum += p[index]
        
        return index
    
    
    def _P_matrices(self, sequence, t):
        
        P = {}
        
        if not self.het_model:
            P[0] = self.U  @  np.diag( np.exp(self.eigvals * t) )  @  self.U_inv
            
        else:
            for site in sequence:
                
                c = site.rate_class
                
                if c not in P:
                    
                    r = site.rate_factor
                    
                    if r > 0.0:
                        P[c] = self.U  @  np.diag( np.exp(self.eigvals * r * t) )  @  self.U_inv
                    elif r == 0:
                        P[c] = np.identity(len(self.subst_model.alphabet))
        
        return P
        
    
    # --------------------------------------------------------------------------
    #                        Jump chain substitution
    # --------------------------------------------------------------------------
    
    def _substitute_jump_chain(self, sequence, t):
        
        total_rate = self._total_subst_rate(sequence, exclude_inserted=True)
        
        # jump chain algorithm for inherited positions
        current_time = 0.0
        while current_time < t:
            
            current_time += np.random.exponential(1/total_rate) if total_rate > 0.0 else float('inf')
            
            if current_time < t:
                
                site, mutation = self._draw_substitution(sequence, total_rate)
                
                total_rate += site.rate_factor * self.subst_model.Q[site._value, site._value]
                total_rate -= site.rate_factor * self.subst_model.Q[mutation, mutation]
                site._value = mutation
        
        # choose random character for insertions
        r = np.random.random(sequence.count_status(State.INSERTION))
        pos = 0
        for site in sequence:
            
            if site.status == State.INSERTION:
                site._value = self._choose_index(self.subst_model.freqs, r[pos])
                pos += 1
                
    
    def _total_subst_rate(self, sequence, exclude_inserted=True):
        
        total = 0.0
        
        for site in sequence:
            
            if exclude_inserted and site.status == State.INSERTION:
                continue
            
            total -= site.rate_factor * self.subst_model.Q[site._value, site._value]
            
        return total
    
    
    def _draw_substitution(self, sequence, total_rate):
        
        chosen_site, chosen_mutation = None, None
        r = np.random.uniform(low=0.0, high=total_rate)
        
        current_sum = 0.0
        for site in sequence:
            
            if site.status == State.INSERTION:
                continue
            
            # use negative value on the diagonal of the rate matrix Q
            site_rate = - site.rate_factor * self.subst_model.Q[site._value, site._value]
            current_sum += site_rate
            
            if current_sum > r:
                chosen_site = site
                
                # rescale the "rest" of r for application on the rate matrix Q
                r = (r - current_sum + site_rate) / site.rate_factor
                
                break
            
        current_sum = 0.0
        for i in range(len(self.subst_model.alphabet)):
            
            if i != site._value:            # skip the entry on diagonal
                current_sum += self.subst_model.Q[site._value, i]
                
                if current_sum > r:
                    chosen_mutation = i
                    break
                    
        return chosen_site, chosen_mutation
    
    
    # --------------------------------------------------------------------------
    #                               Indels
    # --------------------------------------------------------------------------
        
    def _generate_indels(self, sequence, t):
        """Generates indels by a Gillespie process."""
        
        current_time = 0.0
        
        while current_time < t:
            
            ins_rate, del_rate = self.indel_model.get_rates(len(sequence))
            total_rate = ins_rate + del_rate
            
            current_time += np.random.exponential(1/total_rate) if total_rate > 0.0 else float('inf')
            
            if current_time < t:
                
                r = np.random.random()
                if r < ins_rate / total_rate:
                    self._insertion(sequence)
                else:
                    self._deletion(sequence)
    
    
    def _insertion(self, sequence):
        
        d = self.indel_model.draw_length()
        pos = np.random.randint(-1, high=len(sequence))
        
        # initialize insertion before the first item
        if pos == -1:
            current_site = sequence.append_left(None, State.INSERTION,
                                                self.site_counter)
            self.site_counter += 1
            d -= 1
            
        # go to the element at position pos
        else:
            current_site = sequence.element_at(pos)
            
        for _ in range(d):
            current_site = sequence.insert_right_of(current_site,
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

 
if __name__ == "__main__":
    
    from SubstModel import SubstModel
    from IndelModel import IndelModel
    import asymmetree.simulator.TreeSimulator as ts
    
    subst_model = SubstModel("aa", "WAG")
    indel_model = IndelModel(0.01, 0.01, length_model="zipf")
    
    evolver = Evolver(subst_model, indel_model=indel_model)
    print(evolver.subst_model.Q)
    
    T = ts.build_species_tree(5)
    evolver.evolve_along_tree(T, start_length=150)
    
    for node, sequence in evolver.sequences.items():
        print(node.label, subst_model.to_sequence(sequence))
        
    alg_seq = evolver.true_alignment()
    for node, sequence in alg_seq.items():
        print(node.label, sequence)
        
    write_to_file("../../validation/testfile.alignment", alg_seq, al_format='phylip')
    
    
    
