# -*- coding: utf-8 -*-

import numpy as np

from asymmetree.seqevolve.EvolvingSequence import EvoSeq, State
from asymmetree.seqevolve.Alignment import AlignmentBuilder
from asymmetree.file_io.SeqFileIO import write_alignment


__author__ = 'David Schaller'


class Evolver:
    
    def __init__(self, subst_model,
                 indel_model=None, het_model=None,
                 jump_chain=False,
                 **kwargs):
        
        self.subst_model = subst_model
        
        self.indel_model = indel_model
        self.het_model = het_model
        
        if (jump_chain or 
            (het_model and het_model.sitewise) or
            (het_model and het_model.classes > 5)):     # determine best cutoff
            self._jump_chain = True
        else:
            self._jump_chain = False
        
    
    def evolve_along_tree(self, tree, start_length=200, start_seq=None):
        
        self.tree = tree
        self.site_counter = 0
        self.sequences = {}
        
        if start_seq is None:
            root_seq = self._random_sequence(start_length)
        else:
            root_seq = self._initialize_with_sequence(start_seq)
            
        if self.het_model:
            self.het_model.assign(root_seq)
        
        for v in tree.preorder():
            
            if v.parent is None:
                self.sequences[v] = root_seq
            else:
                self.sequences[v] = self._evolve(self.sequences[v.parent], v.dist)
                
        return self.sequences
          
            
    def _random_positions(self, n):
        
        return np.random.choice(len(self.subst_model.alphabet), n,
                                p=self.subst_model.freqs)
    
    
    def _random_sequence(self, n):
        
        seq = EvoSeq()
        
        for x in self._random_positions(n):
            seq.append(x, State.ROOT, self.site_counter)
            self.site_counter += 1
        
        return seq
    
    
    def _initialize_with_sequence(self, sequence):
        
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
        
        if not self._jump_chain:
            self._substitute(child_seq, distance)               # apply transition prob. matrix
        else:
            self._substitute_jump_chain(child_seq, distance)    # apply jump chain
        
        return child_seq
    
    
    def true_alignment(self, include_inner=True, write_to=None,
                       alignment_format='phylip'):
        
        alg_builder = AlignmentBuilder(self.tree, self.sequences,
                                       self.subst_model.alphabet,
                                       include_inner=include_inner)
        
        alignment = alg_builder.build()
        
        if write_to:
            write_alignment(write_to, alignment,
                            alignment_format=alignment_format)
        
        return alignment
    
    
    # --------------------------------------------------------------------------
    #               Substitution by transition probability matrix
    # --------------------------------------------------------------------------
    
    def _substitute(self, sequence, t):
        
        # (cumulative) transition probability matrices
        P = self._cumulative_P_matrices(sequence, t)
        
        r = np.random.random(len(sequence))
        pos = 0
        
        for site in sequence:
            
            if site.status == State.INHERITED:
                # mutate according to matrix of the corresponding class
                site._value = np.argmax(P[site.rate_class][site._value, :] > r[pos])
            
            else:
                # choose random character
                site._value = np.argmax(self.subst_model.freqs_cumulative > r[pos])
            
            pos += 1
    
    
    def _cumulative_P_matrices(self, sequence, t):
        
        P = {}
        
        if not self.het_model:
            P[0] = np.cumsum(self.subst_model.transition_prob_matrix(t),
                             axis=1)

        else:
            for site in sequence:
                
                c = site.rate_class
                
                if c not in P:
                    
                    r = site.rate_factor
                    
                    if r > 0.0:
                        P[c] = np.cumsum(self.subst_model.transition_prob_matrix(r * t),
                                         axis=1)
                    elif r == 0:
                        P[c] = np.cumsum(np.identity(len(self.subst_model.alphabet)),
                                         axis=1)
        
        return P
        
    
    # --------------------------------------------------------------------------
    #                        Substitution by Jump chain
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
        
        # choose random characters for insertions
        r = np.random.random(sequence.count_status(State.INSERTION))
        pos = 0
        for site in sequence:
            
            if site.status == State.INSERTION:
                site._value = np.argmax(self.subst_model.freqs_cumulative > r[pos])
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