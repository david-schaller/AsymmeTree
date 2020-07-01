# -*- coding: utf-8 -*-

from asymmetree.tools.Sampling import Sampler


__author__ = 'David Schaller'


class IndelModel:
    
    
    def __init__(self, insertion_rate, deletion_rate,
                 length_distr=('zipf', 1.821),        # (Chang and Benner 2004)
                 min_length=1,
                 max_length=None,
                 **kwargs):
        
        if insertion_rate < 0.0 or deletion_rate < 0.0:
            raise ValueError("insertion and deletion rates must be non-negative")
        else:
            self._ins_rate = insertion_rate     # insertion rate per site
            self._del_rate = deletion_rate      # deletion rate per site
        
        if isinstance(length_distr, Sampler):
            self.sampler = length_distr
        else:
            self.sampler = Sampler(length_distr,
                                   maximum=max_length,
                                   discrete=True)
        
    
    def get_rates(self, seq_length):
        """Return the current insertion and deletion rate.
        
        Computes the current insertion and deletion rate according to the model
        and the length of the sequence."""
        
        # expected value may be infinite for zipf distribution, in which case
        # the minimal length is used, see 'Sampler'
        return ( (seq_length + 1) * self._ins_rate,
                 (seq_length + self.sampler._exp_val - 1) * self._del_rate )
        
    
    def draw_length(self):
        """Draw the length for an indel."""
        
        return self.sampler.draw()