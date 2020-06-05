# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import zeta


__author__ = 'David Schaller'


class IndelModel:
    
    
    def __init__(self, insertion_rate, deletion_rate,
                 length_model='zipf',
                 max_length=False,
                 **kwargs):
        
        if insertion_rate < 0.0 or deletion_rate < 0.0:
            raise ValueError("insertion and deletion rates must be non-negative")
        else:
            self._ins_rate = insertion_rate     # insertion rate per site
            self._del_rate = deletion_rate      # deletion rate per site
        
        length_model = length_model.lower()
        if length_model in ('zipf', 'negative_binomial'):
            self._length_model = length_model
        else:
            raise ValueError("indel model '{}' is not available".format(length_model))
        
        self._max_length = max_length
        self._params = kwargs
        
        # check the parameter according to length model
        self._check_params()
        
        # compute(/set) mean for the indel length distribution
        self._del_mean = self._distribution_mean()
        
    
    def get_rates(self, seq_length):
        """Return the current insertion and deletion rate.
        
        Computes the current insertion and deletion rate according to the model
        and the length of the sequence."""
        
        return ( (seq_length + 1) * self._ins_rate,
                 (seq_length + self._del_mean - 1) * self._del_rate )
        
    
    def draw_length(self):
        """Draw the length for an indel."""
        
        while True:
            
            if self._length_model == 'zipf':
                d = np.random.zipf( self._params['a'] )
                    
            elif self._length_model == 'negative_binomial':
                d = 1 + np.random.negative_binomial(self._params['r'], 1 - self._params['q'])
            
            if self._max_length is False or d <= self._max_length:
                break
            
        return d
            
    
    def _check_params(self):
        
        if self._max_length is not False:
            
            if not isinstance(self._max_length, int) or self._max_length < 1:
                raise ValueError("maximal indel length must be an int >0")
        
        if self._length_model == 'zipf':
            
            if 'a' not in self._params:
                self._params = {'a': 1.821}     # (Chang and Benner 2004)
                
            elif not isinstance(self._params['a'], float) or self._params['a'] <= 1.0:
                raise ValueError("invalid value for parameter 'a': {}".format(self._params['a']))
                
        if self._length_model == 'negative_binomial':
            
            if 'r' not in self._params:
                self._params['r'] = 1
            if 'q' not in self._params:
                self._params['q'] = 0.5
                
            if not isinstance(self._params['r'], int) or self._params['r'] < 1:
                raise ValueError("parameter 'r' must be an int and >0")
            
            if (not isinstance(self._params['q'], float) or 
                self._params['q'] <= 0 or 
                self._params['q'] >= 1.0):
                raise ValueError("parameter 'q' must be a float >0 and <1")
            
            
    def _distribution_mean(self):
        
        if self._length_model == 'zipf':
            
            if self._params['a'] > 2.0:     # mean value is infinite for a <= 2
                return zeta(self._params['a'] - 1.0) / zeta(self._params['a'])
            else:
                return 0.0                  # maybe replace by some value > 0 ?
                
        if self._length_model == 'negative_binomial':
            return 1 + self._params['r'] * self._params['q'] / (1 - self._params['q'])
            