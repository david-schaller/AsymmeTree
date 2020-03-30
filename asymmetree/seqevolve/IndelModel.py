# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import zeta


class IndelModel:
    
    def __init__(self, insertion_rate, deletion_rate,
                 length_model='zipf', params=None):
        
        if insertion_rate < 0.0 or deletion_rate < 0.0:
            raise ValueError("Insertion and deletion rates must be non-negative!")
        else:
            self._ins_rate = insertion_rate     # insertion rate per site
            self._del_rate = deletion_rate      # deletion rate per site
        
        length_model = length_model.lower()
        if length_model in ('zipf', 'negative_binomial'):
            self._length_model = length_model
        else:
            raise ValueError("Indel model '{}' is not available!".format(length_model))
        
        self._params = params
        
        # check the parameter according to length model
        self._check_params()
        
        # compute(/set) mean for the indel length distribution
        self._del_mean = self._distribution_mean()
        
    
    def get_rates(self, seq_length):
        
        return ( (seq_length + 1) * self._ins_rate,
                 (seq_length + self._del_mean - 1) * self._del_rate )
        
    
    def draw_indel_length(self):
        
        if self._length_model == 'zipf':
            return np.random.zipf( self._params['a'] )
                
        elif self._length_model == 'negative_binomial':
            return 1 + np.random.negative_binomial(self._params['r'], 1 - self._params['q'])
    
    
    def _check_params(self):
        
        if self._length_model == 'zipf':
            
            if not self._params or 'a' not in self._params:
                self._params = {'a': 1.821}     # (Chang and Benner 2004)
                
            elif not isinstance(self._params['a'], float) or self._params['a'] <= 1.0:
                raise ValueError("Invalid value for parameter 'a': {}".format(self._params['a']))
                
        if self._length_model == 'negative_binomial':
            
            if not self._params:
                self._params = {'r': 1, 'q': 0.5}
            
            if 'r' not in self._params:
                self._params['r'] = 1
            if 'q' not in self._params:
                self._params['q'] = 0.5
                
            if not isinstance(self._params['r'], int) or self._params['r'] < 1:
                raise ValueError("Parameter 'r' must be an int and > 0!")
            
            if (not isinstance(self._params['q'], float) or 
                self._params['q'] <= 0 or 
                self._params['q'] >= 1.0):
                raise ValueError("Parameter 'q' must be a float > 0 and < 1!")
            
            
    def _distribution_mean(self):
        
        if self._length_model == 'zipf':
            
            if self._params['a'] > 2.0:     # mean value is infinite for a <= 2
                return zeta(self._params['a'] - 1.0) / zeta(self._params['a'])
            else:
                return 0.0                  # maybe replace by some value > 0 ?
                
        if self._length_model == 'negative_binomial':
            return 1 + self._params['r'] * self._params['q'] / (1 - self._params['q'])
            