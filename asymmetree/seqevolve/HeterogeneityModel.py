# -*- coding: utf-8 -*-

import numpy as np


class HeterogeneityModel:
    
    def __init__(self):
        
        pass
    
    
    def initialize(self, subst_model):
        
        self.subst_model = subst_model
        self._compute_row_sums()
        
        
    def get_rate_factors(self, n):
        
        return [1.0 for _ in range(n)]