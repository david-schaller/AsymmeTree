# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import zeta


__author__ = 'David Schaller'


class Sampler:
    
    
    def __init__(self, params, minimum=None, maximum=None, discrete=False,
                 shift=0):
        """Constructor for class Sampler.
        
        Keyword argument:
            minimum - minimum value to be sampled
            maximum - maximum value to be sampled
            discrete - if True, round sampled values in case of a continuous
                distribution
            shift - shift distribution along x-axis
        """
        
        self._params = params
        self._discrete = discrete
        self._min = minimum
        self._max = maximum
        self._shift = shift
        
        self._check_distribution()
        
    
    def _check_distribution(self):
        
        if self._min is not None and not isinstance(self._min, (int, float)):
            raise ValueError('minimum must be a number')
        if self._max is not None  and not isinstance(self._max, (int, float)):
            raise ValueError('maximum must be a number')
        if not isinstance(self._shift, (int, float)):
            raise ValueError('shift value must be a number')
            
        
        if isinstance(self._params, (int, float)):
            
            self._distr = 'constant'
            self.draw = self._draw_constant
            self._exp_val = (round(self._params) if self._discrete else
                             self._params)
            
        elif isinstance(self._params, (tuple, list)) and len(self._params) >= 2:
            
            if self._params[0] == 'constant':
                
                val = self._params[1]
                if not isinstance(val, (int, float)):
                    raise ValueError('constant value must be a number')
                    
                self._distr = 'constant'
                self.draw = self._draw_constant
                self._exp_val = round(val) if self._discrete else val
                
            elif self._params[0] == 'uniform':
                
                a = self._params[1]
                if not len(self._params) >= 3:
                    raise ValueError('uniform distr. requires 2 parameters')
                b = self._params[2]
                
                if (not isinstance(a, (int, float)) or 
                    not isinstance(b, (int, float)) or b < a):
                    raise ValueError('parameters a and b for uniform distr.'\
                                     'must be numbers and a <= b')
                    
                self._distr = 'uniform'
                self.draw = self._draw_uniform
                self._exp_val = (a + b) / 2
                self._a = a
                self._b = b
                
            elif self._params[0] == 'discrete_uniform':
                
                a = self._params[1]
                if not len(self._params) >= 3:
                    raise ValueError('discrete uniform distr. requires 2 '\
                                     'parameters')
                b = self._params[2]
                
                if (not isinstance(a, int) or 
                    not isinstance(b, int) or b < a):
                    raise ValueError('parameters a and b for discrete uniform '\
                                     'distr. must be ints and a <= b')
                    
                self._distr = 'discrete_uniform'
                self.draw = self._draw_discrete_uniform
                self._exp_val = (a + b) / 2
                self._a = a
                self._b = b
                    
            elif self._params[0] == 'gamma':
                
                shape = self._params[1]
                if not len(self._params) >= 3:
                    raise ValueError('gamma distr. requires 2 parameters')
                scale = self._params[2]
                
                if (not isinstance(shape, float) or shape <= 0.0 or 
                    not isinstance(scale, float) or scale <= 0.0):
                    raise ValueError('scale and shape parameters for gamma '\
                                     'distribution must be floats >0.0')
                    
                self._distr = 'gamma'
                self.draw = self._draw_gamma
                self._exp_val = shape * scale
                self._shape = shape
                self._scale = scale
            
            elif self._params[0] == 'gamma_mean':
                
                mean = self._params[1]
                
                if not isinstance(mean , (int, float)) or mean <= 0.0:
                    raise ValueError('mean of gamma distribution must be '\
                                     'a number >0.0')
                
                self._distr = 'gamma'
                self.draw = self._draw_gamma
                self._exp_val = mean
                self._shape = 1.0
                self._scale = mean / self._shape
                
                
            elif self._params[0] == 'exponential':
                
                rate = self._params[1]
                
                if not isinstance(mean , float) or rate < 0.0:
                    raise ValueError('rate must be a float >=0.0')
                
                self._distr = 'exponential'
                self.draw = self._draw_exponential
                self._exp_val = float('inf') if rate == 0.0 else 1/rate
                self._rate = rate
                
                
            elif self._params[0] == 'zipf':
                
                a = self._params[1]
                    
                if not isinstance(a, float) or a <= 1.0:
                    raise ValueError('parameter a for zipf distr. must be a '\
                                     'float >1')
                
                # expected value is infinite for a <= 2
                if a > 2.0:
                    self._exp_val = zeta(a - 1.0) / zeta(a)
                    
                # use some valid value when exp. value is infinite
                elif self._min is not None:
                    self._exp_val = self._min
                else:
                    self._exp_val = 0.0
                    
                self._distr = 'zipf'
                self.draw = self._draw_zipf
                self._a = a
                    
            elif self._params[0] == 'negative_binomial':
                
                r = self._params[1]
                if not len(self._params) >= 3:
                    raise ValueError('negative binomial distr. requires 2 '\
                                     'parameters')
                q = self._params[2]
                    
                if not isinstance(r, int) or r < 1:
                    raise ValueError('parameter r must be an int and >0')
                if not isinstance(q, float) or q <= 0 or q >= 1.0:
                    raise ValueError('parameter q must be a float >0 and <1')
                    
                self._distr = 'negative_binomial'
                self.draw = self._draw_neg_bin
                self._exp_val = r * q / (1 - q)
                self._r = r
                self._q = q
                
            else:
                raise ValueError("length distribution '{}' is not "\
                                 "supported".format(self._params[0]))
            
        else:
            raise ValueError('could not parse distribution')
        
        self._exp_val += self._shift
        
        # expected value must be between min. and max.
        if ((self._min is not None and self._exp_val < self._min) or
            (self._max is not None and self._exp_val > self._max)):
            raise ValueError('expected value must be >= minimum and'\
                             '<= maximum')
            
                
    def _draw_constant(self):
        
        return self._exp_val
    
    
    def _draw_uniform(self):
        
        while True:
            x = np.random.uniform(low=self._a, high=self._b) + self._shift
            
            if self._discrete:
                x = round(x)
                
            if ((not self._min or x >= self.min) and
                (not self._max or x <= self._max)):
                return x
            
    
    def _draw_discrete_uniform(self):
        
        while True:
            x = np.random.randint(self._a, high=self._b) + self._shift
                
            if ((not self._min or x >= self.min) and
                (not self._max or x <= self._max)):
                return x
                
    
    def _draw_gamma(self):
        
        while True:
            x = np.random.gamma(self._shape, scale=self._scale) + self._shift
            
            if self._discrete:
                x = round(x)
                
            if ((not self._min or x >= self.min) and
                (not self._max or x <= self._max)):
                return x
            
            
    def _draw_exponential(self):
        
        while True:
            if self._rate == 0.0:
                x = float('inf')
            else:
                x = np.random.exponential(scale=1/self._rate) + self._shift
            
            if self._discrete:
                x = round(x)
                
            if ((not self._min or x >= self.min) and
                (not self._max or x <= self._max)):
                return x
            
    
    def _draw_zipf(self):
        
        while True:
            x = np.random.zipf(self._a) + self._shift
                
            if ((not self._min or x >= self.min) and
                (not self._max or x <= self._max)):
                return x
            
    
    def _draw_neg_bin(self):
        
        while True:
            x = np.random.negative_binomial(self._r, 1 - self._q) + self._shift
                
            if ((not self._min or x >= self.min) and
                (not self._max or x <= self._max)):
                return x