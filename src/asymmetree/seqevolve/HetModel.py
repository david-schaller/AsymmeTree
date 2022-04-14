# -*- coding: utf-8 -*-

"""Heterogeneity models for sequence evolution.

References
----------
.. [1] Z. Yang.
   Computational molecular evolution.
   Oxford series in ecology and evolution. Oxford University Press, 2006.
   ISBN 978-0-19-856699-1 978-0-19-856702-8.
"""

import numpy as np

from asymmetree.seqevolve.EvolvingSequence import State


__author__ = 'David Schaller'


class HetModel:
    """Heterogeneity model.
    
    Supports rate heterogeneity based on a Gamma distribution ('+Gamma'-model)
    as well as invariant sites ('+I'-model), see [1].
    
    References
    ----------
    .. [1] Z. Yang.
       Computational molecular evolution.
       Oxford series in ecology and evolution. Oxford University Press, 2006.
       ISBN 978-0-19-856699-1 978-0-19-856702-8.
    """
    
    def __init__(self, alpha, classes=5, sitewise=False, invariant=0.0):
        """
        Parameters
        ----------
        alpha : float
            Parameter of the Gamma distribution (with mean 1) from which the
            rate factors are drawn.
        classes : int, optional
            Number of classes such that sites in the same class share the same
            rate factor. The default is 5.
        sitewise : bool, optional
            If True, ignore the 'classes' attribute and treat each site as its
            own class. The default is False.
        invariant : float, optional
            Expected proportion of invariant sites, i.e., sites where no
            mutations happen. The default is 0.0.
        """
        
        if not isinstance(alpha, float) or alpha <= 0.0:
            raise ValueError('heterogeneity parameter alpha must be a float > 0.0')
        self._alpha = alpha
        
        if not isinstance(classes, int) or classes <= 0:
            raise ValueError('number of classes must be an int >0')
        self.classes = classes
        
        self.sitewise = sitewise
        
        if not isinstance(invariant, float) or invariant < 0.0 or invariant > 1.0:
            raise ValueError('proportion of invariant sites must be in [0.0, 1.0]')
        self._invariant = invariant
        
        if not self.sitewise:
            self._initilize_classes()
        
    
    def _initilize_classes(self):
        
        if self.classes == 1:
            self._class_rates = [1.0]
            
        else:
            self._class_rates = np.random.gamma(self._alpha,
                                                scale=1/self._alpha,
                                                size=self.classes)
        
        
    def assign(self, sequence, exclude_inherited=True):
        """Assign rate classes and rate factors tothe sites of a sequence.
        
        Parameters
        ----------
        sequence : asymmetree.seqevolve.EvolvingSequence.EvoSeq
            The sequence.
        exclude_inherited : bool, optional
            If True, do not assign new classes and rates to inherited sites,
            i.e., only if the sequence corresponds to the root of a tree along
            which is simulated or if new sites are added as a result of an
            insertion event.
        """
        
        n = len(sequence)
        if exclude_inherited:
            n -= sequence.count_status(State.INHERITED)
        
        rate_classes, rate_factors = self._draw(n)
        
        pos = 0
        for site in sequence:
            
            if not exclude_inherited or site.status != State.INHERITED:
                site.rate_class = rate_classes[pos] if rate_classes[pos] is not None else site.site_id
                site.rate_factor = rate_factors[pos]
                pos += 1
                    
                    
    def _draw(self, n):
        
        # mode 1: sitewise heterogeneity
        if self.sitewise:
            
            rate_classes = [None for _ in range(n)]
            rate_factors = np.random.gamma(self._alpha, scale=1/self._alpha, size=n)
        
        # mode 2: one or multiple  classes
        else:
            
            rate_classes = np.random.randint(self.classes, size=n).tolist()
            rate_factors = np.asarray([self._class_rates[c] for c in rate_classes])
            
        if self._invariant:
            
            drawn_variable = np.random.random(n) > self._invariant
            
            for i in range(n):
                if not drawn_variable[i]:
                    rate_classes[i] = 'invariant'
            
            rate_factors *= drawn_variable
        
        return rate_classes, rate_factors