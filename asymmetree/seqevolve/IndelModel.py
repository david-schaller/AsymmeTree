# -*- coding: utf-8 -*-

"""Indel models for sequence evolution.

References
----------
.. [1] Z. Yang.
   Computational molecular evolution.
   Oxford series in ecology and evolution. Oxford University Press, 2006.
   ISBN 978-0-19-856699-1 978-0-19-856702-8.
.. [2] R. A. Cartwright.
   DNA assembly with gaps (Dawg): Simulating sequence evolution.
   In: Bioinformatics, 21(Suppl 3):iii31-iii38, November 2005.
   doi:10.1093/bioinformatics/bti1200.
"""

from asymmetree.tools.Sampling import Sampler


__author__ = 'David Schaller'


class IndelModel:
    """Indel model"""
    
    def __init__(self, insertion_rate, deletion_rate,
                 length_distr=('zipf', 1.821),        # (Chang and Benner 2004)
                 min_length=1,
                 max_length=None,
                 **kwargs):
        """
        Parameters
        ----------
        insertition_rate : float
            Insertion rate w.r.t. units of evolutionary distance.
        deletion_rate : float
            Deletion rate w.r.t. units of evolutionary distance.
        length_distr : tuple, optional
            Distribution of the length of the indels, see documentation for
            available options. The default is a Zipf distribution as observed
            for empirical data in [1].
        min_length : int, optional
            The minimal length of an indel, default is 1.
        max_length : int, optional
            The maximal length of an indel. The default is None in which case
            the distribution is not truncated (on the upper side).
        
        References
        ----------
        .. [1] M. S. Chang and S. A. Benner.
           Empirical Analysis of Protein Insertions and Deletions Determining
           Parameters for the Correct Placement of Gaps in Protein Sequence
           Alignments. 
           In: Journal of Molecular Biology, 341(2):617-631, August 2004.
           doi:10.1016/j.jmb.2004.05.045.
        """
        
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
        and the length of the sequence.
        
        Parameters
        ----------
        seq_length : int
            The length of the sequence.
        
        Returns
        -------
        tuple of two floats
            The total rate for an insertation and deletion event, see [1].
            
        References
        ----------
        .. [1] R. A. Cartwright.
           DNA assembly with gaps (Dawg): Simulating sequence evolution.
           In: Bioinformatics, 21(Suppl 3):iii31-iii38, November 2005.
           doi:10.1093/bioinformatics/bti1200.
        """
        
        # expected value may be infinite for zipf distribution, in which case
        # the minimal length is used, see 'Sampler'
        return ( (seq_length + 1) * self._ins_rate,
                 (seq_length + self.sampler._exp_val - 1) * self._del_rate )
        
    
    def draw_length(self):
        """Draw the length for an indel.
        
        Returns
        -------
        int
            The drawn length for the indel.
        """
        
        return self.sampler()
    