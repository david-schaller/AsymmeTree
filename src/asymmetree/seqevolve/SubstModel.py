# -*- coding: utf-8 -*-

"""Substitution models for nucleotide and amino acid evolution.

References
----------
.. [1] Z. Yang.
   Computational molecular evolution.
   Oxford series in ecology and evolution. Oxford University Press, 2006.
   ISBN 978-0-19-856699-1 978-0-19-856702-8.
"""

import numpy as np
from scipy import linalg

from asymmetree.seqevolve.EmpiricalModels import empirical_models
from asymmetree.file_io.SubstModelIO import parse_paml


__author__ = 'David Schaller'


class SubstModel:
    """Substitution model for nucleotide or amino acid evolution."""
    
    nuc_models = {'JC69', 'K80', 'GTR', 'CUSTOM'}
    aa_models = {'JC69', 'CUSTOM'}
    
    nucleotides = 'ACGT'
    amino_acids = 'ARNDCQEGHILKMFPSTWYV'
    
    
    def __init__(self, model_type, model_name,
                 **kwargs):
        """
        Parameters
        ----------
        model_type : str
            Available options are 'n' for nucleotide and 'a' for amino acid
            sequences.
        model_name : str
            Availbale are e.g. 'JC69' (nuc./a.a.), 'K80' (nuc.), 'GTR' (nuc.),
            'WAG' (a.a.), 'JTT' (a.a.), 'BLOSUM62' (a.a.), 'LG' (a.a.), and
            'DAYHOFF' (a.a.).
            For the option 'CUSTOM', the keyword parameter 'filename' has to
            be specified with a path to a custom model in paml format.
        kwargs : optional
            Some models require additional parameters, e.g. 'K80' requires the
            parameter 'kappa' and 'GTR' requires 'abcdef' and 'f', see
            documention for details.
        
        Raises
        ------
        ValueError
            If an invalid model was specified.
        """
        
        self.model_type = model_type.lower()
        self.model_name = model_name.upper()
        
        self._params = kwargs
        
        if (self.model_type in ('n', 'nuc', 'nucleotide') and
            self.model_name in SubstModel.nuc_models):
            
            self.model_type = 'n'
            self.alphabet = SubstModel.nucleotides
            
        elif (self.model_type in ('a', 'aa', 'amino', 'aminoacid', 'protein') and
              (self.model_name in SubstModel.aa_models or 
               self.model_name in empirical_models)):
            
            self.model_type = 'a'
            self.alphabet = SubstModel.amino_acids
                
        else:
            raise ValueError("model '{}', '{}' is not "\
                             "available".format(model_type, model_name))
        
        self.alphabet_dict = {item: index for index, item in enumerate(self.alphabet)}
        
        self._load_exchangeability_and_freqs()
        self._build_rate_matrix()
        
    
    def _load_exchangeability_and_freqs(self):
        """Load the exchangeability matrix S and the stationary frequencies pi."""
        
        # a custom model (via a paml file) was specified
        if self.model_name == 'CUSTOM':
            
            if self._params is None or 'filename' not in self._params:
                    raise ValueError("custom model requires the parameter "\
                                     "'filename'")
            
            S, freqs = parse_paml(self._params['filename'],
                                  model_type=self.model_type)
            self.S, self.freqs = np.asarray(S), np.asarray(freqs)
        
        # non-empirical nucleotide models
        elif self.model_type == 'n':
            
            if self.model_name == 'JC69':
                self.S, self.freqs = _JC69_nuc()
                
            elif self.model_name == 'K80':
                
                if self._params is None or 'kappa' not in self._params:
                    raise ValueError("model 'K80' requires the parameter 'kappa'")
                
                self.S, self.freqs = _K80_nuc(self._params['kappa'])
            
            elif self.model_name == 'GTR':
                
                if self._params is None or 'abcdef' not in self._params or 'f' not in self._params:
                    raise ValueError("model 'GTR' requires the parameters "\
                                     "'abcdef' and 'f'")
                
                self.S, self.freqs = _GTR_nuc(self._params['abcdef'],
                                              self._params['f'])
        
        # non-empirical and empirical amino acid models
        elif self.model_type == 'a':
            
            if self.model_name == 'JC69':
                self.S, self.freqs = _JC69_aa()
                
            elif self.model_name in empirical_models:
                self.S, self.freqs = empirical_models[self.model_name]()
                
        # make sure stationary frequencies sum to 1
        self.freqs /= np.sum(self.freqs)
        
        # compute cumulative frequencies for evolver
        self.freqs_cumulative = np.cumsum(self.freqs)
    
    
    def _build_rate_matrix(self):
        """Compute the rate matrix Q."""
    
        # rate matrix from exchangeability matrix and stationary frequencies
        Q = self.S @ np.diag(self.freqs)
        
        # compute diagonals, rows sum to 0
        for i in range(Q.shape[0]):
            Q[i, i] = - (np.sum(Q[i, :]) - Q[i, i])
            
        # normalize matrix
        k = - np.dot(self.freqs, np.diagonal(Q))
        Q /= k
        
        self.Q = Q
    
    
    def eigensystem(self):
        
        if not hasattr(self, 'eigenvals'):
            self.eigenvals, self.U, self.U_inv = diagonalize(self.Q, self.freqs)
        
        return self.eigenvals, self.U, self.U_inv
    
    
    def transition_prob_matrix(self, t):
        """Calculate the transition probability matrix P(t).
        
        P(t) = U x e^(Lambda*t) x U^(-1).
        
        Parameters
        ----------
        t : float
            The time / evolutinary distance.
        
        Returns
        -------
        2-dimensional numpy array
        """
        
        # ensure that eigensystem has been computed
        self.eigensystem()
        
        # first multiplication element-wise, since corresponding matrix
        # only has non-zero entries on the diagonal
        
        return (self.U * np.exp(self.eigenvals * t))  @  self.U_inv
        
    
    def to_indices(self, sequence):
        """Convert a sequence into a list of indices the represent the positions
        of the letters in the alphabet of the model.
        
        Parameters
        ----------
        sequence : iterable
            The sequence to be converted.
        
        Returns
        -------
        list
            The list of indices.
        
        Raises
        ------
        ValueError
            If the sequence contains a letter that is not in the alphabet of
            the model.
        """
        
        try:
            result = [self.alphabet_dict[letter] for letter in sequence]
        except KeyError:
            raise ValueError('invalid sequence for the specified model')
            
        return result
    
    
    def to_sequence(self, evoseq):
        """Construct the str representation of a sequence.
        
        Parameters
        ----------
        evoseq : asymmetree.seqevolve.EvolvingSequence.EvoSeq
            The sequence to be converted.
        
        Returns
        -------
        str
            The nucleotide or amino acid sequence depending on the model.
        """
        
        return ''.join(self.alphabet[x._value] for x in evoseq)
    

# --------------------------------------------------------------------------
#                           Module functions
# --------------------------------------------------------------------------
    
    
def diagonalize(Q, freqs):
    """Diagonalize a matrix.
    
    Parameters
    ----------
    Q : 2-dimensional numpy array
        Rate matrix.
    freqs : 1-dimensional numpy array
        stationary frequencies
    
    Returns
    -------
    1-dimensional numpy array
        Eigenvalues.
    2-dimensional numpy array
        Diagonalized matrix.
    2-dimensional numpy array
        Inverse of the diagonalized matrix.
    """
        
    # matrix is already symmetric
    if np.allclose(Q, Q.T):
        
        eigenvals, U = linalg.eigh(Q)
        U_inv = linalg.inv(U)
    
    # matrix is not symmetric but model is time-reversible
    else:
    
        Phi = np.diag(np.sqrt(freqs))
        Phi_inv = linalg.inv(Phi)
        
        B = Phi @ Q @ Phi_inv
        
        eigenvals, R = linalg.eigh(B)

        U = Phi_inv @ R
        U_inv = linalg.inv(R) @ Phi
    
    return eigenvals, U, U_inv
    
    
# --------------------------------------------------------------------------
#                           NUCLEOTIDE MODELS
# --------------------------------------------------------------------------

def _JC69_nuc():
    
    S = np.array([[0, 1, 1, 1],
                  [1, 0, 1, 1],
                  [1, 1, 0, 1],
                  [1, 1, 1, 0]])
    
    freqs = np.array([0.25, 0.25, 0.25, 0.25])
    
    return S, freqs


def _K80_nuc(kappa):
    
    S = np.array([[0,     1,     kappa, 1    ],
                  [1,     0,     1,     kappa],
                  [kappa, 1,     0,     1    ],
                  [1,     kappa, 1,     0    ]])
    
    freqs = np.array([0.25, 0.25, 0.25, 0.25])
    
    return S, freqs


def _GTR_nuc(abcdef, f):
    """Generalized time-reversible model.
    
    Parameterization as in e.g. used in PAML and ALF:
    a:    C <--> T
    b:    A <--> T
    c:    G <--> T
    d:    A <--> C  
    e:    C <--> G
    f:    A <--> G
    """
    
    a, b, c, d, e, f = abcdef
    
    S = np.array([[0,  d,  f,  b],
                  [d,  0,  e,  a],
                  [f,  e,  0,  c],
                  [b,  a,  c,  0]])
    
    freqs = np.asarray(abcdef)
    freqs /= np.sum(freqs)
    
    return S, freqs


# --------------------------------------------------------------------------
#                          AMINO ACID MODELS
# --------------------------------------------------------------------------

def _JC69_aa():
    
    S = np.ones((20, 20))
    S -= np.identity(20)
    
    freqs = np.full((20,), 1/20)
    
    return S, freqs