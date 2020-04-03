# -*- coding: utf-8 -*-

import numpy as np

from asymmetree.seqevolve.EmpiricalModels import empirical_models


class SubstModel:
    
    
    nuc_models = {'JC69', 'K80'}
    aa_models = {'JC69'}
    
    nucleotides = 'ACGT'
    amino_acids = 'ARNDCQEGHILKMFPSTWYV'
    
    
    def __init__(self, model_type, model_name):
        
        self.model_type = model_type.lower()
        self.model_name = model_name.upper()
        
        if (self.model_type in ('nuc', 'nucleotide') and
            self.model_name in SubstModel.nuc_models):
            
            self.model_type = 'nuc'
            self.alphabet = SubstModel.nucleotides
            
        elif (self.model_type in ('aa', 'amino', 'aminoacid', 'protein') and
              (self.model_name in SubstModel.aa_models or 
               self.model_name in empirical_models)):
            
            self.model_type = 'aa'
            self.alphabet = SubstModel.amino_acids
                
        else:
            raise ValueError("Model '{}', '{}' is not available!".format(model_type, model_name))
        
        self.alphabet_dict = {item: index for index, item in enumerate(self.alphabet)}
        
        self._load_exchangeability_and_freqs()
        self._build_rate_matrix()
        
    
    def _load_exchangeability_and_freqs(self):
        """Load the exchangeability matrix S and the stationary frequencies pi."""
        
        if self.model_type == 'nuc':
            
            if self.model_name == 'JC69':
                self.S, self.freqs = _JC69_nuc()
        
        elif self.model_type == 'aa':
            
            if self.model_name == 'JC69':
                self.S, self.freqs = _JC69_aa()
                
            elif self.model_name in empirical_models:
                self.S, self.freqs = empirical_models[self.model_name]()
                
        # make sure stationary frequencies sum to 1
        self.freqs /= np.sum(self.freqs)
    
    
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
            self.eigvals, self.U, self.U_inv = diagonalize(self.Q, self.freqs)
        
        return self.eigvals, self.U, self.U_inv
        
    
    def to_indices(self, sequence):
        
        try:
            result = [self.alphabet_dict[letter] for letter in sequence]
        except KeyError:
            raise ValueError("Invalid sequence for the specified model!")
            
        return result
    
    
    def to_sequence(self, evoseq):
        
        return "".join(self.alphabet[x._value] for x in evoseq)
    

# --------------------------------------------------------------------------
#                           Module functions
# --------------------------------------------------------------------------
    
    
def diagonalize(Q, freqs):
        
    # matrix is already symmetric
    if np.allclose(Q, Q.T):
        
        eigvals, U = np.linalg.eigh(Q)
        U_inv = np.linalg.inv(U)
    
    # matrix is not symmetric but model is time-reversible
    else:
    
        Phi = np.diag(np.sqrt(freqs))
        Phi_inv = np.linalg.inv(Phi)
        
        B = Phi @ Q @ Phi_inv
        
        eigvals, R = np.linalg.eigh(B)

        U = Phi_inv @ R
        U_inv = np.linalg.inv(R) @ Phi
    
    return eigvals, U, U_inv
    
    
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


# --------------------------------------------------------------------------
#                          AMINO ACID MODELS
# --------------------------------------------------------------------------

def _JC69_aa():
    
    S = np.ones((20, 20))
    S -= np.identity(20)
    
    freqs = np.full((20,), 1/20)
    
    return S, freqs