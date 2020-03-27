# -*- coding: utf-8 -*-

import numpy as np


def build_matrix_Q(S, freqs):
    
    Q = S @ np.diag(freqs)
    
    # compute diagonals, rows sum to 0
    for i in range(Q.shape[0]):
        Q[i, i] = - (np.sum(Q[i, :]) - Q[i, i])
        
    # normalize matrix
    k = - np.dot(freqs, np.diagonal(Q))
    Q /= k
    
    return Q


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


if __name__ == "__main__":
    
    S, freqs = _JC69_nuc()
    print(S)
    print(freqs)
    
    Q = build_matrix_Q(S, freqs)
    print(Q)
    
    eigvals, U, U_inv = diagonalize(Q, freqs)