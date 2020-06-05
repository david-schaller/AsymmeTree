# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize

from asymmetree.seqevolve import SubstModel

"""Calculation of sequence distances for various models.

References:
    Ziheng Yang (2006). Computational Molecular Evolution.
    Oxford Series in Ecology and Evolution.
"""


__author__ = 'David Schaller'


MIN_LENGTH = 1e-6
MAX_LENGTH = 20.0
#MAXITER = 100


def maximum_likelihood_distance(seq1, seq2,
                                subst_model=None,
                                model_type=None, model_name=None):
    
    if subst_model is None:
        
        if model_type is None or model_name is None:
            raise ValueError('no substitution model specified')
        else:
            subst_model = SubstModel(model_type, model_name)
            
    seqs = _to_indices(seq1, seq2, subst_model)
    
    x0 = _initial_guess(seqs, subst_model.model_type)
    
    # distance 0.0 if sequences are equal
    if x0 == 0.0:
        return x0
    
    opt_result = scipy.optimize.minimize(_likelihood, x0, args=(seqs, subst_model),
                                         bounds=((MIN_LENGTH, None),),)
    
    if opt_result.success:
        return opt_result.x
    else:
        return float('nan')
    

def _initial_guess(seqs, model_type):
    
    p = np.sum(seqs[0] != seqs[1]) / seqs.shape[1]
    
    x0 = _JC69_transform(p, amino_acid=(model_type=='a'))
    
    x0 = x0 if x0 != float('nan') else MAX_LENGTH / 2
    
    return x0
    

def _likelihood(x, seqs, subst_model):
    
    P = subst_model.transition_prob_matrix(x)
    
    # - log-likelihood
    return - np.sum( np.log( P[seqs[0], seqs[1]] ))
            
            
def _to_indices(seq1, seq2, subst_model):
    
    if len(seq1) != len(seq2):
        raise ValueError("unequal sequence lengths: {} and {}".format(len(seq1), len(seq2)))

    seq1_indeces, seq2_indeces = [], []
    
    for i in range(len(seq1)):
        
        if seq1[i] != '-' and seq2[i] != '-':
            
            seq1_indeces.append( subst_model.alphabet_dict[seq1[i]] )
            seq2_indeces.append( subst_model.alphabet_dict[seq2[i]] )

    return np.asarray([seq1_indeces, seq2_indeces], dtype=int)     


def p_distance(seq1, seq2, exclude_gaps=True):
    """Calculate the p-distance of two aligned sequences.
    
    = length-normalized Hamming distance.
    
    Keyword arguments:
        exclude_gaps - ignore columns with a gap in one sequence,
           gaps in both sequences are always ignored; default=True.
    """
    
    if len(seq1) != len(seq2):
        raise ValueError("unequal sequence lengths: {} and {}".format(len(seq1), len(seq2)))
        
    diffs, valid_columns = 0, 0
    
    for i in range(len(seq1)):
        
        if seq1[i] != '-' and seq2[i] != '-':
            
            valid_columns += 1
            if seq1[i] != seq2[i]:
                diffs += 1
        
        elif not exclude_gaps and (seq1[i] != '-' or seq2[i] != '-'):
            
            valid_columns += 1
            diffs += 1
            
    p = diffs / valid_columns if valid_columns > 0 else float('inf')
    
    return p, valid_columns


def JC69_distance(seq1, seq2, exclude_gaps=True, amino_acid=False, 
                  variance=False):
    """Jukes Cantor 1969 distance."""
    
    p, valid_columns = p_distance(seq1, seq2, exclude_gaps=exclude_gaps)
    
    d = _JC69_transform(p, amino_acid=amino_acid)
    
    if not variance:
        return d
    
    else:
        var_d = _JC69_distance_var(p, valid_columns,
                                   amino_acid=amino_acid)
        return d, var_d


def _JC69_transform(p, amino_acid=False):
    """Jukes Cantor 1969 transformation of p-distance."""
    
    a = 3 if not amino_acid else 19     # numerator
    b = 4 if not amino_acid else 20     # denominator
    
    if p >= a / b:
        d = float('nan')
    
    else:
        d = - (a/b) * np.log(1 - (b/a) * p)
        
    return d


def _JC69_distance_var(p, n, amino_acid=False):
    """Jukes Cantor 1969 distance variance."""
    
    a = 3 if not amino_acid else 19     # numerator
    b = 4 if not amino_acid else 20     # denominator
    
    if p >= a / b:
        var = float('nan')
    
    else:
        var = p * (1 - p) / (n * (1 - (b/a) * p) ** 2)
        
    return var


def K80_distance(seq1, seq2, variance=False):
    """Kimura 1980 distance and transition/transversion ratio."""
    
    S, V, valid_columns = _IV_proportions(seq1, seq2)
    
    d, kappa = _K80_transform(S, V)
    
    if not variance:
        return d, kappa
    
    else:
        var_d = _K80_distance_var(S, V, valid_columns)
        return d, kappa, var_d
    

def _IV_proportions(seq1, seq2):
    """Computes the proportions of transitions and transversions."""
    
    purines = {'A', 'a', 'G', 'g', 'R', 'r'}
    pyrimidines = {'C', 'c', 'T', 't', 'U', 'u', 'Y', 'y'}
    
    if len(seq1) != len(seq2):
        raise ValueError("unequal sequence lengths: {} and {}".format(len(seq1), len(seq2)))
        
    transitions, transversions, valid_columns = 0, 0, 0
    
    for i in range(len(seq1)):
        
        if seq1[i] != '-' and seq2[i] != '-':
            
            valid_columns += 1
            
            if seq1[i] == seq2[i]:
                continue
            
            elif ((seq1[i] in purines and seq2[i] in purines) or
                  (seq1[i] in pyrimidines and seq2[i] in pyrimidines)):
                transitions += 1
            
            elif ((seq1[i] in purines and seq2[i] in pyrimidines) or
                  (seq1[i] in pyrimidines and seq2[i] in purines)):
                transversions += 1
                
            
    S = transitions / valid_columns if valid_columns > 0 else float('nan')
    V = transversions / valid_columns if valid_columns > 0 else float('nan')
    
    return S, V, valid_columns


def _K80_transform(S, V):
    """Kimura 1980 distance and transition/transversion ratio."""
    
    a1 = 1.0 - 2.0 * S - V
    a2 = 1.0 - 2.0 * V
    
    if a1 <= 0.0 or a2 <= 0.0:
        return float('nan'), float('nan')
    
    else:
        d = - 0.5 * np.log(a1) - 0.25 * np.log(a2)
        kappa = 2.0 * np.log(a1) / np.log(a2) - 1.0 if a2 < 1.0 else float('nan')
        
        return d, kappa


def _K80_distance_var(S, V, n):
    """Kimura 1980 distance variance."""
    
    a1 = 1.0 - 2.0 * S - V
    a2 = 1.0 - 2.0 * V
    
    if a1 <= 0.0 or a2 <= 0.0:
        return float('nan')
    
    else:
        a = 1.0 / a1
        b = 0.5 * ( 1.0 / a1 + 1.0 / a2 )
        
        return (a**2 * S + b**2 * V - (a * S + b * V)**2) / n
