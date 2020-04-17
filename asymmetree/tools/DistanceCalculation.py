# -*- coding: utf-8 -*-

from math import log
import numpy as np
from scipy.optimize import minimize

from asymmetree.seqevolve import SubstModel

"""Calculation of sequence distances for various models.

References:
    Ziheng Yang (2006). Computational Molecular Evolution.
    Oxford Series in Ecology and Evolution.
"""


def maximum_likelihood_distance(seq1, seq2,
                                subst_model=None,
                                model_type=None, model_name=None):
    
    if subst_model is None:
        
        if model_type is None or model_name is None:
            raise ValueError('No substitution model specified!')
        else:
            subst_model = SubstModel(model_type, model_name)
            
    seqs = _to_indices(seq1, seq2, subst_model)
    
    # use Jukes-Cantor as starting guess
    if subst_model.model_type == 'n':
        x0 = JC69_distance(seq1, seq2, amino_acid=False)
    elif subst_model.model_type == 'a':
        x0 = JC69_distance(seq1, seq2, amino_acid=True)
    x0 = x0 if x0 != float('nan') else 5.0
        
    opt_result = minimize(_likelihood, x0, args=(seqs, subst_model),
                          bounds=[(0.0, None)],)
    
    if opt_result.success:
        return opt_result.x
    else:
        return float('nan')
    

def _likelihood(x, seqs, subst_model):
    
    P = subst_model.transition_prob_matrix(x)
    
    log_likelihood = 0
    for k in range(seqs.shape[1]):
        p = P[ seqs[0][k], seqs[1][k] ]
        log_likelihood += log(p) if p > 0.0 else float('inf')
        
    return -log_likelihood
            
            
def _to_indices(seq1, seq2, subst_model):
    
    if len(seq1) != len(seq2):
        raise ValueError("Unequal sequence lengths: {} and {}!".format(len(seq1), len(seq2)))

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
        raise ValueError("Unequal sequence lengths: {} and {}!".format(len(seq1), len(seq2)))
        
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
        d = - (a/b) * log(1 - (b/a) * p)
        
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
        raise ValueError("Unequal sequence lengths: {} and {}!".format(len(seq1), len(seq2)))
        
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
    
    try:
        d = - 0.5 * log(1 - 2 * S - V) - 0.25 * log(1 - 2 * V)
        kappa = 2 * log(1 - 2 * S - V) / log(1 - 2 * V) - 1
        
        return d, kappa
    
    except ValueError:
        return float('nan'), float('nan')
    
    except ZeroDivisionError:
        return float('nan'), float('nan')


def _K80_distance_var(S, V, n):
    """Kimura 1980 distance variance."""
    
    try:
        a = 1 / (1 - 2 * S - V)
        b = 0.5 * ( 1 / (1 - 2 * S - V) + 1 / (1 - 2 * V) )
        
        var_d = (a**2 * S + b**2 * V - (a * S + b * V)**2) / n
        
        return var_d
    
    except ValueError:
        return float('nan')
    
    except ZeroDivisionError:
        return float('nan')
