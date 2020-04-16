# -*- coding: utf-8 -*-

from math import log

"""Calculation of sequence distances for various models.

References:
    Ziheng Yang (2006). Computational Molecular Evolution.
    Oxford Series in Ecology and Evolution.
"""


def p_distance(seq1, seq2, exclude_gaps=True):
    
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


def JC69_distance(p, amino_acid=False):
    """Jukes Cantor 1969 distance."""
    
    a = 3 if not amino_acid else 19     # numerator
    b = 4 if not amino_acid else 20     # denominator
    
    if p >= a / b:
        d = float('nan')
    
    else:
        d = - (a/b) * log(1 - (b/a) * p)
        
    return d


def JC69_distance_var(p, n, amino_acid=False):
    """Jukes Cantor 1969 distance variance."""
    
    a = 3 if not amino_acid else 19     # numerator
    b = 4 if not amino_acid else 20     # denominator
    
    if p >= a / b:
        var = float('nan')
    
    else:
        var = p * (1 - p) / (n * (1 - (b/a) * p) ** 2)
        
    return var


def IV_proportions(seq1, seq2):
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


def K80_distance(S, V):
    """Kimura 1980 distance and transition/transversion ratio."""
    
    try:
        d = - 0.5 * log(1 - 2 * S - V) - 0.25 * log(1 - 2 * V)
        kappa = 2 * log(1 - 2 * S - V) / log(1 - 2 * V) - 1
        
        return d, kappa
    
    except ValueError:
        return float('nan'), float('nan')
    
    except ZeroDivisionError:
        return float('nan'), float('nan')


def K80_distance_var(S, V, n):
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