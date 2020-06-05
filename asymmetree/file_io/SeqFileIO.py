# -*- coding: utf-8 -*-

from asymmetree import TreeNode


__author__ = 'David Schaller'


# ------------------------------------------------------------------------------
#                          Auxiliary functions
# ------------------------------------------------------------------------------

def labeled_sequences(sequences):
    """Converts dictionaries to lists and keys into strings.
    
    If the keys are 'TreeNode' instance, then their label is taken. Also works
    if 'sequences' already is a list (of tuples).
    """
    
    result = []
    
    if isinstance(sequences, dict):
        sequences = sequences.items()
    
    for key, sequence in sequences:
        
        if isinstance(key, TreeNode):
            result.append( (key.label, sequence) )
        else:
            result.append( (str(key),  sequence) )
            
    return result


# ------------------------------------------------------------------------------
#                               Fasta files
# ------------------------------------------------------------------------------

def write_fasta(filename, sequences):
    
    sequences = labeled_sequences(sequences)
    
    with open(filename, 'w') as f:
        
        for label, seq in sequences:
            
            f.write('>{}\n'.format(label))
            pos = 0
            while pos < len(seq):
                f.write(seq[pos:min(pos+80, len(seq))])
                pos += 80
                f.write('\n')


# ------------------------------------------------------------------------------
#                       Multiple Sequence Alignments
# ------------------------------------------------------------------------------

def write_alignment(filename, alignment, alignment_format='phylip'):
    
    alignment = labeled_sequences(alignment)
    
    with open(filename, 'w') as f:
        
        if alignment_format == 'phylip':
            _write_phylip(f, alignment)
    
        elif alignment_format == 'clustal':
            _write_clustal(f, alignment)
        
        elif alignment_format == 'pretty':
            _write_pretty(f, alignment)
            
        else:
            raise ValueError("alignment format '{}' is not "\
                             "available".format(alignment_format))
            

def _check_alignment(alignment):
    
    max_length = 0                          # maximal label length
    seq_length = None                       # length of the aligned sequences
    for label, seq in alignment:
        
        if len(label) > max_length:
            max_length = len(label)
            
        if seq_length is None:
            seq_length = len(seq)
        elif seq_length != len(seq):
            raise ValueError('aligned sequences must have the same length')
            
    return max_length, seq_length
            
            
def _write_phylip(f, alignment):
    
    max_length, seq_length = _check_alignment(alignment)
    max_length = max(max_length, 9)
    
    
    f.write('  {} {} i'.format(len(alignment), seq_length))
    
    format_str = '\n{:' + str(max_length+1) + '}'
    current = 0
    
    while current < seq_length:
        
        end = min(seq_length, current+50)
        
        for label, seq in alignment:
            
            if current == 0:
                f.write(format_str.format(label))
            else:
                f.write(format_str.format(''))
                
            f.write(' '.join( seq[i:min(i+10,end)] for i in range(current, end, 10) ))
                
        if end != seq_length:
            f.write('\n')
        
        current += 50


def _write_clustal(f, alignment):
    
    max_length, seq_length = _check_alignment(alignment)
    
    f.write('CLUSTAL W (1.8) multiple sequence alignment\n\n')
    
    format_str = '\n{:' + str(max_length+4) + '}{}'
    current = 0
    
    while current < seq_length:
        
        end = min(seq_length, current+60)
        
        for label, seq in alignment:
            f.write(format_str.format(label, seq[current:end]))
                
        if end != seq_length:
            f.write('\n\n')
        
        current += 60
    
        
def _write_pretty(f, alignment):
    
    _, seq_length = _check_alignment(alignment)
    
    f.write('  1          11         21         31         41       50\n')
    f.write('  |          |          |          |          |        |')
    
    current = 0
    
    while current < seq_length:
        
        end = min(seq_length, current+50)
        
        for label, seq in alignment:
            
            count = end - current - seq[current:end].count('-')
            seq_string = ' '.join( seq[i:min(i+10,end)] for i in range(current, end, 10) )
            f.write('\n  {:54}{:>6} {}'.format(seq_string, count, label))
                
        if end != seq_length:
            f.write('\n\n\n')
        
        current += 50