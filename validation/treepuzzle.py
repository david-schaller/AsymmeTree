# -*- coding: utf-8 -*-

import subprocess

import numpy as np


def parse_distance_matrix(filename):
    
    with open(filename, 'r') as f:
        
        header = f.readline()
        number = int(header.split()[0])
        
        D = np.zeros((number, number))
        labels = []
        
        line = f.readline()
        
        for i in range(number):
            
            row = line.strip()
            line = f.readline()
            while line and line.startswith(' '):
                row += line.rstrip()
                line = f.readline()
            
            row = row.split()
            
            labels.append(row[0])
            
            for j in range(number):
                
                D[i, j] = float(row[j+1])
                
    return D, labels


def rearrange_matrix(D, labels, target_labels):
    
    n = len(labels)
    
    if set(labels) != set(target_labels):
        raise ValueError("Label sets must be equal!")
    
    mapping = []
    for i in range(n):
        mapping.append( labels.index(target_labels[i]) )
    
    D_new = np.zeros_like(D)
    
    for i in range(n):
        for j in range(n):
            
            D_new[i, j] = D[mapping[i], mapping[j]]
    
    return D_new


def calc_distances(infile, model='WAG'):
    
    models = ['WAG', 'DAYHOFF', 'JTT', 'MTREV24', 'BLOSUM62', 'VT']
    
    args = ['puzzle', infile]
    
    communication = b'k\nk\nk\n'
    
    if model.upper() in models:
        for i in range(models.index(model.upper())):
            communication += b'm\n'
            
    else:
        raise ValueError("Model '{}' not available!".format(model))
    
    communication += b'y\n'
    
    proc = subprocess.Popen(args, stdin=subprocess.PIPE)
    
    for _ in range(1):
        proc.communicate(communication)
        

if __name__ == '__main__':
    
    calc_distances('testfile.alignment')
    
    D, labels = parse_distance_matrix('testfile.alignment.dist')
    
    print(D)
    print(labels)
    
    sorted_labels = sorted(labels)
    D = rearrange_matrix(D, labels, sorted_labels)
    
    print(D)
    print(sorted_labels)
    
    