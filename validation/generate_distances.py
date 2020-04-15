# -*- coding: utf-8 -*-

import glob, pickle

import numpy as np

import asymmetree.tools.DistanceCalculation as dc

import treepuzzle


def parse_phylip_alignment(filename):
    
    with open(filename, 'r') as f:
        
        header = f.readline()
        header = header.split()
        n = int(header[0])
        length = int(header[1])
        
        labels = []
        sequences = {}
        
        for i in range(n):
            line = f.readline()
            line = line.split()
            
            label = line[0]
            seq_part = "".join(line[1:])
            
            labels.append(label)
            sequences[label] = seq_part
            
        while len(sequences[labels[0]]) < length:
            f.readline()
            
            for i in range(n):
                line = f.readline()
                sequences[labels[i]] += "".join(line.split())
                
    return labels, sequences


def JC69_distances(labels, seqs):
    
    N = len(labels)
    D = np.zeros((N,N))
    
    for i in range(N-1):
        for j in range(i, N):
            p, n = dc.p_distance(seqs[labels[i]], seqs[labels[j]])
            d = dc.JC69_distance(p)
            
            D[i, j] = d
            D[j, i] = d
    
    return D


def K80_distances(labels, seqs):
    
    N = len(labels)
    D = np.zeros((N,N))
    
    for i in range(N-1):
        for j in range(i, N):
            
            S, V, n = dc.IV_proportions(seqs[labels[i]], seqs[labels[j]])
            d, kappa = dc.K80_distance(S, V)
            
            D[i, j] = d
            D[j, i] = d
    
    return D
            

if __name__ == '__main__':
    
    directory = 'testfiles_JC69_seqs'
    files = glob.glob(directory + '/*.phylip')
    matrices = []
    
    for i in range(len(files)):
        
        labels, seqs = parse_phylip_alignment("{}/{}.phylip".format(directory, i))
        
        D = JC69_distances(labels, seqs)
        matrices.append(D)
        
    with open(directory + '/distances.pickle', 'wb') as f:
        
        pickle.dump(matrices, f)
        
        
        
    directory = 'testfiles_K80_seqs'
    files = glob.glob(directory + '/*.phylip')
    matrices = []
    
    for i in range(len(files)):
        
        labels, seqs = parse_phylip_alignment("{}/{}.phylip".format(directory, i))
        
        D = K80_distances(labels, seqs)
        matrices.append(D)
        
    with open(directory + '/distances.pickle', 'wb') as f:
        
        pickle.dump(matrices, f)

