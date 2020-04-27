# -*- coding: utf-8 -*-

import glob, pickle

import numpy as np

import asymmetree.tools.DistanceCalculation as dc

from validation import aux_functions


def JC69_distances(labels, seqs):
    
    N = len(labels)
    D = np.zeros((N,N))
    
    for i in range(N-1):
        for j in range(i, N):
            
            d = dc.JC69_distance(seqs[labels[i]], seqs[labels[j]])
            
            D[i, j] = d
            D[j, i] = d
    
    return D


def K80_distances(labels, seqs):
    
    N = len(labels)
    D = np.zeros((N,N))
    
    for i in range(N-1):
        for j in range(i, N):
            
            d, kappa = dc.K80_distance(seqs[labels[i]], seqs[labels[j]])
            
            D[i, j] = d
            D[j, i] = d
    
    return D


def all_nuc_distances(directory, model):
    
    files = glob.glob(directory + '/*.phylip')
    matrices = []
    
    for i in range(len(files)):
        
        print("Model '{}', File {} of {}".format(model, i+1, len(files)))
        
        labels, seqs = aux_functions.parse_phylip_alignment("{}/{}.phylip".format(directory, i))
        
        if model == 'JC69':
            D = JC69_distances(labels, seqs)
        elif model == 'K80':
            D = K80_distances(labels, seqs)
            
        matrices.append((labels, D))
        
    with open(directory + '/distances.pickle', 'wb') as f:
        
        pickle.dump(matrices, f)
        

def puzzle_distances(directory, model):
    
    files = glob.glob(directory + '/*.phylip')
    matrices = []
    
    for i in range(len(files)):
        
        print("Model '{}', File {} of {}".format(model, i+1, len(files)))
        
        aux_functions.calc_distances("{}/{}.phylip".format(directory, i),
                                  model=model)
        
        D, labels = aux_functions.parse_distance_matrix("{}/{}.phylip.dist".format(directory, i))
        
        matrices.append((labels, D))
        
    with open(directory + '/distances.pickle', 'wb') as f:
        
        pickle.dump(matrices, f)
        
        
def summarize_distances(models, directories, outfile):
    
    distances = []
    for directory in directories:
        
        with open(directory + '/distances.pickle', 'rb') as f:
            distances.append( pickle.load(f))
            
    true_distances = distances[0]
    
    table = []
    
    for i in range(len(true_distances)):
        
        labels, true_D = true_distances[i]
        
        Ds = []
        for estimated_distances in distances[1:]:
            D_labels, D_unsorted = estimated_distances[i]
            D = aux_functions.rearrange_matrix(D_unsorted, D_labels, labels)
            Ds.append(D)
        
        for k in range(true_D.shape[0]-1):
            for l in range(k+1, true_D.shape[0]):
                
                row = [true_D[k, l]]
                
                for D in Ds:
                    row.append( D[k, l] )
                    
                table.append(row)
                
    with open(outfile, 'w') as f:
        
        f.write(','.join(models))
        
        for row in table:
            f.write('\n' + ','.join(str(value) for value in row))
    
            

if __name__ == '__main__':
    
    models = ['true', 'JC69', 'K80', 'WAG', 'JTT']
    directories = ['testfiles_scenarios'] + ['testfiles_{}_seqs'.format(m) for m in models[1:]]
    outfile = 'results/all_distances.csv'
    
    all_nuc_distances(directories[1], models[1])
    all_nuc_distances(directories[2], models[2])
    
    puzzle_distances(directories[3], models[3])
    puzzle_distances(directories[4], models[4])
    
    summarize_distances(models, directories, outfile)