# -*- coding: utf-8 -*-

import os, pickle, glob

from tralda.datastructures.Tree import Tree

import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import distance_matrix


directory = 'testfiles_scenarios'
number_of_trees = 1000


def simulate():
    
    if not os.path.exists(directory):
        os.mkdir(directory)  
    
    for i in range(number_of_trees):
        
        S = te.simulate_species_tree_age(1.0, model='yule')
        S.serialize('{}/scenario{}.pickle'.format(directory, i))
        
        node_rates, _ = te.autocorrelation_factors(S, variance)
    
       
def load(directory):
    
    files = glob.glob(directory + '/scenario*.pickle')
    trees = []
    
    for i in range(len(files)):
        with open('{}/scenario{}.pickle'.format(directory, i), 'rb') as f:
            T_nx = pickle.load(f)
            T = Tree.parse_nx(*T_nx)
            trees.append(T)
            
    return trees


def true_distances(directory):
    
    trees = load(directory)
    
    matrices = []
    for T in trees:
        
        leaves, D = distance_matrix(T)
        labels = [v.label for v in leaves]
        
        matrices.append((labels, D))
        
    with open(directory + '/distances.pickle', 'wb') as f:
        
        pickle.dump(matrices, f)
    
        
if __name__ == '__main__':
    
    simulate(directory, number_of_trees)
    
    true_distances(directory)