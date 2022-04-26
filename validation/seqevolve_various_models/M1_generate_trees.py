# -*- coding: utf-8 -*-

import os, pickle, glob

from tralda.datastructures.Tree import Tree

import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import distance_matrix


def simulate(directory, number_of_trees, species_per_tree):
    
    if not os.path.exists(directory):
        os.mkdir(directory)  
    
    for i in range(number_of_trees):
        
        S = te.simulate_species_tree(50)
        T_simulator = te.GeneTreeSimulator(S)
        T = T_simulator.simulate()   # dupl./loss/HGT disabled
        
        te.rate_heterogeneity(T, S, autocorr_variance=0.2)
        
        T.serialize('{}/scenario{}.pickle'.format(directory, i))
    
       
def load(directory):
    
    files = glob.glob(directory + '/scenario*.pickle')
    trees = []
    
    for i in range(len(files)):
        trees.append(Tree.load('{}/scenario{}.pickle'.format(directory, i)))
            
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
    
    directory = 'testfiles_scenarios'
    number_of_trees = 100
    species_per_tree = 50
    
    simulate(directory, number_of_trees, species_per_tree)
    
    true_distances(directory)