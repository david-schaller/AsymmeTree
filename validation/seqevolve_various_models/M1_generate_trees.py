# -*- coding: utf-8 -*-

import os, pickle, glob

from asymmetree.datastructures.PhyloTree import PhyloTree
import asymmetree.treeevolve as te


def simulate(directory, number_of_trees, species_per_tree):
    
    if not os.path.exists(directory):
        os.mkdir(directory)  
    
    for i in range(number_of_trees):
        
        S = te.simulate_species_tree(50)
        T_simulator = te.GeneTreeSimulator(S)
        T = T_simulator.simulate((0.0, 0.0, 0.0))   # dupl./loss/HGT disabled
        
        te.assign_rates(T, S, autocorr_variance=0.2)
        
        T_nx = T.to_nx()
        with open('{}/scenario{}.pickle'.format(directory, i), 'wb') as f:
            pickle.dump(T_nx, f)
    
       
def load(directory):
    
    files = glob.glob(directory + '/scenario*.pickle')
    trees = []
    
    for i in range(len(files)):
        with open('{}/scenario{}.pickle'.format(directory, i), 'rb') as f:
            T_nx = pickle.load(f)
            T = PhyloTree.parse_nx(*T_nx)
            trees.append(T)
            
    return trees


def true_distances(directory):
    
    trees = load(directory)
    
    matrices = []
    for T in trees:
        
        leaves, D = T.distance_matrix()
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