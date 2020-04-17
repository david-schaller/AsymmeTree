# -*- coding: utf-8 -*-

import os

from asymmetree.seqevolve import Evolver, SubstModel

from M1_generate_trees import load

seq_length = 500

subst_models = [SubstModel('n', 'JC69'),
                SubstModel('n', 'K80', kappa=2.0),
                SubstModel('a', 'WAG'),
                SubstModel('a', 'JTT'),]

directory = 'testfiles_{}_seqs'

trees = load('testfiles_scenarios')

for subst_model in subst_models:
    
    dir_name = directory.format(subst_model.model_name)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name) 
        
    evolver = Evolver(subst_model, jump_chain=False)
    
    for i in range(len(trees)):
        
        evolver.evolve_along_tree(trees[i], start_length=seq_length)
        
        evolver.true_alignment(include_inner=False,
                               write_to="{}/{}.phylip".format(dir_name, i))
        
        

    
