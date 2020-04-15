# -*- coding: utf-8 -*-

import os

from asymmetree.seqevolve import Evolver, SubstModel

from generate_trees import load

subst_models = [SubstModel('n', 'JC69'),
                SubstModel('n', 'K80', params={'kappa': 0.6}),
                SubstModel('a', 'WAG'),
                SubstModel('a', 'JTT'),]

directory = 'testfiles_{}_seqs'

trees = load()

for subst_model in subst_models:
    
    dir_name = directory.format(subst_model.model_name)
    if not os.path.exists(dir_name):
        os.mkdir(dir_name) 
        
    evolver = Evolver(subst_model, jump_chain=False)
    
    for i in range(len(trees)):
        
        evolver.evolve_along_tree(trees[i], start_length=500)
        
        evolver.true_alignment(include_inner=False,
                               write_to="{}/{}.phylip".format(dir_name, i))
        
        

    
