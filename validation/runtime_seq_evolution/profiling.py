# -*- coding: utf-8 -*-

import cProfile

import asymmetree.treeevolve as te
import asymmetree.seqevolve as se


tree = te.simulate_species_tree(50, model='innovation',
                                planted=False)

subst_model = se.SubstModel('a', 'WAG')

# het_model = se.HetModel(1.0, classes=5)
het_model = None

evolver = se.Evolver(subst_model, het_model=het_model)

cProfile.run('evolver.evolve_along_tree(tree, start_length=10000)')
