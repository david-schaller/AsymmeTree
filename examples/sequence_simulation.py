# -*- coding: utf-8 -*-

import asymmetree.seqevolve as se
import asymmetree.treeevolve as te


__author__ = 'David Schaller'

# specifify models
subst_model = se.SubstModel('a', 'CUSTOM',
                            filename='../resources/subst_matrices/WAG.paml')
indel_model = se.IndelModel(0.01, 0.01, length_model='zipf')

# initialize evolver
evolver = se.Evolver(subst_model, indel_model=indel_model, jump_chain=False)
print(evolver.subst_model.Q)

# simulate along a tree
T = te.simulate_species_tree(5)
evolver.evolve_along_tree(T, start_length=150)

for node, sequence in evolver.sequences.items():
    print(node.label, subst_model.to_sequence(sequence))
    
alg_seq = evolver.true_alignment(write_to='testfile.alignment')
for node, sequence in alg_seq.items():
    print(node.label, sequence)