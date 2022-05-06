# -*- coding: utf-8 -*-

import asymmetree.seqevolve as se
import asymmetree.treeevolve as te


__author__ = 'David Schaller'

# specify models
subst_model = se.SubstModel('a', 'CUSTOM',
                            filename='../resources/subst_matrices/WAG.paml')
indel_model = se.IndelModel(0.01, 0.01, length_distr=('zipf', 1.821))
#indel_model = se.IndelModel(0.01, 0.01, length_distr=('negative_binomial', 1, 0.5))

# initialize evolver
evolver = se.Evolver(subst_model, indel_model=indel_model, gillespie=False)
print(evolver.subst_model.Q)

# simulate along a tree
T = te.species_tree_N_age(5, 1.0)
evolver.evolve_along_tree(T, start_length=150)

evolver.write_sequences('testfile.fasta', include_inner=True)
for node, sequence in evolver.sequences.items():
    print(node.label, subst_model.to_sequence(sequence))

alg_seq = evolver.true_alignment(write_to='testfile.alignment')
for node, sequence in alg_seq.items():
    print(node.label, sequence)
