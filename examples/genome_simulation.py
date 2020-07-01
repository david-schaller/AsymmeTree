# -*- coding: utf-8 -*-

from asymmetree.treeevolve import simulate_species_tree
from asymmetree.genome import GenomeSimulator
from asymmetree.seqevolve import SubstModel, IndelModel, HetModel


__author__ = 'David Schaller'


species_tree = simulate_species_tree(10, model='innovation')

subst_model = SubstModel('a', 'JTT')
indel_model = IndelModel(0.01, 0.01, length_distr=('zipf', 1.821))
het_model = None

genome_sim = GenomeSimulator(species_tree, outdir='testfile_genome')

genome_sim.simulate_gene_trees(50, dupl_rate=1.0, loss_rate=0.5,
                               base_rate=('gamma', 1.0, 1.0),
                               prohibit_extinction='per_species')

genome_sim.simulate_sequences(subst_model,
                              indel_model=indel_model,
                              het_model=het_model,
                              length_distr=('constant', 200))