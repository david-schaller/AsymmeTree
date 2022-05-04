# -*- coding: utf-8 -*-


import asymmetree.treeevolve as te
from asymmetree.visualize.TreeVis import visualize, assign_colors


# --------------------------------------------------------------------------
#                            SPECIES TREE
# --------------------------------------------------------------------------

S = te.simulate_species_tree(6, planted=True, contraction_probability=0.2)
# S = te.simulate_species_tree(6, model='yule', birth_rate=1.0, 
#                              rescale_to_height=1.0)


# --------------------------------------------------------------------------
#                             GENE TREE
# --------------------------------------------------------------------------

T_simulator = te.GeneTreeSimulator(S)
T = T_simulator.simulate(dupl_rate=1.0,
                         loss_rate=1.0,
                         hgt_rate=1.0,
                         dupl_polytomy=0.5,
                         replace_prob=0.5,
                         transfer_distance_bias='inverse')


# --------------------------------------------------------------------------
#                         RATE HETEROGENEITY
# --------------------------------------------------------------------------

te.rate_heterogeneity(T, S, inplace=True)


# --------------------------------------------------------------------------
#                           VISUALIZATION
# --------------------------------------------------------------------------

species_colors, gene_colors = assign_colors(S, T)

visualize(S, color_dict=species_colors, save_as='testfile_speciestree.pdf')
visualize(T, color_dict=gene_colors, save_as='testfile_genetree.pdf')
