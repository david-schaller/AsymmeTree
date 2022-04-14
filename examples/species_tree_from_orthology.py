# -*- coding: utf-8 -*-

import asymmetree.treeevolve as te
from asymmetree.analysis.BestMatches import orthology_from_tree, bmg_from_tree
from asymmetree.paraphylo import TreeReconstructor
from asymmetree.tools.PhyloTreeTools import (to_newick, parse_newick,
                                             reconstruct_timestamps)


__author__ = 'David Schaller'


# species tree

#S = ts.simulate_species_tree(10)
S = parse_newick('(((((16:0.38786287055727103,(18:0.2071277923445058,19:0.2071277923445058)17:0.18073507821276524)12:0.10075553853805931,(14:0.13224182895383052,15:0.13224182895383052)13:0.3563765801414998)4:0.07517286794939665,(6:0.5373882998574596,(8:0.4434182448023457,(10:0.04929450217312242,11:0.04929450217312242)9:0.3941237426292233)7:0.0939700550551139)5:0.02640297718726732)2:0.2512472266526016,3:0.8150385036973286)1:0.18496149630267142)0:0.0;')
reconstruct_timestamps(S)

print('------------- original species tree -------------')
print(to_newick(S))

tr = TreeReconstructor(cotree_mode='best')


# gene families

for i in range(100):
    TGT_simulator = te.GeneTreeSimulator(S)
    TGT = TGT_simulator.simulate(dupl_rate=1.0, loss_rate=1.0)
    TGT = te.assign_rates(TGT, S, base_rate=1,
                          autocorr_variance=0.2,
                          rate_increase=('gamma', 0.5, 2.2),
                          CSN_weights=(1, 1, 1))
    PGT = te.prune_losses(TGT)
    
#     #add orthology graph
#    ortho_graph = orthology_from_tree(PGT)
#    tr.add_ortho_graph(ortho_graph)
    
    # add RBMG
    _, rbmg = bmg_from_tree(PGT, supply_rbmg=True)
    tr.add_ortho_graph(rbmg)
    

# estimation

S_estimate = tr.build_species_tree(mode='mincut')
print('------------- MINCUT species tree -------------')
print(tr.newick_with_support())

S_estimate2 = tr.build_species_tree(mode='BPMF')
print('------------- BPMF species tree -------------')
print(tr.newick_with_support())

S_estimate3 = tr.build_species_tree(mode='greedy')
print('------------- GREEDY species tree -------------')
print(tr.newick_with_support())