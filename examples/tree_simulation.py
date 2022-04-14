# -*- coding: utf-8 -*-

import asymmetree.treeevolve as te
from asymmetree.analysis.BestMatches import lrt_from_tree
from asymmetree.tools.PhyloTreeTools import (to_newick,)

D = 1.0
L = 1.0
H = 0.0


# --------------------------------------------------------------------------
#                            SPECIES TREE
# --------------------------------------------------------------------------

S = te.simulate_species_tree(10, planted=True, non_binary_prob=0.2)
print('------------- S -------------')
print(to_newick(S))

# --------------------------------------------------------------------------
#                             GENE TREE
# --------------------------------------------------------------------------

TGT_simulator = te.GeneTreeSimulator(S)
TGT = TGT_simulator.simulate(dupl_rate=D, loss_rate=L, hgt_rate=H,
                             prohibit_extinction='per_species')

TGT = te.assign_rates(TGT, S, base_rate=1,
                      autocorr_variance=0.2,
                      rate_increase=('gamma', 0.5, 2.2),
                      CSN_weights=(1, 1, 1))
print('------------- TGT -------------')
print(to_newick(TGT))
print('all species have at least one copy:', TGT_simulator._assert_no_extinction(TGT))

# --------------------------------------------------------------------------
#                       OBSERVABLE GENE TREE
# --------------------------------------------------------------------------

PGT = te.prune_losses(TGT)
print('------------- PGT -------------')
print(to_newick(PGT))

# --------------------------------------------------------------------------
#                       LEAST RESOLVED TREE
# --------------------------------------------------------------------------

lrt = lrt_from_tree(PGT)
print('------------- LRT -------------')
print(to_newick(lrt))