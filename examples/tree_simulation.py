# -*- coding: utf-8 -*-

import asymmetree.treeevolve as te
from asymmetree.analysis.BestMatches import lrt_from_tree
from asymmetree.tools.PhyloTreeTools import (to_newick,)

# --------------------------------------------------------------------------
#                            SPECIES TREE
# --------------------------------------------------------------------------

S = te.simulate_species_tree(10, planted=True, contraction_probability=0.0,
                             contraction_proportion=0.2,
                             contraction_bias='exponential')
print('------------- S -------------')
print(to_newick(S))

# --------------------------------------------------------------------------
#                             GENE TREE
# --------------------------------------------------------------------------

TGT_simulator = te.GeneTreeSimulator(S)
TGT = TGT_simulator.simulate(dupl_rate=1.0, loss_rate=1.0, hgt_rate=0.2,
                             prohibit_extinction='per_species')

# --------------------------------------------------------------------------
#                         RATE HETEROGENEITY
# --------------------------------------------------------------------------

TGT = te.rate_heterogeneity(TGT, S, base_rate=1.0,
                            autocorr_variance=0.2,
                            rate_increase=('gamma', 0.5, 2.2),
                            CSN_weights=(1, 1, 1))
print('------------- TGT -------------')
print(to_newick(TGT))
print('all species have at least one copy:', 
      TGT_simulator._assert_no_extinction(TGT))

# --------------------------------------------------------------------------
#                          PRUNED GENE TREE
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
