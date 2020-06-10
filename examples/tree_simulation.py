# -*- coding: utf-8 -*-

import asymmetree.treeevolve as te
from asymmetree.best_matches import lrt_from_observable_tree, Quartets


DLH_rates=(1.0, 1.0, 0.0)

# --------------------------------------------------------------------------
#                            SPECIES TREE
# --------------------------------------------------------------------------

S = te.simulate_species_tree(10, planted=True, non_binary_prob=0.2)
print('------------- S -------------')
print(S.to_newick())

# --------------------------------------------------------------------------
#                             GENE TREE
# --------------------------------------------------------------------------

TGT_simulator = te.GeneTreeSimulator(S)
TGT = TGT_simulator.simulate(DLH_rates=DLH_rates,
                             prohibit_extinction='per_species')

TGT = te.assign_rates(TGT, S, base_rate=1,
                      autocorr_variance=0.2,
                      gamma_param=(0.5, 1.0, 2.2),
                      CSN_weights=(1, 1, 1))
print('------------- TGT -------------')
print(TGT.to_newick())
print('all species have at least one copy:', TGT_simulator._assert_no_extinction(TGT))

# --------------------------------------------------------------------------
#                       OBSERVABLE GENE TREE
# --------------------------------------------------------------------------

OGT = te.observable_tree(TGT)
print('------------- OGT -------------')
print(OGT.to_newick())

# --------------------------------------------------------------------------
#                       LEAST RESOLVED TREE
# --------------------------------------------------------------------------

lrt = lrt_from_observable_tree(OGT)
print('------------- LRT -------------')
print(lrt.to_newick())

# --------------------------------------------------------------------------
#                              SCENARIO
# --------------------------------------------------------------------------

scenario = te.Scenario(S, TGT, DLH_rates, OGT=OGT)      # wrap everything in a scenario
D = scenario.get_distance_matrix()                      # compute the distance matrix


# --------------------------------------------------------------------------
#                 QUARTET METHOD FOR BMG INFERENCE
# --------------------------------------------------------------------------


qu = Quartets(scenario, D, voting_mode='majority',      # inititialize 'Quartets' instance
              closest_outgroups=True)
qu.build_graphs()                                       # build BMG

print('no. of edges in the true graph:',
      scenario.bmg_subtrees.size())

print('no. of edges in inferred graph:',
      qu.bmg.size())

for e in scenario.bmg_subtrees.edges:                # print missing edges
    if not qu.bmg.has_edge(*e):
        print('edge {} is missing'.format(e))