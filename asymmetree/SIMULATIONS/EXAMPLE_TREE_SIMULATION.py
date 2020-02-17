# -*- coding: utf-8 -*-

import asymmetree.simulator.TreeSimulator as ts
import asymmetree.simulator.TreeImbalancer as tm
from asymmetree.simulator.Scenario import Scenario
from asymmetree.best_matches.Quartets import Quartets
import asymmetree.best_matches.LRTConstructor as lrt

# --------------------------------------------------------------------------
#                            RATES FOR
#      duplication, loss and horizontal gene transfer events
# --------------------------------------------------------------------------

DLH_rates = (1,1,0)

# --------------------------------------------------------------------------
#                            SPECIES TREE
# --------------------------------------------------------------------------

S = ts.build_species_tree(10, planted=True, non_binary=0.2)
print("------------- S -------------")
print(S.to_newick())

# --------------------------------------------------------------------------
#                             GENE TREE
# --------------------------------------------------------------------------

TGT = ts.build_gene_tree(S, DLH_rates)
TGT = tm.imbalance_tree(TGT, S, baseline_rate=1,
                        lognormal_v=0.2,
                        gamma_param=(0.5, 1.0, 2.2),
                        weights=(1, 1, 1),
                        copy_tree=False)
print("------------- TGT -------------")
print(TGT.to_newick(distance_only=False))

# --------------------------------------------------------------------------
#                       OBSERVABLE GENE TREE
# --------------------------------------------------------------------------

OGT = ts.observable_tree(TGT)
print("------------- OGT -------------")
print(OGT.to_newick(distance_only=False))

# --------------------------------------------------------------------------
#                       LEAST RESOLVED TREE
# --------------------------------------------------------------------------

LRT = lrt.LRT_from_observable_tree(OGT)
print("------------- LRT -------------")
print(LRT.to_newick(distance_only=False))

# --------------------------------------------------------------------------
#                              SCENARIO
# --------------------------------------------------------------------------

scenario = Scenario(S, TGT, DLH_rates, OGT=OGT)         # wrap everything in a scenario
D = scenario.get_distance_matrix()                      # compute the distance matrix


# --------------------------------------------------------------------------
#                 QUARTET METHOD FOR BMG INFERENCE
# --------------------------------------------------------------------------


qu = Quartets(scenario, D, voting_mode="majority",      # inititialize 'Quartets' instance
              closest_outgroups=True)
qu.build_graphs()                                       # build BMG

print("No. of edges in the true graph:",
      scenario.BMG_subtrees.size())

print("No. of edges in inferred graph:",
      qu.BMG.size())

for u, v in scenario.BMG_subtrees.edges:                # print missing edges
    if not qu.BMG.has_edge(u,v):
        print("Edge ({},{}) is missing!".format(u,v))