from __future__ import annotations

from tralda.utils.graph_tools import contingency_table
from tralda.utils.graph_tools import graphs_equal

import asymmetree.treeevolve as te
from asymmetree.analysis import (
    undirected_fitch,
    rs_transfer_edges,
    below_equal_above,
    ldt_graph,
    RsScenarioConstructor,
)
from asymmetree.utils.phylogenetic_trees import (
    to_newick,
)

S = te.species_tree_n_age(10, 1.0)
TGT = te.dated_gene_tree(S, dupl_rate=1.0, loss_rate=0.5, hgt_rate=0.5)
PGT = te.prune_losses(TGT)

print("--- S ---\n", to_newick(S))
print(to_newick(S, distance=False, label_inner=False))
print("--- PGT ---\n", to_newick(PGT))

ldt, above, equal = below_equal_above(PGT, S)
fitch = undirected_fitch(PGT, rs_transfer_edges(PGT, S))
n = ldt.order()
print("Genes:", n, "Total relations:", int(n * (n - 1) / 2))
print("< {}\n= {}\n> {}".format(ldt.size(), equal.size(), above.size()))

rs_scen_constr = RsScenarioConstructor(ldt)
result = rs_scen_constr.run()

if result:
    S2, T2 = result
    print("--- S2 ---\n", to_newick(S2, distance=False))
    print("--- T2 ---\n", to_newick(T2, distance=False))
    ldt2 = ldt_graph(T2, S2)
    print(ldt2.order(), ldt2.size(), graphs_equal(ldt, ldt2))

    print("--- fitch ---")
    fitch2 = undirected_fitch(T2, rs_transfer_edges(T2, S2))
    print("Order: {} vs {}".format(fitch.order(), fitch2.order()))
    print("Size: {} vs {}".format(fitch.size(), fitch2.size()))
    print(contingency_table(fitch, fitch2))
else:
    print(False)
