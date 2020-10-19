# -*- coding: utf-8 -*-

import itertools
import networkx as nx

from asymmetree.tools.TreeTools import LCA


__author__ = 'David Schaller'


def below_equal_above(T, S):
    
    L_T = [l for l in T.leaves()]
    L_T_ids = [l.ID for l in L_T]
    L_S = {l.ID: l for l in S.leaves()}
    
    below = nx.Graph()
    below.add_nodes_from(L_T_ids)
    
    equal = nx.Graph()
    equal.add_nodes_from(L_T_ids)
    
    above = nx.Graph()
    above.add_nodes_from(L_T_ids)
    
    lca_T = LCA(T)
    lca_S = LCA(S)
    
    for a, b in itertools.combinations(L_T, 2):
        
        t_ab = lca_T.get(a, b).tstamp
        t_AB = lca_S.get(L_S[a.color], L_S[b.color]).tstamp
        
        if t_ab < t_AB:
            below.add_edge(a.ID, b.ID)
        elif t_ab == t_AB:
            equal.add_edge(a.ID, b.ID)
        else:
            above.add_edge(a.ID, b.ID)
    
    
    return below, above, equal


if __name__ == '__main__':
    
    import asymmetree.treeevolve as te
    
    S = te.simulate_species_tree(10)
    TGT = te.simulate_dated_gene_tree(S, dupl_rate=1.0, loss_rate=0.5,
                                      hgt_rate=0.5)
    OGT = te.observable_tree(TGT)
    
    print('--- S ---\n', S.to_newick())
    print('--- OGT ---\n', OGT.to_newick())
    
    below, above, equal = below_equal_above(OGT, S)
    n = below.order()
    print('Genes:', n, 'Total relations:', int(n * (n-1) / 2))
    print('< {}\n= {}\n> {}'.format(below.size(), equal.size(), above.size()))