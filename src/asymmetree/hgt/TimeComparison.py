# -*- coding: utf-8 -*-

import itertools
import networkx as nx

from asymmetree.datastructures.Tree import LCA


__author__ = 'David Schaller'


def below_equal_above(T, S, lca_T=None, lca_S=None):
    
    L_T = [l for l in T.leaves()]
    L_T_ids = [l.ID for l in L_T]
    L_S = {l.ID: l for l in S.leaves()}
    
    below = nx.Graph()
    below.add_nodes_from(L_T_ids)
    
    equal = nx.Graph()
    equal.add_nodes_from(L_T_ids)
    
    above = nx.Graph()
    above.add_nodes_from(L_T_ids)
    
    if not isinstance(lca_T, LCA):
        lca_T = LCA(T)
    if not isinstance(lca_S, LCA):
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


def ldt_graph(T, S, lca_T=None, lca_S=None):
    """Later-divergence-time graph.
    
    Returns a graph with the leaves of the gene tree T as vertex set, and 
    edges ab if and only if a and b diverged later than the corresponding
    species A and B in the species tree S."""
    
    ldt, _, _ = below_equal_above(T, S, lca_T=lca_T, lca_S=lca_S)
    return ldt


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