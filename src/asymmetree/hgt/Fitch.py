# -*- coding: utf-8 -*-


import itertools

import networkx as nx

from asymmetree.datastructures.Tree import LCA
from asymmetree.tools.GraphTools import independent_sets


__author__ = 'David Schaller'


def true_transfer_edges(T):
    """Returns a set containing v if (u, v) is labeled as a transfer edge."""
    
    return {v for _, v in T.edges() if v.transferred}


def rs_transfer_edges(T, S, lca_S=None):
    """Transfer edges in T according to the relaxed scenario definition.
    
    An edge (u,v) in T is an (rs-)transfer edge if u and v are mapped to
    incomparable nodes/edges in the species tree S.
    """
    
    if not isinstance(lca_S, LCA):
        lca_S = LCA(S)
        
    transfer_edges = set()
    
    for u, v in T.edges():
        if not lca_S.are_comparable(u.color, v.color):
            transfer_edges.add(v)
    
    return transfer_edges


def fitch(tree, transfer_edges, supply_undirected=False, lca_T=None):
    """Returns the (directed) Fitch graph.
    
    Keyword arguments:
        supply_undirected - additionally return the undirected Fitch graph,
            default is False
        lca_T - instance of LCA corresponding to the tree, default is False
            in which case a new instance is created and used
    """
    
    if not isinstance(lca_T, LCA):
        lca_T = LCA(tree)
        
    if not isinstance(transfer_edges, (set, dict)):
        transfer_edges = set(transfer_edges)
    
    leaves = tree.supply_leaves()
    fitch = nx.DiGraph()
    
    # store for each leaf the first transfer edge on the way to the root
    first_transfer = {}
    
    for x in leaves:
        fitch.add_node(x.ID, label=x.label, color=x.color)
        
        current = x
        while current:
            if current in transfer_edges:
                first_transfer[x] = current
                break
            current = current.parent
    
    for x, y in itertools.permutations(leaves, 2):
        
        if (y in first_transfer and
            lca_T.ancestor_not_equal(lca_T(x, y), first_transfer[y])):
            fitch.add_edge(x.ID, y.ID)
    
    if not supply_undirected:
        return fitch
    else:
        return fitch, fitch.to_undirected()
    

def undirected_fitch(tree, transfer_edges, lca_T=None):
    """Returns the undirected Fitch graph.
    
    Keyword arguments:
        lca_T - instance of LCA corresponding to the tree, default is False
            in which case a new instance is created and used
    """
    
    return fitch(tree, transfer_edges, supply_undirected=True, lca_T=lca_T)[1]


def _rs_fitch_aux_graph(G, color_set, ind_sets, ignore_set):
    
    aux_graph = nx.Graph()
    aux_graph.add_nodes_from(color_set)
    
    for i, ind_set in enumerate(ind_sets):
        
        # skip the independent set that shall not be included
        if ignore_set == i:
            continue
        
        for x, y in itertools.combinations(ind_set, 2):
            X = G.nodes[x]['color']
            Y = G.nodes[y]['color']
            
            if X != Y:
                aux_graph.add_edge(X, Y)
    
    return aux_graph


def is_rs_fitch(G, color_set=None):
    """Checks whether a given graph is an rs-Fitch graph.
    
    I.e. whether the graph is the Fitch graph of some relaxed scenario (rs).
    """
    
    if color_set is None:
        color_set = set()
        for v in G.nodes():
            color_set.add(G.nodes[v]['color'])
    
    ind_sets = independent_sets(G)
    
    # not a complete multipartite graph
    if ind_sets is False:
        return False
    
    # trivial number of independent sets
    k = len(ind_sets)
    if k <= 1:
        return True
    
    # more than one independent set
    for ignore_set in range(-1, k):
        
        # ignore_set = -1 means that no set is ignored
        aux_graph = _rs_fitch_aux_graph(G, color_set, ind_sets, ignore_set)
        if not nx.is_connected(aux_graph):
            return True
    
    # if nothing was True
    return False


if __name__ == '__main__':
    
    import asymmetree.treeevolve as te
    import asymmetree.tools.GraphTools as gt
    
    S = te.simulate_species_tree(10)
    TGT = te.simulate_dated_gene_tree(S, dupl_rate=1.0, loss_rate=0.5,
                                      hgt_rate=0.5)
    OGT = te.observable_tree(TGT)
    
    print('--- S ---\n', S.to_newick())
    print('--- OGT ---\n', OGT.to_newick())
    
    transf1 = true_transfer_edges(OGT)
    transf2 = rs_transfer_edges(OGT, S)
    
    print(transf1)
    print(transf2)
    print(transf1.issuperset(transf2))
    print(transf1-transf2)
    
    fitch_d, fitch_u = fitch(OGT, transf2, supply_undirected=True)
    n = fitch_d.order()
    
    print(gt.independent_sets(fitch_u))
    
    # print(fitch_d.edges())
    # print(fitch_d.size())
    # print(fitch_u.edges())
    # print(fitch_u.size())
    
    # from asymmetree.visualize.GeneTreeVis import GeneTreeVis
    # GeneTreeVis(TGT)
    # GeneTreeVis(OGT)
    
    print('is rs-Fitch graph:',
          is_rs_fitch(fitch_u, color_set=[x.ID for x in S.leaves()]))
    
    G = nx.Graph()
    G.add_node('a', color=1)
    G.add_node('a2', color=1)
    G.add_node('b', color=2)
    G.add_node('b2', color=2)
    G.add_edges_from([('a', 'a2'), ('a', 'b2'), ('b', 'a2'), ('b', 'b2')])
    print('is rs-Fitch graph:', is_rs_fitch(G, color_set=[1,2]))
    