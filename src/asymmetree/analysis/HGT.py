# -*- coding: utf-8 -*-


import itertools

import networkx as nx

from tralda.datastructures import Tree, TreeNode, LCA
from tralda.tools.GraphTools import independent_sets
from tralda.cograph.Cograph import (to_cotree, paths_of_length_2)
from tralda.supertree import Build

from asymmetree.tools.PhyloTreeTools import (add_planted_root,
                                             phylo_tree_attributes,
                                             assign_missing_labels)
from asymmetree.treeevolve.SpeciesTree import (simulate_timing,
                                               distance_from_timing)


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                     Transfer edges and Fitch graphs
# --------------------------------------------------------------------------

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
    
    leaves = [l for l in tree.leaves()]
    fitch = nx.DiGraph()
    
    # store for each leaf the first transfer edge on the way to the root
    first_transfer = {}
    
    for x in leaves:
        fitch.add_node(x.label, color=x.color)
        
        current = x
        while current:
            if current in transfer_edges:
                first_transfer[x] = current
                break
            current = current.parent
    
    for x, y in itertools.permutations(leaves, 2):
        
        if (y in first_transfer and
            lca_T.ancestor_not_equal(lca_T(x, y), first_transfer[y])):
            fitch.add_edge(x.label, y.label)
    
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


# --------------------------------------------------------------------------
#                   Later-divergence-time (LDT) graphs
# --------------------------------------------------------------------------

def below_equal_above(T, S, lca_T=None, lca_S=None):
    
    L_T = [l for l in T.leaves()]
    L_S = {l.label: l for l in S.leaves()}
    
    below = nx.Graph()
    equal = nx.Graph()
    above = nx.Graph()
    for u in L_T:
        below.add_node(u.label, color=u.color)
        equal.add_node(u.label, color=u.color)
        above.add_node(u.label, color=u.color)
    
    if not isinstance(lca_T, LCA):
        lca_T = LCA(T)
    if not isinstance(lca_S, LCA):
        lca_S = LCA(S)
    
    for a, b in itertools.combinations(L_T, 2):
        
        t_ab = lca_T.get(a, b).tstamp
        t_AB = lca_S.get(L_S[a.color], L_S[b.color]).tstamp
        
        if t_ab < t_AB:
            below.add_edge(a.label, b.label)
        elif t_ab == t_AB:
            equal.add_edge(a.label, b.label)
        else:
            above.add_edge(a.label, b.label)
    
    
    return below, above, equal


def ldt_graph(T, S, lca_T=None, lca_S=None):
    """Later-divergence-time graph.
    
    Returns a graph with the leaves of the gene tree T as vertex set, and 
    edges ab if and only if a and b diverged later than the corresponding
    species A and B in the species tree S."""
    
    ldt, _, _ = below_equal_above(T, S, lca_T=lca_T, lca_S=lca_S)
    return ldt


class RsScenarioConstructor:
    
    
    def __init__(self, colored_cograph, color_set=None):
        
        self.G = colored_cograph
        
        self.L = [x for x in self.G.nodes()]
        
        if color_set is None:
            color_set = set()
            for v in self.L:
                color_set.add(self.G.nodes[v]['color'])
        elif not isinstance(color_set, (set, dict)):
            color_set = set(color_set)
        self.color_set = color_set
        
        
    def run(self):
        
        self.S = self._species_tree()
        if not self.S:
            return False
        
        self.epsilon = float('inf')
        for u, v in self.S.edges():
            self.epsilon = min(self.epsilon,
                               abs(u.tstamp - v.tstamp))
        self.epsilon /= 3.0
        
        self.S_leaves = self.S.leaf_dict()
        
        # top-level call on full leaf set and rho_S
        self.T = Tree(self._build_gene_tree(self.L, self.S.root.children[0]))
        
        # --- post-processing of T ---
        
        # planted root
        add_planted_root(self.T)
        self.T.root.tstamp = self.S.root.tstamp
        self.T.root.color = self.S.root.label
        
        # suppress all vertices with a single child except the planted root
        to_suppress = []
        for u in self.T.preorder():
            if u.parent and len(u.children) == 1:
                to_suppress.append((u.parent, u))
        self.T.contract(to_suppress)
        
        distance_from_timing(self.T)
        
        return self.S, self.T
        
    
    def _species_tree(self):
        
        cotree = to_cotree(self.G)
        
        if not cotree:
            return False
        
        triples = set()
        
        for a, b, c in paths_of_length_2(cotree):
            
            A = self.G.nodes[a.label]['color']
            B = self.G.nodes[b.label]['color']
            C = self.G.nodes[c.label]['color']
            
            if A != B and A != C and B != C:
                
                # sort to avoid redundancy in triple set
                if A > C:
                    A, C = C, A
                
                triples.add( (A, C, B) )
        
        build = Build(triples, self.color_set, mincut=False)
        S = build.build_tree()
        
        if not S:
            return False
        else:
            phylo_tree_attributes(S, inplace=True)
            add_planted_root(S)
            assign_missing_labels(S)
            simulate_timing(S)
            distance_from_timing(S)
            return S
    
    
    def _build_gene_tree(self, L, u_S):
        
        u_T = TreeNode(label='', event='D',
                       color=(u_S.parent.label, u_S.label),
                       tstamp=u_S.tstamp + self.epsilon)
        
        # u_S is a leaf
        if not u_S.children:
            for x in L:
                leaf = TreeNode(label=x,
                                color=self.G.nodes[x]['color'],
                                tstamp=0.0)
                u_T.add_child(leaf)
        
        # u_S is an inner vertex
        else:
            
            # maps color to the respective child of u_S
            color_to_v_S = {}
            for v_S in u_S.children:
                for x in self.S_leaves[v_S]:
                    color_to_v_S[x.label] = v_S
                    
            # connected components C of G[L']
            for C in nx.connected_components(self.G.subgraph(L)):
                
                # choose v_S_star s.t. L(S(v_S_star)) \cap sigma(C) is non-empty
                v_S_star = color_to_v_S[self.G.nodes[next(iter(C))]['color']]
                v_T = TreeNode(label='', event='H',
                               color=(u_S.label, v_S_star.label),
                               tstamp=u_S.tstamp - self.epsilon)
                u_T.add_child(v_T)
                
                # aux. graph for equivalence classes
                aux_graph = nx.Graph()
                aux_graph.add_nodes_from(C)
                
                for x, y in itertools.combinations(C, 2):
                    if (color_to_v_S[self.G.nodes[x]['color']] is
                        color_to_v_S[self.G.nodes[y]['color']]):
                        aux_graph.add_edge(x, y)
                
                # for each equivalence class K that is a subset of C
                for K in nx.connected_components(aux_graph):
                    
                    # choose v_S s.t. \sigma(K) \subseteq L(S(v_S))
                    v_S = color_to_v_S[self.G.nodes[next(iter(K))]['color']]
                    
                    w_K = self._build_gene_tree(K, v_S)
                    v_T.add_child(w_K)
                    
                    if v_S is not v_S_star:
                        w_K.transferred = 1
        
        return u_T
