# -*- coding: utf-8 -*-

import itertools
import networkx as nx

from asymmetree.datastructures.Tree import LCA
from asymmetree.datastructures.PhyloTree import PhyloTree, PhyloTreeNode
from asymmetree.cograph.Cograph import Cotree
from asymmetree.tools.Build import Build
from asymmetree.treeevolve.SpeciesTree import (simulate_timing,
                                               distance_from_timing)


__author__ = 'David Schaller'


def below_equal_above(T, S, lca_T=None, lca_S=None):
    
    L_T = [l for l in T.leaves()]
    L_S = {l.ID: l for l in S.leaves()}
    
    below = nx.Graph()
    equal = nx.Graph()
    above = nx.Graph()
    for u in L_T:
        below.add_node(u.ID, label=u.label, color=u.color)
        equal.add_node(u.ID, label=u.label, color=u.color)
        above.add_node(u.ID, label=u.label, color=u.color)
    
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
        
        self.S.supply_leaves()
        
        # top-level call on full leaf set and rho_S
        self.T = PhyloTree(self._build_gene_tree(self.L, self.S.root.children[0]))
        
        # --- post-processing of T ---
        
        # planted root
        self.T.add_planted_root()
        self.T.root.tstamp = self.S.root.tstamp
        self.T.root.color = self.S.root.ID
        
        # suppress all vertices with a single child except the planted root
        to_suppress = []
        for u in self.T.preorder():
            if u.parent and len(u.children) == 1:
                to_suppress.append((u.parent, u))
        self.T.contract(to_suppress)
        
        self.T.reconstruct_IDs()
        distance_from_timing(self.T)
        
        return self.S, self.T
        
    
    def _species_tree(self):
        
        cotree = Cotree.cotree(self.G)
        
        if not cotree:
            return False
        
        triples = set()
        
        for a, b, c in cotree.paths_of_length_2():
            
            A = self.G.nodes[a.ID]['color']
            B = self.G.nodes[b.ID]['color']
            C = self.G.nodes[c.ID]['color']
            
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
            S.add_planted_root()
            S.reconstruct_IDs()
            simulate_timing(S)
            distance_from_timing(S)
            return S
    
    
    def _build_gene_tree(self, L, u_S):
        
        u_T = PhyloTreeNode(-1, label='D',
                            color=(u_S.parent.ID, u_S.ID),
                            tstamp=u_S.tstamp + self.epsilon)
        
        # u_S is a leaf
        if not u_S.children:
            for x in L:
                leaf = PhyloTreeNode(x, label=str(x),
                                     color=self.G.nodes[x]['color'],
                                     tstamp=0.0)
                u_T.add_child(leaf)
        
        # u_S is an inner vertex
        else:
            
            # maps color to the respective child of u_S
            color_to_v_S = {}
            for v_S in u_S.children:
                for x in v_S.leaves:
                    color_to_v_S[x.ID] = v_S
                    
            # connected components C of G[L']
            for C in nx.connected_components(self.G.subgraph(L)):
                
                # choose v_S_star s.t. L(S(v_S_star)) \cap sigma(C) is non-empty
                v_S_star = color_to_v_S[self.G.nodes[next(iter(C))]['color']]
                v_T = PhyloTreeNode(-1, label='H',
                                    color=(u_S.ID, v_S_star.ID),
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


if __name__ == '__main__':
    
    import asymmetree.treeevolve as te
    import asymmetree.tools.GraphTools as gt
    from asymmetree.hgt.Fitch import undirected_fitch, rs_transfer_edges
    
    S = te.simulate_species_tree(10)
    TGT = te.simulate_dated_gene_tree(S, dupl_rate=1.0, loss_rate=0.5,
                                      hgt_rate=0.5)
    OGT = te.observable_tree(TGT)
    
    print('--- S ---\n', S.to_newick())
    print(S.to_newick(distance=False, label_inner=False))
    print('--- OGT ---\n', OGT.to_newick())
    
    ldt, above, equal = below_equal_above(OGT, S)
    fitch = undirected_fitch(OGT, rs_transfer_edges(OGT, S))
    n = ldt.order()
    print('Genes:', n, 'Total relations:', int(n * (n-1) / 2))
    print('< {}\n= {}\n> {}'.format(ldt.size(), equal.size(), above.size()))
    
    rs_scen_constr = RsScenarioConstructor(ldt)
    result = rs_scen_constr.run()
    
    if result:
        S2, T2 = result
        print('--- S2 ---\n', S2.to_newick(distance=False))
        print('--- T2 ---\n', T2.to_newick(distance=False))
        ldt2 = ldt_graph(T2, S2)
        print(ldt2.order(), ldt2.size(), gt.graphs_equal(ldt, ldt2))
        
        print('--- fitch ---')
        fitch2 = undirected_fitch(T2, rs_transfer_edges(T2, S2))
        print('Order: {} vs {}'.format(fitch.order(), fitch2.order()))
        print('Size: {} vs {}'.format(fitch.size(), fitch2.size()))
        print(gt.contingency_table(fitch, fitch2))
    else:
        print(False)