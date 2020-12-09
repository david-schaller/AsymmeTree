# -*- coding: utf-8 -*-

"""Algorithms for building trees from rooted triple sets."""


import itertools
import networkx as nx

from asymmetree.datastructures.PhyloTree import PhyloTree, PhyloTreeNode
from asymmetree.datastructures.Partition import Partition

from asymmetree.tools.Partitioning import (Karger,
                                           greedy_bipartition,
                                           gradient_walk_bipartition)


__author__ = 'David Schaller'


def aho_graph(R, L, weighted=False, triple_weights=None):
    """Construct the auxiliary graph (Aho graph) for BUILD.
        
    Edges {a,b} are optionally weighted by the number of occurrences, resp.,
    sum of weights of triples of the form ab|x or ba|x.
    """
    
    G = nx.Graph()
    G.add_nodes_from(L)

    for a, b, c in R:
        if not G.has_edge(a, b):
            G.add_edge(a, b, weight=0)
        
        if weighted:
            if triple_weights:
                G[a][b]['weight'] += triple_weights[a, b, c]
            else:
                G[a][b]['weight'] += 1
    
    return G


def _triple_connect(G, t):
    
    G.add_edge(t[0], t[1])
    G.add_edge(t[0], t[2])
    

def mtt_partition(L, R, F):
    
    # auxiliary graph initialized as Aho graph
    G = aho_graph(R, L, weighted=False)
    
    # auxiliary partition
    P = Partition(nx.connected_components(G))
    
    if len(P) == 1:
        return P, G
    
    # aux. set of forbidden triples
    S = {t for t in F if P.separated_xy_z(*t)}
    
    # lookup of forb. triples to which u belongs
    L = {u: [] for u in L}
    for t in F:
        for u in t:
            L[u].append(t)
    
    while S:
        t = S.pop()
        _triple_connect(G, t)
        
        # merge returns the smaller of the two merged sets
        smaller_set = P.merge(t[0], t[2])
        
        # update S by traversing the L(u)
        for u in smaller_set:
            for t in L[u]:
                if t in S and not P.separated_xy_z(*t):
                    S.remove(t)
                    _triple_connect(G, t)
                elif t not in S and P.separated_xy_z(*t):
                    S.add(t)
    
    return P, G


class Build:
    """BUILD algorithm (Aho et al. 1981)."""
    
    def __init__(self, R, L, mincut=False, 
                 weighted_mincut=False, triple_weights=None):
        
        self.R = R
        self.L = L
        self.mincut = mincut
        self.weighted_mincut = weighted_mincut
        self.triple_weights = triple_weights
    
    
    def build_tree(self, return_root=False, print_info=False):
        """Build a tree displaying all triples in R if possible.
        
        Keyword arguments:
            return_root - if True, return 'PhyloTreeNode' instead of
                'PhyloTree' instance
            print_info - print information about inconsistencies
        """
        
        self.cut_value = 0
        self.cut_list = []
        self.print_info = print_info
        
        root = self._aho(self.L, self.R)
        if not root:
            if self.print_info: print('no such tree exists')
            return None
        else:
            return root if return_root else PhyloTree(root)
    
    
    def _aho(self, L, R):
        """Recursive Aho-algorithm."""
        
        # trivial cases: only one or two leaves left in L
        if len(L) == 1:
            leaf = L.pop()
            return PhyloTreeNode(leaf, label=str(leaf))
        elif len(L) == 2:
            node = PhyloTreeNode(-1)
            for _ in range(2):
                leaf = L.pop()
                node.add_child(PhyloTreeNode(leaf, label=str(leaf)))
            return node
            
        help_graph = aho_graph(R, L, weighted=self.weighted_mincut,
                                     triple_weights=self.triple_weights)
        conn_comps = self._connected_components(help_graph)
        
        # return False if less than 2 connected components
        if len(conn_comps) <= 1:
            if self.print_info: print('Connected component:\n', conn_comps)
            return False
        
        # otherwise proceed recursively
        node = PhyloTreeNode(-1)            # place new inner node
        for cc in conn_comps:
            Li = set(cc)                    # list/dictionary --> set
            Ri = []
            for t in R:                     # construct triple subset
                if Li.issuperset(t):
                    Ri.append(t)
            Ti = self._aho(Li, Ri)          # recursive call
            if not Ti:
                return False                # raise False to previous call
            else:
                node.add_child(Ti)          # add root of the subtree
   
        return node
    
    
    def _connected_components(self, aho_graph):
        """Determines the connected components of the graph.
        
        And optionally executes a min cut if there is only one component."""
        
        conn_comps = list(nx.connected_components(aho_graph))
        if (not self.mincut) or len(conn_comps) > 1:
            return conn_comps
        else:
            # Stoer–Wagner algorithm
            cut_value, partition = nx.stoer_wagner(aho_graph)
            self.cut_value += cut_value
            if len(partition[0]) < len(partition[1]):
                smaller_comp = partition[0]
            else:
                smaller_comp = partition[1]
            for edge in aho_graph.edges():
                if ((edge[0] in smaller_comp and edge[1] not in smaller_comp)
                    or
                    (edge[1] in smaller_comp and edge[0] not in smaller_comp)):
                    self.cut_list.append(edge)
            return partition
        

class Build2:
    """BUILD / MTT algorithm with minimal cost bipartition."""
    
    def __init__(self, R, L, F=None,
                 allow_inconsistency=True,
                 bipart_method='mincut',
                 cost_function=None, cost_function_args=None,
                 weighted_mincut=False, triple_weights=None):
        
        self.R = R
        self.L = L
        
        # forbidden triples --> activates MTT if non-empty
        self.F = F
        
        # allow inconsistencies or return False?
        self.allow_inconsistency = allow_inconsistency
        
        if bipart_method in ('mincut', 'karger', 'greedy', 'gradient_walk'):
            self.bipart_method = bipart_method
        else:
            raise ValueError("unknown bipartition method "\
                             "'{}'".format(bipart_method))
        
        self.cost_function = cost_function
        self.cost_function_args = cost_function_args
            
        # parameters if bipartition method is mincut
        self.weighted_mincut = weighted_mincut
        self.triple_weights = triple_weights
    
    
    def build_tree(self, return_root=False):
        """Build a tree displaying all triples in R if possible.
        
        Keyword arguments:
            return_root - if True, return 'PhyloTreeNode' instead of
                'PhyloTree' instance
        """
        
        self.total_cost = 0
        
        if self.F:
            root = self._mtt(self.L, self.R, self.F)
        else:
            root = self._aho(self.L, self.R)
            
        return root if return_root else PhyloTree(root)
    
    
    def _trivial_case(self, L):
        
        if len(L) == 1:
            leaf = L.pop()
            return PhyloTreeNode(leaf, label=str(leaf))
        
        elif len(L) == 2:
            node = PhyloTreeNode(-1)
            for _ in range(2):
                leaf = L.pop()
                node.add_child(PhyloTreeNode(leaf, label=str(leaf)))
            return node
    
    
    def _aho(self, L, R):
        """Recursive Aho-algorithm."""
        
        # trivial case: one or two leaves left in L
        if len(L) <= 2:
            return self._trivial_case(L)
            
        aux_graph = aho_graph(R, L, weighted=self.weighted_mincut,
                                    triple_weights=self.triple_weights)
        partition = list(nx.connected_components(aux_graph))
        
        if len(partition) < 2:
            if not self.allow_inconsistency:
                return False
            else:
                partition = self._bipartition(L, aux_graph)
        
        node = PhyloTreeNode(-1)            # place new inner node
        for s in partition:
            Li, Ri = set(s), []
            for t in R:                     # construct triple subset
                if Li.issuperset(t):
                    Ri.append(t)
            Ti = self._aho(Li, Ri)          # recursive call
            if not Ti:
                return False                # raise False to previous call
            else:
                node.add_child(Ti)          # add roots of the subtrees
   
        return node
    
    
    def _mtt(self, L, R, F):
        """Recursive MTT algorithm."""
        
        # trivial case: one or two leaves left in L
        if len(L) <= 2:
            return self._trivial_case(L)
        
        partition, aux_graph = mtt_partition(L, R, F)
        
        if len(partition) < 2:
            if not self.allow_inconsistency:
                return False
            else:
                partition = self._bipartition(L, aux_graph)
        
        node = PhyloTreeNode(-1)            # place new inner node
        for s in partition:
            Li, Ri, Fi = set(s), [], []
            for Xi, X in ((Ri, R), (Fi, F)):
                for t in X:
                    if Li.issuperset(t):
                        Xi.append(t)
            Ti = self._mtt(Li, Ri, Fi)      # recursive call
            if not Ti:
                return False                # raise False to previous call
            else:
                node.add_child(Ti)          # add roots of the subtrees
   
        return node
    
    
    def _bipartition(self, L, aux_graph):
        
        best_cost, best_bp = float('inf'), None
        
        if self.bipart_method == 'mincut':
            # Stoer–Wagner algorithm
            best_cost, best_bp = nx.stoer_wagner(aux_graph)
        
        elif self.bipart_method == 'karger':
            karger = Karger(aux_graph)
            
            for _, bp in karger.generate():
                cost = self.cost_function(bp, *self.cost_function_args)
                
                if cost < best_cost:
                    best_cost, best_bp = cost, bp
        
        elif self.bipart_method == 'greedy':
            
            for _ in range(5):
                cost, bp = greedy_bipartition(L, self.cost_function,
                                              args=self.cost_function_args)
                if cost < best_cost:
                    best_cost, best_bp = cost, bp
        
        elif self.bipart_method == 'gradient_walk':
            
            for _ in range(5):
                cost, bp = gradient_walk_bipartition(L, self.cost_function,
                               args=self.cost_function_args)
                if cost < best_cost:
                    best_cost, best_bp = cost, bp
        
        self.total_cost += best_cost
        return best_bp


def greedy_BUILD(R, L, triple_weights=None, return_root=False):
    """Greedy heuristic for triple consistency.
    
    Add triples one by one and checks consistency via BUILD.
    
    Keyword arguments:
        triple_weights - weights for the triples; default is None in which
            case all triples are uniformly weighted
        return_root - if True, return 'PhyloTreeNode' instead of
            'PhyloTree' instance
    """
        
    if triple_weights:
        triples = sorted(R,
                         key=lambda triple: triple_weights[triple],
                         reverse=True)
    else:
        triples = R
            
    consistent_triples = []
    root = None
    
    for t in triples:
        consistent_triples.append(t)
        build = Build(consistent_triples, L, mincut=False)
        new_root = build.build_tree(return_root=True)
        if new_root:
            root = new_root
        else:
            consistent_triples.pop()
    
    return root if return_root else PhyloTree(root)
    

def best_pair_merge_first(R, L, triple_weights=None, return_root=False):
    """Wu’s (2004) Best-Pair-Merge-First (BPMF) heuristic.
    
    Modified version by Byrka et al. (2010) and added weights.
    
    Keyword arguments:
        triple_weights - weights for the triples; default is None in which
            case all triples are uniformly weighted
        return_root - if True, return 'PhyloTreeNode' instead of
            'PhyloTree' instance
    """
    
    # initialization
    nodes = {PhyloTreeNode(leaf, label=str(leaf)): {leaf} for leaf in L}
    leaf_to_node = {}
    
    for node in nodes:
        leaf_to_node[node.ID] = node
    
    # merging
    for i in range(len(L)-1):
        
        score = {(S_i, S_j): 0
                 for S_i, S_j in itertools.combinations(nodes.keys(), 2)}
        
        for x, y, z in R:
            
            w = triple_weights[(x,y,z)] if triple_weights else 1
            
            S_i, S_j, S_k = (leaf_to_node[x],
                             leaf_to_node[y],
                             leaf_to_node[z])
            
            if (S_i is not S_j) and (S_i is not S_k) and (S_j is not S_k):
                
                if (S_i, S_j) in score:
                    score[(S_i, S_j)] += 2 * w
                else:
                    score[(S_j, S_i)] += 2 * w
                    
                if (S_i, S_k) in score:
                    score[(S_i, S_k)] -= w
                else:
                    score[(S_k, S_i)] -= w
                    
                if (S_j, S_k) in score:
                    score[(S_j, S_k)] -= w
                else:
                    score[(S_k, S_j)] -= w
        
        current_max = float('-inf')
        S_i, S_j = None, None
        
        for pair, pair_score in score.items():
            
            if pair_score > current_max:
                current_max = pair_score
                S_i, S_j = pair
        
        # create new node S_k connecting S_i and S_j
        S_k = PhyloTreeNode(-1)
        S_k.add_child(S_i)
        S_k.add_child(S_j)
        
        nodes[S_k] = nodes[S_i] | nodes[S_j]    # set union
        for leaf in nodes[S_k]:
            leaf_to_node[leaf] = S_k
        
        del nodes[S_i]
        del nodes[S_j]
    
    if len(nodes) != 1:
        raise RuntimeError('more than 1 node left')
    
    root = next(iter(nodes))
    return root if return_root else PhyloTree(root)