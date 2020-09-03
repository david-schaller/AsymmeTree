# -*- coding: utf-8 -*-

"""Algorithms for building trees from rooted triple sets."""


import itertools
import networkx as nx

from asymmetree.datastructures.PhyloTree import PhyloTree, PhyloTreeNode


__author__ = 'David Schaller'


def aho_graph(R, L, weighted=False, triple_weights=None):
    """Construct the auxiliary graph (Aho graph) for BUILD.
        
    Edges {a,b} are optionally weighted by the number of occurrences, resp.,
    sum of weights of triples of the form ab|x or ba|x.
    """
    
    G = nx.Graph()
    G.add_nodes_from(L)

    for t1, t2, t3 in R:
        if not G.has_edge(t1, t2):
            G.add_edge(t1, t2, weight=0)
        
        if weighted:
            if triple_weights:
                G[t1][t2]['weight'] += triple_weights[t1, t2, t3]
            else:
                G[t1][t2]['weight'] += 1
    
    return G


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
        
        # trivial case: only one leaf left in L
        if len(L) == 1:
            leaf = L.pop()
            return PhyloTreeNode(leaf, label=str(leaf))
            
        help_graph = aho_graph(R, L, weighted=self.weighted_mincut,
                                     triple_weights=self.triple_weights)
        conn_comps = self._connected_components(help_graph)
        
        # return False if less than 2 connected components
        if len(conn_comps) <= 1:
            if self.print_info: print('Connected component:\n', conn_comps)
            return False
        
        # otherwise proceed recursively
        child_nodes = []
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
                child_nodes.append(Ti)
                
        node = PhyloTreeNode(-1)            # place new inner node
        for Ti in child_nodes:
            node.add_child(Ti)              # add roots of the subtrees
   
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