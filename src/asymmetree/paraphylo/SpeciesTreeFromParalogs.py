# -*- coding: utf-8 -*-

import itertools

import networkx as nx

from asymmetree.cograph.CographEditor import CographEditor
from asymmetree.datastructures.PhyloTree import PhyloTree, PhyloTreeNode
from asymmetree.tools.Build import Build


__author__ = 'David Schaller'


class TreeReconstructor:
    
    
    def __init__(self, edit_runs_per_cc=10, cotree_mode='best'):
        
        self.R = {}             # (weighted) species triples
        self.L = set()          # set of species leaves
        
        self._edit_runs_per_cc = edit_runs_per_cc   # runs of cograph editing heuristic
                                                    # per connected component
        self.S = None           # the reconstructed species tree
        
        if cotree_mode == 'best':                   # only use the best cotree
            self._cotree_mode = 'best'
        elif cotree_mode == 'all':                  # use all cotrees
            self._cotree_mode = 'all'
        else:
            raise ValueError("invalid argument '{}' for cotree usage".format(cotree_mode))
        
        
    def add_ortho_graph(self, ortho_graph):
        
        for cc in nx.connected_components(ortho_graph):
            
            G = nx.Graph()
            color_dict = {}
            
            for u in cc:
                G.add_node(u)
                for v in ortho_graph.neighbors(u):
                    G.add_edge(u, v)
                u_color = ortho_graph.nodes[u]['color']
                color_dict[u] = u_color
                self.L.add(u_color)
                
            # apply minimal cograph editing (heuristic)
            ce = CographEditor(G)
            cotree = ce.cograph_edit(run_number=self._edit_runs_per_cc)
            
            # add the informative triples (root is a speciation)
            if self._cotree_mode == 'best':
                self._informative_triples(cotree, color_dict)
                
            elif self._cotree_mode == 'all':
                weight_per_tree = 1 / len(ce.cotrees)
                for cotree in ce.cotrees:
                    self._informative_triples(cotree, color_dict,
                                              weight=weight_per_tree)
            
            
    def _informative_triples(self, cotree, color_dict, weight=1.0):
        """Add informative species triples from a (co)tree to R."""
        
        cotree.supply_leaves()
        
        for u in cotree.preorder():
            if u.label == 'parallel' or not u.children:
                continue
            
            for v1, v2 in itertools.permutations(u.children, 2):
                for x, y in itertools.combinations(v1.leaves, 2):
                    for z in v2.leaves:
                        X, Y, Z = (color_dict[x.ID],
                                   color_dict[y.ID],
                                   color_dict[z.ID])
                        
                        if X != Y and X != Z and Y != Z:
                            self._add_triple(X, Y, Z, weight)
                            
        
    def _add_triple(self, a, b, c, weight):
        """Add a triple ab|c (= ba|c)."""
        
        if a <= b:
            triple = (a, b, c)
        else:
            triple = (b, a, c)
            
        if triple in self.R:
            self.R[triple] += weight
        else:
            self.R[triple] = weight
            
            
    def build_species_tree(self, mode='bpmf', weighted=True):
        """Build a species tree from the informative triples."""
        
        if mode.lower() == 'bpmf':
            root = self._BPMF(weighted=weighted)
        elif mode.lower() in ('mincut', 'min-cut'):
            build = Build(self.R, self.L, mincut=True,
                          weighted_mincut=weighted, triple_weights=self.R)
            root = build.build_tree(return_root=True)
        elif mode.lower() == 'greedy':
            root = self._GREEDY(weighted=weighted)
        else:
            raise ValueError("mode '{}' is not valid".format(mode))
            
        if not root:
            raise RuntimeError('could not build a species tree')
        
        self.S = PhyloTree(root)
        
        return self.S
        
    
    def _BPMF(self, weighted=True):
        """Wu’s Best-Pair-Merge-First heuristic.
        
        Modified version by Byrka et al. 2010 and added weights."""
        
        # initialization
        nodes = {PhyloTreeNode(leaf, label=str(leaf)): {leaf}
                 for leaf in self.L}
        leaf_to_node = {}
        
        for node in nodes:
            leaf_to_node[node.ID] = node
        
        # merging
        for i in range(len(self.L)-1):
            
            score = {(S_i, S_j): 0
                     for S_i, S_j in itertools.combinations(nodes.keys(), 2)}
            
            for x, y, z in self.R.keys():
                
                w = self.R[(x,y,z)]
                
                S_i, S_j, S_k = (leaf_to_node[x],
                                 leaf_to_node[y],
                                 leaf_to_node[z])
                
                if (S_i is not S_j) and (S_i is not S_k) and (S_j is not S_k):
                    
                    if (S_i, S_j) in score:
                        score[(S_i, S_j)] += 2 * w if weighted else 2
                    else:
                        score[(S_j, S_i)] += 2 * w if weighted else 2
                        
                    if (S_i, S_k) in score:
                        score[(S_i, S_k)] -= 1 * w if weighted else 1
                    else:
                        score[(S_k, S_i)] -= 1 * w if weighted else 1
                        
                    if (S_j, S_k) in score:
                        score[(S_j, S_k)] -= 1 * w if weighted else 1
                    else:
                        score[(S_k, S_j)] -= 1 * w if weighted else 1
            
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
            
        assert len(nodes) == 1, 'more than 1 node left'
        
        return next(iter(nodes))
    
    
    def _GREEDY(self, weighted=True):
        
        if weighted:
            triples = sorted(self.R.keys(),
                             key=lambda triple: self.R[triple],
                             reverse=True)
        else:
            triples = self.R.keys()
                
        consistent_triples = []
        root = None
        
        for t in triples:
            consistent_triples.append(t)
            build = Build(consistent_triples, self.L, mincut=False)
            new_root = build.build_tree(return_root=True)
            if new_root:
                root = new_root
            else:
                consistent_triples.pop()
        
        return root
    
    
    def _max_consistent_triple_set(self):
        
        if self.S is None:
            raise RuntimeError('species tree has not been built yet')
            
        R_total = self.S.get_triples()      # all triples of the tree
        R_max_cons = {}                     # max. consistent subset of R (i.e. heuristic)
        
        for l1, l2, l3 in R_total:
            a = l1.ID
            b = l2.ID
            c = l3.ID
            
            if a <= b:
                triple = (a, b, c)
            else:
                triple = (b, a, c)
                
            if triple in self.R:
                R_max_cons[triple] = self.R[triple]
        
        return R_max_cons
    
    
    def _total_support(self, R_max_cons):
        """Compute the total support s of the tree as in ParaPhylo."""
        
        numerator = 0
        denominator = 0
        
        for triple, weight in R_max_cons.items():
            
            numerator += weight
            denominator += weight
            
            a, b, c = triple
            triple2 = (a, c, b) if a <= c else (c, a, b)
            triple3 = (b, c, a) if b <= c else (c, b, a)
            
            if triple2 in self.R:
                denominator += self.R[triple2]
            if triple3 in self.R:
                denominator += self.R[triple3]
        
        support = numerator / denominator if denominator > 0 else 0
        return support
    
    
    def _subtree_support(self):
        """Compute the subtree support s_v for the inner nodes as in ParaPhylo."""
        
        support = {}                        # internal node v --> support for T(v)
        
        self.S.supply_leaves()
        all_leaves = set(self.S.root.leaves)
        
        for v in self.S.preorder():
            if not v.children or v is self.S.root:
                continue
            
            numerator = 0
            denominator = 0
            outspecies = all_leaves.difference(v.leaves)
            
            for l1, l2 in itertools.combinations(v.leaves, 2):
                for l3 in outspecies:
                    a = l1.ID
                    b = l2.ID
                    c = l3.ID
                    triple1 = (a, b, c) if a <= b else (b, a, c)
                    triple2 = (a, c, b) if a <= c else (c, a, b)
                    triple3 = (b, c, a) if b <= c else (c, b, a)
                    
                    if triple1 in self.R:
                        numerator += self.R[triple1]
                        denominator += self.R[triple1]
                    if triple2 in self.R:
                        denominator += self.R[triple2]
                    if triple3 in self.R:
                        denominator += self.R[triple3]
                        
            support[v] = numerator / denominator if denominator > 0 else 0
        
        return support
    
    
    def newick_with_support(self, v=None, supports=None):
        """Recursive PhyloTree --> Newick (str) function."""
        
        if v is None:
            supports = self._subtree_support()
            supports[self.S.root] = self._total_support(self._max_consistent_triple_set())
            return self.newick_with_support(v=self.S.root, supports=supports) + ';'
        
        elif not v.children:
            return str(v.label)
        
        else:
            s = ''
            for child in v.children:
                s += self.newick_with_support(v=child, supports=supports) + ','
            return '({}){}'.format(s[:-1], supports[v])