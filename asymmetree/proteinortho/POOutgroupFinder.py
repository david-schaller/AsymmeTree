# -*- coding: utf-8 -*-

import itertools

import numpy as np

from asymmetree.tools.PhyloTree import PhyloTree, PhyloTreeNode


class OutgroupFinder:
    
    def __init__(self, tree_file):
        
        self.tree_file = tree_file
        
        self._prepare()
        
        
    def _prepare(self):
        
        self._parse_tree()
        
        self._compute_lcas_in_S()
        self._compute_outgroup_species()
      
        
    def _parse_tree(self):
        
        with open(self.tree_file, "r") as f:
            newick = f.readline()
        
        self.S = PhyloTree.parse_newick(newick)
        
        
    def _compute_lcas_in_S(self):
        
        self.S.supply_leaves()
        N = len(self.S.root.leaves)
        self.leaf_index = {leaf.label: i for i, leaf in enumerate(self.S.root.leaves)}
        
        self.lca = np.zeros((N,N), dtype=PhyloTreeNode)
        self.subtree_species = {}
        
        for v in self.S.preorder():
            self.subtree_species[v] = [leaf.label for leaf in v.leaves]
            
            if not v.children:
                self.lca[self.leaf_index[v.label],self.leaf_index[v.label]] = v
            elif len(v.children) >= 2:
                for c1, c2 in itertools.combinations(v.children, 2):
                    for x in c1.leaves:
                        x_index = self.leaf_index[x.label]
                        for y in c2.leaves:
                            y_index = self.leaf_index[y.label]
                            self.lca[x_index, y_index] = v
                            self.lca[y_index, x_index] = v
    
    
    def _compute_outgroup_species(self):
        
        self.lca_to_outgroups = {key: set() for key in self.subtree_species.keys()}
        
        for v in self.S.preorder():
            if v.parent:
                self.lca_to_outgroups[v].update(self.lca_to_outgroups[v.parent])
                
            for c1, c2 in itertools.permutations(v.children, 2):
                self.lca_to_outgroups[c1].update(self.subtree_species[c2])
                
    
    def __call__(self, spec_x, spec_Y, candidates, limit=10):
        
        # set of candidates might be empty (?)
        if not candidates:
            return False
        
        outgroup_spec = self.lca_to_outgroups[self.lca[self.leaf_index[spec_x],
                                                       self.leaf_index[spec_Y]]]
        
        # set of potential outgroup species can be empty (when lca is the root)
        if not outgroup_spec:
            return False
        
        available_outgroup_spec = outgroup_spec.intersection(candidates)
        
        total_outgroup_nr = 0
        for spec in available_outgroup_spec:
            total_outgroup_nr += len(candidates[spec])
            
        # set of actual outgroup genes can be empty
        if total_outgroup_nr == 0:
            return False
        
        outgroups = []
        
        # pick one outgroup from each species in each round
        i, outgroup_counter, stop = 0, 0, min(limit, total_outgroup_nr)
        while outgroup_counter < stop:
            
            for spec in available_outgroup_spec:
                if i < len(candidates[spec]):
                    outgroups.append( candidates[spec][i] )
                    outgroup_counter += 1
                    if outgroup_counter == stop:
                        break
            i += 1
        
        return outgroups