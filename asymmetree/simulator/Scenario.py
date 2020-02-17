# -*- coding: utf-8 -*-

"""
Scenario.

Wrapper class for species and gene tree. Compute statistics, BMG/RBMG.
"""

import networkx as nx

from asymmetree.simulator import TreeSimulator as ts
from asymmetree.simulator.DistanceMatrix import distance_matrix
from asymmetree.best_matches import TrueBMG


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


class Scenario:
    
    def __init__(self, S, TGT, DLH_rates, OGT=None):
        self.S = S
        self.TGT = TGT
        self.DLH_rates = DLH_rates
        
        self.OGT = OGT if OGT else ts.observable_tree(TGT)
        
        self._genes()
        self._count_events()
        self._sort_species_to_subtrees()
        self._sort_genes_to_species()
        
        self.TOG = TrueBMG.true_orthology_graph(self.OGT)
        self.BMG, self.RBMG = TrueBMG.best_match_graphs(self.OGT)
        
        self.BMG_subtrees, self.RBMG_subtrees = self.reduce_to_subtrees(self.BMG, self.RBMG)
        
        
    def _genes(self):
        
        self.OGT.supply_leaves()
        color_dict = {}
        for gene in self.OGT.root.leaves:
            if gene.color not in color_dict:
                color_dict[gene.color] = []
            color_dict[gene.color].append(gene)
        self.genes = []                                 # color-sorted gene list
        for color, gene_list in color_dict.items():
            for gene in gene_list:
                self.genes.append(gene)
                                                        # maps gene to index in matrix
        self.gene_index = {gene: i for i, gene in enumerate(self.genes)}
    
    
    def _count_events(self):
        
        self.event_counts = [self.S.number_of_species, 0, 0, 0, 0]    # count S / D / L / H /
                                                                      # ancient duplications                
        for v in self.TGT.preorder():
            if v.label == "D":
                self.event_counts[1] += 1
                if self.S.root.ID == v.color[0]:
                    self.event_counts[4] += 1
            elif v.label == "*":
                self.event_counts[2] += 1
            elif v.label == "H":
                self.event_counts[3] += 1
    
    
    def _sort_species_to_subtrees(self):
        
        # sort the species into the subtrees of the first speciation event of S
        self.subtree_list, self.subtree_index = [], {}
        for u in self.S.preorder():                          # exclude root of planted tree
            if len(u.children) > 1:
                for i in range(len(u.children)):
                    self.subtree_list.append([gene.ID for gene in self.S.preorder(node=u.children[i]) 
                                         if len(gene.children) == 0])
                    self.subtree_index.update({item: i for item in self.subtree_list[i]})
                break
        
        # max. number of outgroup genes (in T) for each subtree of S
        self.outgroup_dict = {i: [] for i in range(len(self.subtree_list))}
        for gene in self.genes:
            for i in self.outgroup_dict.keys():
                if self.subtree_index[gene.color] != i:
                    self.outgroup_dict[i].append(gene)
    
    
    def _sort_genes_to_species(self):
        
        self.color_dict = {item: [] for item in self.subtree_index.keys()}
        for gene in self.genes:
            self.color_dict[gene.color].append(gene)
    
    
    def get_data(self):
        
        return (self.S, self.genes, self.gene_index,
                self.subtree_list, self.subtree_index,
                self.outgroup_dict, self.color_dict)
    
    
    def rates_and_counts(self):
        """Return event rates and counts."""
        
        return list(self.DLH_rates) + self.event_counts
    
    
    def reduce_to_subtrees(self, full_BMG, full_RBMG):
        """Return a subgraph of the true RBMG with edges {u,v} for which the
        corresponding species (colors) are in the same subtree of root(S)."""
        BMG_subtrees = nx.DiGraph()
        RBMG_subtrees = nx.Graph()
        
        for G_subtrees, true_G in [(BMG_subtrees, full_BMG),
                                   (RBMG_subtrees, full_RBMG)]:
            G_subtrees.add_nodes_from(true_G.nodes(data=True))
            for u, v in true_G.edges:
                u_col = true_G.nodes[u]['color']
                v_col = true_G.nodes[v]['color']
                
                if self.subtree_index[u_col] == self.subtree_index[v_col]:
                    G_subtrees.add_edge(u, v)
        
        return BMG_subtrees, RBMG_subtrees
    
    
    def get_distance_matrix(self):
        _, _, D = distance_matrix(self.OGT, leaves=self.genes,
                                  leaf_index=self.gene_index)
        return D
    
    
    def possible_edges_BMG(self):
        """Return the number of possible edges in the BMG i.e. for gene pairs for
        which an outgroup is available in S."""
        
        counts = [0 for i in range(len(self.subtree_list))]
        
        for g in self.genes:
            counts[self.subtree_index[g.color]] += 1
        
        possible_edges = 0
        for count in counts:
            possible_edges += (count * (count-1))
        
        return possible_edges