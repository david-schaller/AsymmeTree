# -*- coding: utf-8 -*-

import unittest

import networkx as nx

from tralda.tools.GraphTools import is_subgraph
from tralda.cograph import to_cotree

import asymmetree.treeevolve as te
import asymmetree.analysis as analysis


__author__ = 'David Schaller'


class TestHGT(unittest.TestCase):

    
    def test_ldt_fitch(self):
        
        S = te.species_tree_n_age(20, 1.0, model='yule')
    
        # true gene tree (with losses)
        TGT = te.dated_gene_tree(S, dupl_rate=1.0,
                                 loss_rate=0.5,
                                 hgt_rate=0.2)
        
        # pruned gene tree
        PGT = te.prune_losses(TGT)
        
        # finally we can extract the LDT and Fitch graph
        ldt = analysis.ldt_graph(PGT, S)
        transfer_edges = analysis.rs_transfer_edges(PGT, S)
        fitch = analysis.undirected_fitch(PGT, transfer_edges)
        
        cotree = to_cotree(ldt)
        
        self.assertTrue( is_subgraph(ldt, fitch) and cotree )
    
    
    def test_replacing_hgt(self):
        
        N = 20
        
        S = te.species_tree_n_age(N, 1.0, model='yule')
    
        # true gene tree
        TGT = te.dated_gene_tree(S, dupl_rate=0.0,
                                 loss_rate=0.0,
                                 hgt_rate=1.0,
                                 prohibit_extinction='per_species',
                                 replace_prob=1.0,)
        
        # pruned gene tree
        PGT = te.prune_losses(TGT)
        
        leaves = [v for v in PGT.leaves()]
        reconciliations = {v.reconc for v in leaves}
        
        # print(TGT.to_newick())
        # print(PGT.to_newick())
        
        self.assertTrue(len(reconciliations) == N and len(leaves) == N)
        
    
    def test_transfer_distance_bias(self):
        
        n = 20
        
        S = te.species_tree_n_age(n, 1.0, model='yule')
    
        # true gene tree (with losses)
        TGT = te.dated_gene_tree(S, dupl_rate=0.5,
                                 loss_rate=0.5,
                                 hgt_rate=1.0,
                                 prohibit_extinction='per_species',
                                 replace_prob=0.5,
                                 transfer_distance_bias='inverse')
        
        # pruned gene tree
        PGT = te.prune_losses(TGT)
        
        leaves = [v for v in PGT.leaves()]
        reconciliations = {v.reconc for v in leaves}
        
        self.assertTrue(len(reconciliations) == n)
        
    
    def test_rs_edges(self):
        
        S = te.species_tree_n_age(10, 1.0, model='yule')
        TGT = te.dated_gene_tree(S, dupl_rate=1.0, loss_rate=0.5, hgt_rate=0.5)
        PGT = te.prune_losses(TGT)
        
        transf1 = analysis.true_transfer_edges(PGT)
        transf2 = analysis.rs_transfer_edges(PGT, S)
        
        self.assertTrue( transf1.issuperset(transf2) )
    
    
    def test_rs_fitch_counterexamle(self):
        
        G = nx.Graph()
        G.add_node('a', color=1)
        G.add_node('a2', color=1)
        G.add_node('b', color=2)
        G.add_node('b2', color=2)
        G.add_edges_from([('a', 'a2'), ('a', 'b2'), ('b', 'a2'), ('b', 'b2')])
        
        self.assertFalse( analysis.is_rs_fitch(G, color_set=[1,2]) )
            

if __name__ == '__main__':
    
    unittest.main()
