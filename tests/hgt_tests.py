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
        
        S = te.simulate_species_tree(20, model='innovation')
    
        # true gene tree (with losses)
        TGT = te.simulate_dated_gene_tree(S, dupl_rate=1.0,
                                             loss_rate=0.5,
                                             hgt_rate=0.2)
        
        # observable gene tree
        OGT = te.observable_tree(TGT)
        
        # finally we can extract the LDT and Fitch graph
        ldt = analysis.ldt_graph(OGT, S)
        transfer_edges = analysis.rs_transfer_edges(OGT, S)
        fitch = analysis.undirected_fitch(OGT, transfer_edges)
        
        cotree = to_cotree(ldt)
        
        self.assertTrue( is_subgraph(ldt, fitch) and cotree )
    
    
    def test_replacing_hgt(self):
        
        N = 20
        
        S = te.simulate_species_tree(N, model='innovation')
    
        # true gene tree
        TGT = te.simulate_dated_gene_tree(S, dupl_rate=0.0,
                                             loss_rate=0.0,
                                             hgt_rate=1.0,
                                             prohibit_extinction='per_species',
                                             replace_prob=1.0,)
        
        # observable gene tree
        OGT = te.observable_tree(TGT)
        
        leaves = [v for v in OGT.leaves()]
        colors = {v.color for v in leaves}
        
        # print(TGT.to_newick())
        # print(OGT.to_newick())
        
        self.assertTrue(len(colors) == N and len(leaves) == N)
        
    
    def test_transfer_distance_bias(self):
        
        N = 20
        
        S = te.simulate_species_tree(N, model='innovation')
    
        # true gene tree (with losses)
        TGT = te.simulate_dated_gene_tree(S, dupl_rate=0.5,
                                             loss_rate=0.5,
                                             hgt_rate=1.0,
                                             prohibit_extinction='per_species',
                                             replace_prob=0.5,
                                             transfer_distance_bias='inverse')
        
        # observable gene tree
        OGT = te.observable_tree(TGT)
        
        leaves = [v for v in OGT.leaves()]
        colors = {v.color for v in leaves}
        
        self.assertTrue(len(colors) == N)
        
    
    def test_rs_edges(self):
        
        S = te.simulate_species_tree(10)
        TGT = te.simulate_dated_gene_tree(S, dupl_rate=1.0, loss_rate=0.5,
                                          hgt_rate=0.5)
        OGT = te.observable_tree(TGT)
        
        transf1 = analysis.true_transfer_edges(OGT)
        transf2 = analysis.rs_transfer_edges(OGT, S)
        
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