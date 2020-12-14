# -*- coding: utf-8 -*-

import unittest


import asymmetree.treeevolve as te
import asymmetree.hgt as hgt
import asymmetree.tools.GraphTools as gt
from asymmetree.cograph import Cotree


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
        ldt = hgt.ldt_graph(OGT, S)
        transfer_edges = hgt.rs_transfer_edges(OGT, S)
        fitch = hgt.undirected_fitch(OGT, transfer_edges)
        
        cotree = Cotree.cotree(ldt)
        
        self.assertTrue( gt.is_subgraph(ldt, fitch) and cotree )
            

if __name__ == '__main__':
    
    unittest.main()