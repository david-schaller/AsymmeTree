# -*- coding: utf-8 -*-

import unittest

from asymmetree.datastructures import PhyloTree
import asymmetree.best_matches as bm


__author__ = 'David Schaller'


class TestBMG(unittest.TestCase):

    
    def test_lrt(self):
        
        N, colors = 30, 5
        repeats = 20
        
        for _ in range(repeats):
            
            tree = PhyloTree.random_colored_tree(N, colors)
            bmg = bm.bmg_from_tree(tree)
            
            lrt1 = bm.lrt_from_observable_tree(tree)
            lrt2 = bm.lrt_from_colored_graph(bmg, mincut=False)
            
            self.assertTrue( lrt1.equal_topology(lrt2) )
            

if __name__ == '__main__':
    
    unittest.main()