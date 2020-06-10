# -*- coding: utf-8 -*-

import unittest

from asymmetree import PhyloTree
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
            lrt_constr = bm.LRTConstructor(bmg, mincut=False)
            lrt2 = lrt_constr.build_tree()
            
            self.assertTrue( lrt1.compare_topology(lrt2) )
            

if __name__ == '__main__':
    
    unittest.main()