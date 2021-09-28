# -*- coding: utf-8 -*-

import unittest

import asymmetree.analysis as analysis

from asymmetree.tools.PhyloTreeTools import random_colored_tree


__author__ = 'David Schaller'


class TestBMG(unittest.TestCase):

    
    def test_lrt(self):
        
        N, colors = 30, 5
        repeats = 20
        
        for _ in range(repeats):
            
            tree = random_colored_tree(N, colors)
            bmg = analysis.bmg_from_tree(tree)
            
            lrt1 = analysis.lrt_from_observable_tree(tree)
            lrt2 = analysis.lrt_from_colored_graph(bmg, mincut=False)
            
            self.assertTrue( lrt1.equal_topology(lrt2) )
            

if __name__ == '__main__':
    
    unittest.main()