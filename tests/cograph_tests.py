# -*- coding: utf-8 -*-

import unittest
import networkx as nx

from asymmetree.cograph import Cotree, linear_cograph_detection, CographEditor
from asymmetree.tools.GraphTools import random_graph


__author__ = 'David Schaller'


class TestCographPackage(unittest.TestCase):

    
    def test_is_cograph(self):
        
        cotree = Cotree.random_cotree(100)
        cograph = cotree.to_cograph()
        
        self.assertTrue( linear_cograph_detection(cograph) )
        
        ce = CographEditor(cograph)
        ce.cograph_edit(run_number=10)
        self.assertEqual(ce.best_cost, 0)
        
        
    def test_is_no_cograph(self):
        
        graph = nx.Graph()
        graph.add_edge('a', 'b')
        graph.add_edge('b', 'c')
        graph.add_edge('c', 'd')
        
        self.assertFalse(linear_cograph_detection(graph))
        
        ce = CographEditor(graph)
        ce.cograph_edit(run_number=10)
        self.assertGreater(ce.best_cost, 0)
        
    
    def test_editing(self):
        
        graph = random_graph(100, p=0.3)
        ce = CographEditor(graph)
        ce.cograph_edit(run_number=10)
        cograph = ce.best_T.to_cograph()
        
        self.assertTrue( linear_cograph_detection(cograph) )
        

if __name__ == '__main__':
    
    unittest.main()