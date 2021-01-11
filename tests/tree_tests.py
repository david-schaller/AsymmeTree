# -*- coding: utf-8 -*-

import unittest, os

from asymmetree.datastructures import PhyloTree


__author__ = 'David Schaller'


class TestTrees(unittest.TestCase):
    
    
    def test_newick(self):
        
        N, colors = 30, 5
        repeats = 20
        
        for _ in range(repeats):
            
            tree = PhyloTree.random_colored_tree(N, colors)
            newick = tree.to_newick()
            tree2 = PhyloTree.parse_newick(newick)
            
            # colors of tree must be converted to 'str'
            label_color = [(v.label, str(v.color)) for v in tree.supply_leaves()]
            label_color2 = [(v.label, v.color) for v in tree2.supply_leaves()]
            
            self.assertListEqual(label_color, label_color2)

    
    def test_serialization(self):
        
        N, colors = 30, 5
        repeats = 20
        
        for _ in range(repeats):
            
            tree = PhyloTree.random_colored_tree(N, colors)
            
            tree1 = tree.copy()
            
            tree.serialize('testfile_tree.pickle')
            tree.serialize('testfile_tree.json')
            tree2 = PhyloTree.load('testfile_tree.pickle')
            tree3 = PhyloTree.load('testfile_tree.json')
            os.remove('testfile_tree.pickle')
            os.remove('testfile_tree.json')
            
            tree_nodes = [v.ID for v in tree.preorder()]
            tree1_nodes = [v.ID for v in tree1.preorder()]
            tree2_nodes = [v.ID for v in tree2.preorder()]
            tree3_nodes = [v.ID for v in tree3.preorder()]
            
            self.assertListEqual(tree_nodes, tree1_nodes)
            self.assertListEqual(tree_nodes, tree2_nodes)
            self.assertListEqual(tree_nodes, tree3_nodes)
            

if __name__ == '__main__':
    
    unittest.main()