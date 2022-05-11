# -*- coding: utf-8 -*-

import unittest, os

from tralda.datastructures.Tree import Tree

from asymmetree.tools.PhyloTreeTools import (random_colored_tree,
                                             to_newick,
                                             parse_newick)


__author__ = 'David Schaller'


class TestTrees(unittest.TestCase):
    
    
    def test_newick(self):
        
        N, colors = 30, 5
        repeats = 20
        
        for _ in range(repeats):
            
            tree = random_colored_tree(N, colors)
            newick = to_newick(tree)
            tree2 = parse_newick(newick)
            
            label_color = [(v.label, v.reconc) for v in tree.leaves()]
            label_color2 = [(v.label, v.reconc) for v in tree2.leaves()]
            
            self.assertListEqual(label_color, label_color2)

    
    def test_serialization(self):
        
        N, colors = 30, 5
        repeats = 20
        
        for _ in range(repeats):
            
            tree = random_colored_tree(N, colors)
            
            tree1 = tree.copy()
            
            tree.serialize('testfile_tree.pickle')
            tree.serialize('testfile_tree.json')
            tree2 = Tree.load('testfile_tree.pickle')
            tree3 = Tree.load('testfile_tree.json')
            os.remove('testfile_tree.pickle')
            os.remove('testfile_tree.json')
            
            tree_nodes = [v.label for v in tree.preorder()]
            tree1_nodes = [v.label for v in tree1.preorder()]
            tree2_nodes = [v.label for v in tree2.preorder()]
            tree3_nodes = [v.label for v in tree3.preorder()]
            
            self.assertListEqual(tree_nodes, tree1_nodes)
            self.assertListEqual(tree_nodes, tree2_nodes)
            self.assertListEqual(tree_nodes, tree3_nodes)
            

if __name__ == '__main__':
    
    unittest.main()
