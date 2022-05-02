# -*- coding: utf-8 -*-

import unittest

import asymmetree.treeevolve as te


__author__ = 'David Schaller'


class TestTreeEvolvePackage(unittest.TestCase):

    
    def test_species_tree(self):
        
        N = 30
        
        for model in ('innovation', 'yule', 'BDP', 'EBDP'):
            
            species_tree = te.simulate_species_tree(N, model=model,
                                                    contraction_proportion=0.2,
                                                    contraction_bias='inverse',
                                                    bias_strength=1.5)
            
            self.assertTrue(species_tree._assert_integrity())
            
            leaves = [l for l in species_tree.leaves() if l.event != 'L']
            self.assertEqual(len(leaves), N)
            
            
    def test_no_extinction(self):
        
        N = 10
        repeats = 20
        
        for _ in range(repeats):
            
            species_tree = te.simulate_species_tree(N, model='innovation',
                                                    contraction_probability=0.2)
            
            gene_tree = te.simulate_dated_gene_tree(species_tree,
                            dupl_rate=1.0, loss_rate=1.0, hgt_rate=0.5,
                            prohibit_extinction='per_species')
            
            # check that there is no extinction in any species
            color_dict = {l.label: [] for l in species_tree.preorder()
                          if not l.children and l.event != 'L'}
            
            for v in gene_tree.preorder():
                if not v.children and v.event != 'L':
                    color_dict[v.color].append(v.label)
                    
            for leaf_list in color_dict.values():
                self.assertTrue(leaf_list)
                
            gene_tree2 = te.simulate_dated_gene_tree(species_tree,
                             dupl_rate=1.0, loss_rate=1.0, hgt_rate=0.5,
                             prohibit_extinction='per_family')
            
            # check that there is no extinction in all species
            self.assertTrue([l for l in gene_tree2.leaves()])
            

if __name__ == '__main__':
    
    unittest.main()