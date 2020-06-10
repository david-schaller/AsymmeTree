# -*- coding: utf-8 -*-

import unittest

import asymmetree.treeevolve as te


__author__ = 'David Schaller'


class TestTreeEvolvePackage(unittest.TestCase):

    
    def test_species_tree(self):
        
        N = 30
        
        for model in ('innovation', 'yule', 'BDP', 'EBDP'):
            
            species_tree = te.simulate_species_tree(N, model=model,
                                                    non_binary=0.2)
            
            self.assertTrue(species_tree._assert_integrity())
            
            leaves = species_tree.supply_leaves(exclude_losses=True)
            self.assertEqual(len(leaves), N)
            
            
    def test_no_extinction(self):
        
        N = 10
        repeats = 20
        
        for _ in range(repeats):
            
            species_tree = te.simulate_species_tree(N, model='innovation',
                                                    non_binary=0.2)
            
            gene_tree = te.simulate_dated_gene_tree(species_tree,
                            DLH_rates=(1.0, 1.0, 0.5),
                            prohibit_extinction='per_species')
            
            # check that there is no extinction in any species
            color_dict = {l.ID: [] for l in species_tree.preorder()
                          if not l.children}
            
            for v in gene_tree.preorder():
                if not v.children and not v.is_loss():
                    color_dict[v.color].append(v.ID)
                    
            for leaf_list in color_dict.values():
                self.assertTrue(leaf_list)
                
            gene_tree2 = te.simulate_dated_gene_tree(species_tree,
                             DLH_rates=(1.0, 1.0, 0.5),
                             prohibit_extinction='per_family')
            
            # check that there is no extinction in all species
            self.assertTrue(gene_tree2.supply_leaves())
            

if __name__ == '__main__':
    
    unittest.main()