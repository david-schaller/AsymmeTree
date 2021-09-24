# -*- coding: utf-8 -*-

import unittest

from tralda.datastructures.Tree import TreeNode

import asymmetree.seqevolve as se
import asymmetree.treeevolve as te


__author__ = 'David Schaller'


class TestSeqEvolvePackage(unittest.TestCase):
    
    
    def test_custom_matrix(self):
        
        with self.assertRaises(ValueError):
            se.SubstModel('a', 'CUSTOM')
        
        subst_model = se.SubstModel('a', 'CUSTOM', filename='custom.paml')
        
        self.assertEqual(subst_model.Q.shape, (20, 20))
        
        
    def test_simulation(self):
        
        species_tree = te.simulate_species_tree(10, model='innovation')

        subst_model = se.SubstModel('a', 'JTT')
        indel_model = se.IndelModel(0.01, 0.01, length_distr=('zipf', 1.821))
        het_model = se.HetModel(2.0, classes=10, invariant=0.1)
        
        evolver = se.Evolver(subst_model,
                             indel_model=indel_model,
                             het_model=het_model)
        evolver.evolve_along_tree(species_tree, start_length=150)
        
        for node, sequence in evolver.sequences.items():
            self.assertIsInstance(node, TreeNode)
            self.assertIsInstance(sequence, se.EvoSeq)
        
        alignment_length = -1
        for node, aligned_seq in evolver.true_alignment().items():
            self.assertIsInstance(node, TreeNode)
            self.assertIsInstance(aligned_seq, str)
            
            # test whether all aligned sequences have the same length
            if alignment_length == -1:
                alignment_length = len(aligned_seq)
            else:
                self.assertEqual(len(aligned_seq), alignment_length)
            

if __name__ == '__main__':
    
    unittest.main()