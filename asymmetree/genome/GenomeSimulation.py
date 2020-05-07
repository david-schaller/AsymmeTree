# -*- coding: utf-8 -*-

import os

from asymmetree import PhyloTree
import asymmetree.treeevolve as te
import asymmetree.seqevolve as se


class GenomeSimulator:
    
    def __init__(self, species_tree,
                 output_directory=None):
        
        if not isinstance(species_tree, PhyloTree):
            raise TypeError("Species tree must be of type 'PhyloTree'!")
        
        self.species_tree = species_tree
        self.outdir = output_directory
        
        if self.outdir:
            self._check_outdir()
            self.species_tree
            
        self.true_gene_trees = []
        self.observable_gene_trees = []
        
    
    def _check_outdir(self):
        
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
        elif os.path.exists(self.outdir) and not os.path.isdir(self.outdir):
            raise FileExistsError("'{}' is not a directory".format(self.outdir))
            
    
    def simulate_gene_trees(self, N, **kwargs):
        
        simulator = te.GeneTreeSimulator(self.species_tree)
        
        for i in range(N):
            
            true_gene_tree = simulator.simulate()
    
    