# -*- coding: utf-8 -*-

import os

from asymmetree import PhyloTree
import asymmetree.treeevolve as te
import asymmetree.seqevolve as se


__author__ = "David Schaller"


class GenomeSimulator:
    
    def __init__(self, species_tree,
                 output_directory=None):
        
        if not isinstance(species_tree, PhyloTree):
            raise TypeError("species tree must be of type 'PhyloTree'")
        
        self.S = species_tree
        self.outdir = output_directory
        
        if self.outdir:
            self._check_outdir()
            self.species_tree.serialize(self._path('species_tree.json'),
                                        mode='json')
            
        self.true_gene_trees = []
        self.observable_gene_trees = []
        
    
    def _check_outdir(self):
        
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
        elif os.path.exists(self.outdir) and not os.path.isdir(self.outdir):
            raise FileExistsError("'{}' is not a directory".format(self.outdir))
            
    
    def _path(self, *args):
        
        return os.path.join(self.outdir, *args)
            
    
    def simulate_gene_trees(self, N, **kwargs):
        
        self.number_of_families = N
        self.true_gene_trees = []
        self.observable_gene_trees = []
        
        simulator = te.GeneTreeSimulator(self.S)
        _, autocorr_factors = te.autocorrelation_factors(self.S)
        
        for i in range(N):
            
            TGT = simulator.simulate(**kwargs)
            te.imbalance_tree(TGT, self.S,
                              baseline_rate=1.0,
                              autocorr_factors=autocorr_factors,
                              **kwargs)
            self.true_gene_trees.append(TGT)
            
            OGT = te.observable_tree(TGT)
            self.observable_gene_trees.append(OGT)
            
            if self.outdir:
                TGT.serialize(self._path('true_gene_trees',
                                         'gene_tree{}.json'.format(i)),
                              mode='json')
                
    
    def simulate_sequences(self, subst_model, root_genome=None, **kwargs):
        
        if root_genome is not None:
            if len(root_genome) != len(self.number_of_families):
                raise ValueError('no. of sequences in root genome does not'\
                                 'match no of gene families')
    
    