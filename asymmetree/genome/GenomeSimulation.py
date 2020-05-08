# -*- coding: utf-8 -*-

import os

import numpy as np

from asymmetree import PhyloTree
import asymmetree.treeevolve as te
import asymmetree.seqevolve as se
from asymmetree.file_io.SeqFileIO import write_alignment, write_fasta


__author__ = "David Schaller"


class GenomeSimulator:
    
    def __init__(self, species_tree, outdir=None):
        
        if not isinstance(species_tree, PhyloTree):
            raise TypeError("species tree must be of type 'PhyloTree'")
        
        self.S = species_tree
        self.outdir = outdir
        
        if self.outdir:
            self._check_outdir()
            self.S.serialize(self._path('species_tree.json'), mode='json')
            
        self.true_gene_trees = []
        self.observable_gene_trees = []
        
    
    def _check_outdir(self):
        
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
        elif os.path.exists(self.outdir) and not os.path.isdir(self.outdir):
            raise FileExistsError("'{}' is not a directory".format(self.outdir))
            
        for directory in ('true_gene_trees', 'fasta_files', 'alignments'):
            path = os.path.join(self.outdir, directory)
            if not os.path.exists(path):
                os.makedirs(path)
            
    
    def _path(self, *args):
        
        return os.path.join(self.outdir, *args)
            
    
    def simulate_gene_trees(self, N, **kwargs):
        
        self.number_of_families = N
        self.true_gene_trees = []
        self.observable_gene_trees = []
        
        # sequences should be emptied here if methods were called before
        if hasattr(self, 'sequence_dicts'):
            self.sequence_dicts.clear()
        
        simulator = te.GeneTreeSimulator(self.S)
        
        # autocorrelation between genes of the same or related species
        autocorr_variance = kwargs.pop('autocorr_variance', 0.0)
        _, autocorr_factors = te.autocorrelation_factors(self.S,
                                                         autocorr_variance)
        
        for i in range(self.number_of_families):
            
            TGT = simulator.simulate(**kwargs)
            te.imbalance_tree(TGT, self.S,
                              baseline_rate=1.0,
                              autocorr_factors=autocorr_factors,
                              **kwargs)
            self.true_gene_trees.append(TGT)
            
            OGT = te.observable_tree(TGT)
            self.observable_gene_trees.append(OGT)
            
            if self.outdir:
                filename = self._path('true_gene_trees',
                                      'gene_tree{}.json'.format(i))
                TGT.serialize(filename, mode='json')
                
    
    def simulate_sequences(self, subst_model,
                           root_genome=None,
                           length_distribution=('constant', 200),
                           min_length=10,
                           write_fastas=True,
                           write_alignments=True,
                           **kwargs):
        
        self.subst_model = subst_model
        
        if hasattr(self, 'sequence_dicts'):
            self.sequence_dicts.clear()
        else:
            self.sequence_dicts = []
        
        if root_genome:
            if len(root_genome) != len(self.number_of_families):
                raise ValueError('no. of sequences in root genome does not'\
                                 'match no of gene families')
        else:
            self._check_distribution(length_distribution, min_length)
        
        evolver = se.Evolver(subst_model, **kwargs)
        
        for i in range(self.number_of_families):
            
            OGT = self.observable_gene_trees[i]
            
            if root_genome:
                evolver.evolve_along_tree(OGT, start_seq=root_genome[i].upper())
                
            else:
                start_length = self._get_sequence_length()
                evolver.evolve_along_tree(OGT, start_length=start_length)
                
            self.sequence_dicts.append(evolver.sequences)
            
            # write one alignment file per gene family
            if self.outdir and write_alignments:
                self._write_alignment(i)
        
        # write one fasta file per species
        if self.outdir and write_fastas:
            self._write_fastas(include_inner=False)
            
            
    def _write_alignment(self, family_id):
        
        if not self.outdir:
            raise RuntimeError('no output directory specified for alignments')        
        
        alg_builder = se.AlignmentBuilder(self.observable_gene_trees[family_id],
                                          self.sequence_dicts[family_id],
                                          self.subst_model.alphabet,
                                          include_inner=False)
        
        alignment = []
        for node, sequence in alg_builder.build().items():
            label = self._compose_label(node, family_id)
            alignment.append( (label, sequence) )
            
        basename = 'alignment{}.phylip'.format(family_id)
        filename = self._path('alignments', basename)
        write_alignment(filename, alignment, al_format='phylip')
    
                
    def _write_fastas(self, include_inner=False):
        
        if not self.outdir:
            raise RuntimeError('no output directory specified for fasta files')
        
        # list of all species IDs
        if include_inner:
            species = [v.ID for v in self.S.preorder()]
        else:
            species = [v.ID for v in self.S.preorder() if not v.children]
        
        # labeled sequences sorted by color/species
        sorted_seqs = {s: [] for s in species}
        
        for i in range(self.number_of_families):
            for node, evoseq in self.sequence_dicts[i].items():
                
                # skip inner nodes
                if not include_inner and node.children:
                    continue
                
                label = self._compose_label(node, i)
                sequence = self.subst_model.to_sequence(evoseq)
                sorted_seqs[node.color].append( (label, sequence) )
                
        for spec, sequences in sorted_seqs.items():
            basename = '{}.f{}a'.format(spec, self.subst_model.model_type)
            filename = self._path('fasta_files', basename)
            write_fasta(filename, sequences)
                
                
    def _check_distribution(self, length_distribution, min_length):
        
        if not isinstance(min_length, int) or min_length < 0:
            raise ValueError('minimal sequence length must be an int >=0')
        self._min_length = min_length
            
        if length_distribution[0] == 'constant':
            seq_length = length_distribution[1]
            if not isinstance(seq_length, int) or min_length > seq_length:
                raise ValueError('constant sequence length must be an int >= '\
                                 'minimal sequence length')
                
            self._length_distribution = 'constant'
            self._seq_length = seq_length
            
        elif length_distribution[0] == 'gamma':
            shape = length_distribution[1]
            scale = length_distribution[2]
            
            if (not isinstance(shape, float) or shape <= 0.0 or 
                not isinstance(scale, float) or scale <= 0.0):
                raise ValueError('scale and shape parameters for gamma '\
                                 'distribution must be floats >0.0')
            
            exp_val = round(shape * scale, 3)
            if exp_val < min_length:
                raise ValueError('expected value for gamma distribution '\
                                 '({}) must be >= minimal length'.format(exp_val))
                
            self._shape = shape
            self._scale = scale
            
        else:
            raise ValueError("length distribution '{}' is not "\
                             "supported".format(length_distribution))
            
    
    def _get_sequence_length(self):
        
        if self._length_distribution == 'constant':
            return self._seq_length
        else:
            length = -1
            while length < self._min_length:
                length = round(np.random.gamma(self._shape, scale=self._scale))
            return length
        
    
    def _compose_label(self, node, family_id):
        
        return 'fam{}gene{}spec{}'.format(family_id, node.label, node.color)
            
            
    
    