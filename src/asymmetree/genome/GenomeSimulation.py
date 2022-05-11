# -*- coding: utf-8 -*-

"""Simulation of genomes."""

import os

from tralda.datastructures.Tree import Tree

import asymmetree.treeevolve as te
import asymmetree.seqevolve as se
from asymmetree.file_io.SeqFileIO import write_alignment, write_fasta
from asymmetree.tools.Sampling import Sampler


__author__ = 'David Schaller'


class GenomeSimulator:
    
    def __init__(self, species_tree, outdir=None):
        """
        Parameters
        ----------
        species_tree : Tree
            The species tree along which gene trees are simulated.
        outdir : str, optional
            The path to a directory into which the results of the simulation
            (serialized trees, one fasta file per genome, one 'true' alignement
            per gene family). The default is None, in which case nothing is
            saved to file.
        """
        
        if not isinstance(species_tree, Tree):
            raise TypeError("species tree must be of type 'Tree'")
        
        self.S = species_tree
        self.outdir = outdir
        
        if self.outdir:
            self._check_outdir()
            self.S.serialize(self._path('species_tree.json'), mode='json')
            
        self.true_gene_trees = []
        self.pruned_gene_trees = []
        
    
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
            
    
    def simulate_gene_trees(self, n,
                            **kwargs):
        """Simulate gene trees.
        
        Simulates the 'true gene trees' which still contain the loss events as
        well as the 'pruned gene trees' in which the loss branches are
        removed.
        
        Parameters
        ----------
        n : int
            Number of gene trees to be simulated.
        kwargs : optional
            See function 'gene_trees' in 'treeevolve' subpackage or 
            documentation.
        """
        
        self.number_of_families = n
        
        self.true_gene_trees = te.gene_trees(self.S, n=n, **kwargs)
        
        self.pruned_gene_trees = [te.prune_losses(tree)
                                  for tree in self.true_gene_trees]
        
        # sequences should be emptied here if methods were called before
        if hasattr(self, 'sequence_dicts'):
            self.sequence_dicts.clear()
        
        if self.outdir:
            for i in range(n):
                filename = self._path('true_gene_trees',
                                      'gene_tree{}.json'.format(i))
                self.true_gene_trees[i].serialize(filename, mode='json')
        
        return self.true_gene_trees, self.pruned_gene_trees
                
    
    def simulate_sequences(self, subst_model,
                           root_genome=None,
                           length_distr=('constant', 200),
                           min_length=10,
                           max_length=None,
                           write_fastas=True,
                           write_alignments=True,
                           **kwargs):
        """Simulate sequences along the (pruned) gene trees.
        
        Parameters
        ----------
        subst_model : SubstModel
            Substitution model.
        root_genome : list, optional
            List of sequences for the roots of the gene trees, must contain the
            same number of str sequences as trees that were simulated. 
            The sequences must be compatible with the specified substitution 
            model (model_type='n'/'a'). The default is None, in which case new
            sequences are generated at random.
        length_distr : tuple, optional
            Distribution of the length of the root sequences if root genome is
            not supplied. See documentation for available options.
        min_length : int, optional
            Minimal length at which the distribution of lengths is truncated;
            must be less than the mean of this distribution. The default is 10.
        max_length : int, optional
            Maximal length at which the distribution of lengths is truncated;
            must be less than the mean of this distribution. The default is
            None, in which case the distribution is not truncated.
        write_fastas : bool, optional
            If True and an output directory was speci ed, write the sequences
            (one file per species) into the directory 'fasta_files' in the
            output directory. The default is True.
        write_alignments : bool, optional
            If True and an output directory was speci ed, write the true 
            alignments (one file per gene tree) into the directory 'alignments'
            in the output directory. The default is True.
        kwargs : optional
            See class 'Evolver' in the 'seqevolve' subpackage or documentation.
        """
        
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
            if isinstance(length_distr, Sampler):
                self.sampler = length_distr
            else:
                self.sampler = Sampler(length_distr,
                                       minimum=min_length,
                                       maximum=max_length,
                                       discrete=True)
        
        evolver = se.Evolver(subst_model, **kwargs)
        
        for i in range(self.number_of_families):
            
            PGT = self.pruned_gene_trees[i]
            
            if root_genome:
                evolver.evolve_along_tree(PGT, start_seq=root_genome[i].upper())
                
            else:
                evolver.evolve_along_tree(PGT, start_length=self.sampler())
                
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
        
        alg_builder = se.AlignmentBuilder(self.pruned_gene_trees[family_id],
                                          self.sequence_dicts[family_id],
                                          self.subst_model.alphabet,
                                          include_inner=False)
        
        alignment = []
        for node, sequence in alg_builder.build().items():
            label = self._compose_label(node, family_id)
            alignment.append( (label, sequence) )
            
        basename = 'alignment{}.phylip'.format(family_id)
        filename = self._path('alignments', basename)
        write_alignment(filename, alignment, alignment_format='phylip')
    
                
    def _write_fastas(self, include_inner=False):
        
        if not self.outdir:
            raise RuntimeError('no output directory specified for fasta files')
        
        # list of all species labels
        if include_inner:
            species = [v.label for v in self.S.preorder()]
        else:
            species = [v.label for v in self.S.preorder() if not v.children]
        
        # labeled sequences sorted by species
        sorted_seqs = {s: [] for s in species}
        
        for i in range(self.number_of_families):
            for node, evoseq in self.sequence_dicts[i].items():
                
                # skip inner nodes
                if not include_inner and node.children:
                    continue
                
                label = self._compose_label(node, i)
                sequence = self.subst_model.to_sequence(evoseq)
                sorted_seqs[node.reconc].append( (label, sequence) )
                
        for spec, sequences in sorted_seqs.items():
            basename = '{}.f{}a'.format(spec, self.subst_model.model_type)
            filename = self._path('fasta_files', basename)
            write_fasta(filename, sequences)
        
    
    def _compose_label(self, node, family_id):
        
        return 'fam{}gene{}spec{}'.format(family_id, node.label, node.reconc)
