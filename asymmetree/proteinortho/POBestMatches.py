# -*- coding: utf-8 -*-

import os

import networkx as nx
from Bio import SeqIO

try:
    from .POParser import parse_best_match_candidates
    from .POOutgroupFinder import OutgroupFinder
except ModuleNotFoundError:
    from proteinortho.POParser import parse_best_match_candidates
    from proteinortho.POOutgroupFinder import OutgroupFinder


class POBestMatches:
    
    def __init__(self, fasta_files, candidate_dir,
                 tree_file):
        
        self.BMG = nx.DiGraph()
        self.fasta_files = fasta_files
        self.candidate_dir = candidate_dir
        
        self.tree_file = tree_file
        
    
    def __call__(self):
        
        self._check_file_existence()
        self._build_fasta_index()
        
        self.outgroup_finder = OutgroupFinder(self.tree_file)
    
    
    def _check_file_existence(self):
        
        self.species = {}
        self.candidate_files = {}
        
        # check if all fasta files exist
        for filename in self.fasta_files:
            
            if os.path.isfile(filename):
                self.species[os.path.basename(filename)] = filename
            else:
                raise FileNotFoundError("Could not find file '{}'.".format(filename))
        
        # check if the best match candidate files exist
        for spec in self.species.keys():
            
            filename = os.path.join(self.candidate_dir, spec + ".bm_candidates")
            if os.path.isfile(filename):
                self.candidate_files[spec] = filename
            else:
                raise FileNotFoundError("Could not find file '{}'.".format(filename))
                
        # check if the tree file exists
        if self.tree_file is not None and not os.path.isfile(self.tree_file):
            raise FileNotFoundError("Could not find tree file '{}'.".format(self.tree_file))
    
    
    def _build_fasta_index(self):
        
        self.fasta_index = {}
        
        for spec, filename in self.species.items():
            self.fasta_index[spec] = SeqIO.index(filename, "fasta")
            
            for label in self.fasta_index[spec].keys():
                self.BMG.add_node("{}_{}".format(spec, label))

    
    def _get_seq_from_label(self, label, spec):
        
#        # original label in fasta file (and index) without 'SPECIES_'
#        label = label[len(spec)+1:]
        
        return self.fasta_index[spec][label].seq
    
    
    def _best_matches_for_species(self, spec_x):
        
        candidates = parse_best_match_candidates(self.candidate_files[spec_x])
        
        for x in candidates.keys():
            for spec_Y, Y in candidates[x].items():
                
                outgroups = self.outgroup_finder(spec_x, spec_Y, candidates[x])
    
    
            

if __name__ == "__main__":
    
    import glob
    
    directory = "test_files_2"
    
    fasta_files = glob.glob(os.path.join(directory, "*.faa"))
    candidate_dir = os.path.join(directory, "test.bm_candidates")
    tree_file = os.path.join(directory, "true_species_tree")
    
    po_bm = POBestMatches(fasta_files, candidate_dir, tree_file)
    po_bm()
    
    po_bm._best_matches_for_species("8.faa")