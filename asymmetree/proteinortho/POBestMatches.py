# -*- coding: utf-8 -*-

import os, itertools, random

import numpy as np
import networkx as nx
from Bio import SeqIO


from asymmetree.POParser import parse_best_match_candidates
from asymmetree.POOutgroupFinder import OutgroupFinder
from asymmetree.PODistanceCalculator import distance_2seqs

class POBestMatches:
    
    def __init__(self, fasta_files, candidate_dir,
                 tree_file,
                 epsilon=0.5,
                 voting_mode="majority"):
        
        self.BMG = nx.DiGraph()
        self.fasta_files = fasta_files
        self.candidate_dir = candidate_dir
        
        self.tree_file = tree_file
        
        self.epsilon = epsilon
        
        if voting_mode in ("majority", "weighted"):
            self.voting_mode = voting_mode
        else:
            raise ValueError(f"Invalid voting mode for quartets '{voting_mode}'")
        
    
    def __call__(self):
        
        self._check_file_existence()
        self._build_fasta_index()
        
        self.distance_cache = {}
        
        self.outgroup_finder = OutgroupFinder(self.tree_file)
        
        for spec in self.species.keys():
            print("----- Species '{}' -----".format(spec))
            self._best_matches_for_species(spec)
            
        return self.BMG
    
    
    def _check_file_existence(self):
        
        self.species = {}
        self.candidate_files = {}
        
        # check if all fasta files exist
        for filename in self.fasta_files:
            
            if os.path.isfile(filename):
                self.species[os.path.basename(filename)] = filename
            else:
                raise FileNotFoundError(f"Could not find file '{filename}'.")
        
        # check if the best match candidate files exist
        for spec in self.species.keys():
            
            filename = os.path.join(self.candidate_dir, spec + ".bm_candidates")
            if os.path.isfile(filename):
                self.candidate_files[spec] = filename
            else:
                raise FileNotFoundError(f"Could not find file '{filename}'.")
                
        # check if the tree file exists
        if self.tree_file is not None and not os.path.isfile(self.tree_file):
            raise FileNotFoundError(f"Could not find tree file '{self.tree_file}'.")
    
    
    def _build_fasta_index(self):
        
        self.fasta_index = {}
        
        for spec, filename in self.species.items():
            self.fasta_index[spec] = SeqIO.index(filename, "fasta")
            
            for gene_specifier in self.fasta_index[spec].keys():
                self.BMG.add_node("{}_{}".format(spec, gene_specifier),
                                  gene_specifier=gene_specifier,
                                  color=spec)

    
    def _get_seq_from_label(self, label, spec):
        
#        # original label in fasta file (and index) without 'SPECIES_'
#        label = label[len(spec)+1:]
        
        return self.fasta_index[spec][label].seq
    
    
    def _best_matches_for_species(self, spec_x):
        
        candidates = parse_best_match_candidates(self.candidate_files[spec_x])
        
        i = 1
        for x in candidates.keys():
            
            print("----- Processing gene {}/{} -----".format(i, len(candidates)))
            i += 1
            
            x_label = "{}_{}".format(spec_x, x)
            
            for spec_Y, Y in candidates[x].items():
                
                # trivial case: |Y| = 1
                if len(Y) == 1:
                    best_matches = [ Y[0][1] ]
                else:
                    outgroups = self.outgroup_finder(spec_x, spec_Y, candidates[x],
                                                     limit=10)
                    
                    if outgroups:
                        best_matches = self._quartet_method(x, spec_x, Y, spec_Y,
                                                            outgroups)
                    else:
                        best_matches = self._ext_best_hits(x, spec_x, Y, spec_Y)
                
                for bm in best_matches:
                    y_label = "{}_{}".format(spec_Y, bm)
                    self.BMG.add_edge(x_label, y_label)
                    
                    
    def _get_distance(self, a, spec_a, b, spec_b):
        
        # pairwise distances are symmetric --> compute just once
        if (a, spec_a) < (b, spec_b):
            identifier = (a, spec_a, b, spec_b)
        else:
            identifier = (b, spec_b, a, spec_a)
        
        if identifier in self.distance_cache:
            return self.distance_cache[identifier]
        else:
            distance = distance_2seqs(self.fasta_index[spec_a][a].seq,
                                      self.fasta_index[spec_b][b].seq)
#            distance = random.random()
            self.distance_cache[identifier] = distance
            return distance
        
    
    def _ext_best_hits(self, x, spec_x, Y, spec_Y):
        
        best_matches = []
        
        distances, min_distance = [], float("inf")
        for _, y, _, _, _ in Y:
            distance_to_y = self._get_distance(x, spec_x, y, spec_Y)
            distances.append(distance_to_y)
            if distance_to_y < min_distance:
                min_distance = distance_to_y
                
        for i in range(len(Y)):
            if distances[i] < (1 + self.epsilon) * min_distance:
                best_matches.append(Y[i][1])
        
        return best_matches
    
    
    def _quartet_method(self, x, spec_x, Y, spec_Y, Z):
        
        H = nx.DiGraph()
#        f = 1 - self.insecurity_factor
        
        for y1, y2 in itertools.combinations(Y, 2):
            
            y1, y2 = y1[1], y2[1]
            
            votes = [0, 0, 0, 0]
            for spec_z, z, _, _, _ in Z:
                q, weight = self._supported_quartet(x,      spec_x,
                                                    y1, y2, spec_Y,
                                                    z,      spec_z)
                if self.voting_mode == "weighted sum":
                    votes[q] += weight
                elif self.voting_mode == "majority":
                    votes[q] += 1
            
            sorted_votes = sorted((vote, i) for i, vote in enumerate(votes))
            
            if sorted_votes[2][0] > 0.9 * sorted_votes[3][0]:   # not sure if the two 
                quartet = 3                                     # top votes are too close
                                                                # default is star topology
            else:
                quartet = sorted_votes[3][1]
            
            if quartet == 0:                                    # 0: xy1 | y2z
                H.add_edge(y2, y1)
            elif quartet == 1:                                  # 1: xy2 | y1z
                H.add_edge(y1, y2)
            else:                                               # 2: xz | y1y2   or   3: star
                H.add_edge(y1, y2)
                H.add_edge(y2, y1)
        
        # find the (unique) strongly connected component without out-edges
        sccs = {i: scc for i, scc in enumerate(nx.strongly_connected_components(H))}
        
        for i, scc in sccs.items():
            for y in scc:
                H.nodes[y]["scc"] = i
        while sccs:
            for i, scc in sccs.items():
                out_edges = False
                for y1 in scc:
                    for y2 in H.successors(y1):
                        if H.nodes[y1]["scc"] != H.nodes[y2]["scc"]:
                            out_edges = True
                            break
                    if out_edges:
                        break
                if not out_edges:
                    return list(scc)
    
    
    def _supported_quartet(self, x,      spec_x,
                                 y1, y2, spec_y,
                                 z,      spec_z):
        
        xy1_y2z = (self._get_distance(x,  spec_x, y1, spec_y) +
                   self._get_distance(y2, spec_y, z,  spec_z))
        xy2_y1z = (self._get_distance(x,  spec_x, y2, spec_y) +
                   self._get_distance(y1, spec_y, z,  spec_z))
        xz_y1y2 = (self._get_distance(x,  spec_x, z,  spec_z) +
                   self._get_distance(y1, spec_y, y2, spec_y))
        
        if xy1_y2z < xy2_y1z and xy1_y2z < xz_y1y2:     # 0: xy1 | y2z
            quartet = 0
        elif xy2_y1z < xy1_y2z and xy2_y1z < xz_y1y2:   # 1: xy2 | y1z
            quartet = 1
        elif xz_y1y2 < xy1_y2z and xz_y1y2 < xy2_y1z:   # 2: xz | y1y2
            quartet = 2
        else:                                           # 3: star
            quartet = 3
            
        d0, d1, d2 = sorted([xy1_y2z, xy2_y1z, xz_y1y2])
        weight = (1 - d0/d2) * np.exp(d1 - d2) if d2 > 0 else 0.0
        
        return quartet, weight
    
    
def print_BMG(BMG, filename):
    
    with open(filename, "w") as f:
        
        for u, v in BMG.edges:
            
            f.write("{}\t{}\t{}\t{}\n".format(BMG.nodes[u]["gene_specifier"],
                                              BMG.nodes[v]["gene_specifier"],
                                              BMG.nodes[u]["color"],
                                              BMG.nodes[v]["color"],))