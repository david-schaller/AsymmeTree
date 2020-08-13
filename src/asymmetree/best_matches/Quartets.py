# -*- coding: utf-8 -*-

"""
Quartet approach for best match inference.
"""

import itertools, random,  subprocess, time, os

import numpy as np
import networkx as nx

from asymmetree.tools.GraphTools import symmetric_part
from asymmetree.file_io.ScenarioFileIO import parse_bmg_edges, matrix_to_phylip, species_to_genes, write_newick


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                           PYTHON IMPLEMENTATION
#
# --------------------------------------------------------------------------
# Outgroup selection methods:
#  (1) outgroups from subtrees under the root of S
#  (2) use all possible outgoups for a pair of colors
#      (corrected closest outgroups)
# --------------------------------------------------------------------------

class Quartets:
    """Implementation of the Quartet approach for best match inference."""
    
    def __init__(self, scenario, distance_matrix,
                 voting_mode='majority',
                 use_distant_genes=False,
                 insecurity_factor=0.0,
                 y_candidates=None,
                 closest_outgroups=False,
                 incongruence_threshold=0.2):
        """Inititialize Quartet method.
        
        Keyword arguments:
            voting_mode -- 'majority' or 'weighted' voting for quartet
                           support.
            use_distant_genes -- use the most distant outgroups instead
                                 of a random choice (only root outgroups)
            insecurity_factor -- threshold for inference of the star
                                 topology
            y_candidates -- restrict the Y candidates (dict x: Y)
            closest_outgroups -- apply the corrected closest outgroups
                                 method
            incongruence_threshold -- threshold for discarding
                                      outgroup species
        """
        
        self.scenario = scenario
        
        (self.S, self.genes, self.gene_index,
         self.subtree_list, self.subtree_index,
         self.outgroup_dict, self.color_dict) = scenario.get_data()
        
        self.D = distance_matrix
        self.voting_mode = voting_mode                          # majority voting or weighted sum
        self.use_distant_genes = use_distant_genes              # use most distant outgroups instead of random choice
        self.insecurity_factor = insecurity_factor              # relative distance between top votes
                                                                # to decide for star topology instead
        
        self.y_candidates = y_candidates                        # Graph or DiGraph that contains
                                                                # best match candidates
        
        self.closest_outgroups = closest_outgroups              # use the closest corrected outgroups method
        self.incongruence_threshold = incongruence_threshold    # threshold ancient dupl. correction
    
    
    def compute_lcas_in_S(self):
        
        self.S.supply_leaves()
        N = len(self.S.root.leaves)
        self.S_leaf_index = {leaf.ID: i for i, leaf in enumerate(self.S.root.leaves)}
        
        self.lca_S = np.zeros((N,N), dtype=int)
        self.subtree_species = {}
        self.subtree_genes = {}
        
        for v in self.S.preorder():
            self.subtree_species[v.ID] = [leaf.ID for leaf in v.leaves]
            self.subtree_genes[v.ID] = []
            for leaf_id in self.subtree_species[v.ID]:
                self.subtree_genes[v.ID].extend(self.color_dict[leaf_id])
            
            if not v.children:
                self.lca_S[self.S_leaf_index[v.ID],self.S_leaf_index[v.ID]] = v.ID
            elif len(v.children) >= 2:
                for c1, c2 in itertools.combinations(v.children, 2):
                    for x in c1.leaves:
                        x_index = self.S_leaf_index[x.ID]
                        for y in c2.leaves:
                            y_index = self.S_leaf_index[y.ID]
                            self.lca_S[x_index, y_index] = v.ID
                            self.lca_S[y_index, x_index] = v.ID
    
    
    def compute_outgroup_species(self):
        
        self.lca_to_outgroups = {key: set() for key in self.subtree_species}
        
        root_done = False
        for v in self.S.preorder():
            
            if v.parent:
                self.lca_to_outgroups[v.ID].update(self.lca_to_outgroups[v.parent.ID])
                
            for c1, c2 in itertools.permutations(v.children, 2):
                
                if self.incongruence_threshold >= 0 and root_done:
                    outgroup_candidates = set(self.subtree_species[c2.ID])
                    genes_c1 = self.subtree_genes[c1.ID]
                    species_c2 = self.subtree_species[c2.ID]
                    if len(genes_c1) > 1 and len(species_c2) > 1:
                        for s1, s2 in itertools.combinations(species_c2, 2):
                            genes_c2 = self.color_dict[s1] + self.color_dict[s2]
                            
                            votes = [0,0]                                           # [congruent, incongruent]
                            for i in range(20):
                                a, b = random.sample(genes_c1, 2)
                                c, d = random.sample(genes_c2, 2)
                                quartet, _ = self.supported_quartet(a, b, c, d)
                                if quartet == 0:                                    # ab | cd
                                    votes[0] += 1                                   # (congruent)
                                else:                                               # ac | bd or ad | bc or star
                                    votes[1] += 1                                   # (incongruent)
                                    
                            if votes[1] / sum(votes) >= self.incongruence_threshold:
                                outgroup_candidates.discard(s1)
                                outgroup_candidates.discard(s2)
                            
                    self.lca_to_outgroups[c1.ID].update(outgroup_candidates)
                    
                else:                                               # outgroups on the root are fix
                    self.lca_to_outgroups[c1.ID].update(self.subtree_species[c2.ID])
                root_done = True

    
    def build_matrix_I(self):
        """Build an index matrix that stores for each gene (row) the indeces
        of the other genes in increasing order by distance."""
        
        N = self.D.shape[0]
        self.I = np.tile(np.arange(N), (N,1))
        D, I = self.D, self.I
        
        def quicksort_row(row, l, r):
            if r <= l:
                return
            i, j = l, r
            piv = D[row, I[row, (l+r)//2]]
            while i <= j:
                while D[row, I[row, i]] < piv:
                    i += 1
                while D[row, I[row, j]] > piv:
                    j -= 1
                if i <= j:
                    I[row,i], I[row,j] = I[row,j], I[row,i]
                    i += 1
                    j -= 1
            quicksort_row(row, l, j)
            quicksort_row(row, i, r)
        
        for row in range(N):
            quicksort_row(row, 0, N-1)
        
    
    def get_outgroups_root_only(self, x, limit=10):
        
        i, x_subtree = self.gene_index[x], self.subtree_index[x.color]
        limit = min(limit, len(self.outgroup_dict[x_subtree]))
        
        if self.use_distant_genes:
            outgroups = []
            j = len(self.genes)-1
            while len(outgroups) < limit and j >= 0:
                gene = self.genes[self.I[i,j]]
                if x_subtree != self.subtree_index[gene.color]:
                    outgroups.append(gene)
                j -= 1
        else:
            outgroups = random.sample(self.outgroup_dict[x_subtree], limit)
        
        return outgroups
    
    
    def get_outgroups_closest(self, x, colorY, limit=10):
        
        outgroups_S = self.lca_to_outgroups[self.lca_S[self.S_leaf_index[x.color],
                                                       self.S_leaf_index[colorY]
                                                       ]]
        
        outgroups = []
        i, j, N = self.gene_index[x], 0, self.D.shape[0]
        while len(outgroups) < limit and j < N:
            if self.genes[self.I[i,j]].color in outgroups_S:
                outgroups.append(self.genes[self.I[i,j]])
            j += 1
        
        return outgroups
    
    
    def supported_quartet(self, x, y1, y2, z):
        
        x, y1, y2, z = (self.gene_index[x], self.gene_index[y1],
                        self.gene_index[y2], self.gene_index[z])
        
        xy1_y2z = self.D[x, y1] + self.D[y2, z ]
        xy2_y1z = self.D[x, y2] + self.D[y1, z ]
        xz_y1y2 = self.D[x, z ] + self.D[y1, y2]
        
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
    
    
    def find_best_matches(self, x, Y, Z):
        
        H = nx.DiGraph()
        f = 1 - self.insecurity_factor
        
        for y1, y2 in itertools.combinations(Y, 2):
            votes = [0, 0, 0, 0]
            for z in Z:
                q, weight = self.supported_quartet(x, y1, y2, z)
                if self.voting_mode == "weighted sum":
                    votes[q] += weight                          # use the weight
                else:
                    votes[q] += 1                               # majority vote
            
            sorted_votes = sorted((vote, i) for i, vote in enumerate(votes))
            
            if sorted_votes[2][0] > f * sorted_votes[3][0]:     # not sure if the two 
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
    
    
    def build_graphs_root_only(self):
        
        if self.use_distant_genes:
            self.build_matrix_I()  
        
        self.bmg = nx.DiGraph()
        for g in self.genes:
            self.bmg.add_node(g.ID, label=g.label, color=g.color)
        
        # consider all genes of other color as potential best matches
        if self.y_candidates is None:
            
            for x in self.genes:
                Z = self.get_outgroups_root_only(x)
                
                for species in self.subtree_list[self.subtree_index[x.color]]:
                    Y = self.color_dict[species]
                    if species == x.color:
                        continue
                    elif len(Y) == 1:
                        self.bmg.add_edge(x.ID, Y[0].ID)
                    elif len(Y) > 1:
                        best_matches = self.find_best_matches(x, Y, Z)
                        for y in best_matches:
                            self.bmg.add_edge(x.ID, y.ID)
                            
        # consider only a given subset as best matches
        else:
            id_dict = {g.ID: g for g in self.genes}
            for x in self.genes:
                Z = self.get_outgroups_root_only(x)
                Y = sorted([id_dict[gene_id] for gene_id in self.y_candidates.neighbors(x.ID)],
                           key=lambda y: y.color)               # sort by color
                
                i = 0
                while i < len(Y):                               # handle colors separately
                    j = i
                    while j < len(Y) and Y[i].color == Y[j].color:
                        j += 1
                    if self.subtree_index[x.color] == self.subtree_index[Y[i].color]:
                        if i == j-1:
                            self.bmg.add_edge(x.ID, Y[i].ID)
                        else:
                            best_matches = self.find_best_matches(x, Y[i:j], Z)
                            for y in best_matches:
                                self.bmg.add_edge(x.ID, y.ID)
                    i = j
        
        self.rbmg = symmetric_part(self.bmg)
        return self.bmg, self.rbmg
    
    
    def build_graphs_closest_outgroups(self):
        
        self.compute_lcas_in_S()
        self.compute_outgroup_species()
        self.build_matrix_I()
        
        self.bmg = nx.DiGraph()
        for g in self.genes:
            self.bmg.add_node(g.ID, label=g.label, color=g.color)
        
        # consider all genes of other color as potential best matches
        if self.y_candidates is None:
            
            for x in self.genes:
                for species in self.subtree_list[self.subtree_index[x.color]]:
                    Y = self.color_dict[species]
                    if species == x.color:
                        continue
                    elif len(Y) == 1:
                        self.bmg.add_edge(x.ID, Y[0].ID)
                    elif len(Y) > 1:
                        Z = self.get_outgroups_closest(x, Y[0].color)
                        best_matches = self.find_best_matches(x, Y, Z)
                        for y in best_matches:
                            self.bmg.add_edge(x.ID, y.ID)
                            
        # consider only a given subset as best matches
        else:
            id_dict = {g.ID: g for g in self.genes}
            for x in self.genes:
                Y = sorted([id_dict[gene_id] for gene_id in self.y_candidates.neighbors(x.ID)],
                           key=lambda y: y.color)               # sort by color
                
                i = 0
                while i < len(Y):                               # handle colors separately
                    j = i
                    while j < len(Y) and Y[i].color == Y[j].color:
                        j += 1
                    if self.subtree_index[x.color] == self.subtree_index[Y[i].color]:
                        if i == j-1:
                            self.bmg.add_edge(x.ID, Y[i].ID)
                        else:
                            Z = self.get_outgroups_closest(x, Y[i].color)
                            best_matches = self.find_best_matches(x, Y[i:j], Z)
                            for y in best_matches:
                                self.bmg.add_edge(x.ID, y.ID)
                    i = j
                    
        self.rbmg = symmetric_part(self.bmg)
        return self.bmg, self.rbmg
    
    
    def build_graphs(self):
        
        if not self.closest_outgroups:
            self.build_graphs_root_only()
        else:
            self.build_graphs_closest_outgroups()
            


# --------------------------------------------------------------------------
#                            EXTERNAL C++ PROGRAM
#
# --------------------------------------------------------------------------

def quartet_qinfer(scenario,
                   matrix_filename, species_filename, tree_filename,
                   epsilon=-1,
                   closest_outgroups=False,
                   incongruence_threshold=0.2,
                   benchmark_file=None,
                   binary_path=None):
    """Compute BMG and RBMG from a distances matrix D using 'qinfer'.
    
    Keyword arguments:
        epsilon -- epsilon for relative BM candidate threshold:
                   y in Y if D(x,y) <= (1+epsilon) * min d(x,y'),
                   default=10E-8 (for limited float precision).
        closest_outgroups -- use the closest outgroups corrected
                             for ancient duplications
        incongruence_threshold -- threshold for discarding
                                  outgroup species
        benchmark_file -- activate benchmarking in 'qinfer' and
                          specify the filename
        binary_path -- path to 'qinfer' binary (if not available
                       within path)
    """
    
    if not binary_path:
        qinfer_command = "qinfer"
    elif os.path.exists(binary_path):
        qinfer_command = binary_path
    else:
        raise FileNotFoundError("path to qinfer binary file '{}' does not exist".format(binary_path))
    
    output = -1
    
    command = [qinfer_command, matrix_filename, species_filename, tree_filename]
    
    if epsilon != -1:
        command.append( "--epsilon=" + str(epsilon) )
    if closest_outgroups:
        command.append( "--all-outgroups=" + str(incongruence_threshold) )
    if benchmark_file is not None:
        command.append( "--benchmark=" + benchmark_file )
    
    # call 'qinfer' and measure execution time
    start = time.time()
    
    try:
        output = subprocess.run(command, stdout=subprocess.PIPE)
    except:
        raise FileNotFoundError("calling qinfer failed")
        
    exec_time = time.time() - start
    
    if output == -1:
        raise Exception("no output from qinfer")
    
    bmg = parse_bmg_edges(output.stdout.decode(), scenario)
    rbmg = symmetric_part(bmg)
    
    return bmg, rbmg, exec_time


def quartet_from_scenario(scenario, epsilon=-1,
                          closest_outgroups=False,
                          incongruence_threshold=0.2,
                          benchmark_file=None):
    """Compute BMG and RBMG from a distances matrix D using 'qinfer'.
    
    Keyword arguments:
        epsilon -- epsilon for relative BM candidate threshold:
                   y in Y if D(x,y) <= (1+epsilon) * min d(x,y'),
                   default=10E-8 (for limited float precision).
        closest_outgroups -- use the closest outgroups corrected
                             for ancient duplications
        incongruence_threshold -- threshold for discarding
                                  outgroup species
        benchmark_file -- activate benchmarking in 'qinfer' and
                          specify the filename
    """

    matrix_filename = "temp.phylip"
    species_filename = "temp_species.txt"
    tree_filename = "temp_tree.txt"
    
    matrix = scenario.get_distance_matrix()
    matrix_to_phylip(matrix_filename, scenario.genes, matrix)
    species_to_genes(species_filename, scenario)
    write_newick(tree_filename, scenario.S)
    
    bmg, rbmg, exec_time = quartet_qinfer(scenario,
                                          matrix_filename, species_filename, tree_filename,
                                          epsilon=epsilon,
                                          closest_outgroups=closest_outgroups,
                                          incongruence_threshold=incongruence_threshold,
                                          benchmark_file=benchmark_file)
    
    os.remove(matrix_filename)
    os.remove(species_filename)
    os.remove(tree_filename)
    
    return bmg, rbmg, exec_time