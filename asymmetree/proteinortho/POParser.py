# -*- coding: utf-8 -*-

""" 
Parser for ProteinOrtho files.

Reads two different versions of ProteinOrtho result files (with/without
'-synteny' option) and builds a graph from the data (color --> 
species/filename).

References:
    - Marcus Lechner, Sven Findeiß, Lydia Steiner, Manja Marz, Peter F.
      Stadler, Sonja J. Prohaska. Proteinortho: Detection of (Co-)orthologs
      in large-scale analysis. BMC Bioinformatics 2011, 12:124
"""

import time
from multiprocessing import Pool

import networkx as nx

from asymmetree.proteinortho.SpeciesTreeFromParalogs import TreeReconstructor   


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


def parse_po_graph(filename):
    
    G = nx.Graph()
    
    with open(filename, "r") as f:
        
        # skip the first two lines
        f.readline(); f.readline
        
        col_a, col_b= "unknown", "unknown"
        line = f.readline().strip()
        
        while line:
            
            line = line.strip()
            
            # new color pair begins
            if line and line.startswith("#") and not line.startswith("# Scores:"):
                colors = line.split()
                if len(colors) == 3:
                    col_a, col_b = colors[1], colors[2]
                
            elif line and not line.startswith("#"):
                edge = line.split()
                id1 = "{}_{}".format(col_a, edge[0])
                id2 = "{}_{}".format(col_b, edge[1])
                if id1 not in G:
                    G.add_node(id1, color=col_a)
                if id2 not in G:
                    G.add_node(id2, color=col_b)
                    
                # file type 1 with simscore (synteny was activated)
                if len(edge) == 8:
                    G.add_edge(id1, id2,
                               evalue_ab = float(edge[2]), 
                               bitscore_ab = float(edge[3]),
                               evalue_ba = float(edge[4]),
                               bitscore_ba = float(edge[5]),
                               same_strand = int(edge[6]),
                               simscore = float(edge[7]))
                
                # file type 2 without simscore
                elif len(edge) == 6:
                    G.add_edge(id1, id2,
                               evalue_ab = float(edge[2]), 
                               bitscore_ab = float(edge[3]),
                               evalue_ba = float(edge[4]),
                               bitscore_ba = float(edge[5]))
                else:
                    print("Check file format!", edge)
                    continue
                
            line = f.readline()
            
    return G


def parse_best_match_candidates(filename):
    
    candidates = {}
    
#    # extract the species from filename
#    species_match = re.search(r"([^/]+)\.bm_candidates$", filename)
#    if species_match:
#        species = species_match.group(1)
#    else:
#        raise ValueError(f"Could not extract species from filename '{filename}'.")
    
    with open(filename, "r") as f:
        
        line = f.readline()
        
        while line:
            
            line = line.strip()
            if not line:
                line = f.readline()
                continue
            
            data = line.split()
            
            if len(data) != 6:
                print(f"Check file format: {filename}, line {line}")
                continue
                
            color      = data[0]
            query_id   = data[1]
            subject_id = data[2]
            e_value    = float(data[3])
            bitscore   = float(data[4])
            identity   = float(data[5]) / 100
            
            if subject_id not in candidates:
                candidates[subject_id] = {}
            if color not in candidates[subject_id]:
                candidates[subject_id][color] = []
            candidates[subject_id][color].append( (color, query_id, e_value, bitscore, identity) )
                
            line = f.readline()
            
    return candidates


def parse_best_matches(filename):
    
    G = nx.DiGraph()
    
    with open(filename, "r") as f:
        
        line = f.readline()
        
        while line:
            
            line = line.strip().split()
            
            if len(line) == 4:
                
                col1, col2 = line[2], line[3]
                id1 = "{}_{}".format(col1, line[0])
                id2 = "{}_{}".format(col2, line[1])
                if id1 not in G:
                    G.add_node(id1, color=col1)
                if id2 not in G:
                    G.add_node(id2, color=col2)
                G.add_edge(id1, id2)
            
            line = f.readline()
            
    return G


def reconstruct(params):
    
    G, cotree_mode, triple_mode = params
    
    start_time = time.time()
    tr = TreeReconstructor(cotree_mode=cotree_mode)
    tr.add_ortho_graph(G)
    tr.build_species_tree(mode=triple_mode)
    newick = tr.newick_with_support()
    
    return newick, time.time() - start_time


def reconstruct_trees_and_write(infile, outfile, cotree_modes, triple_modes,
                                parallel_processing=True):
    
    G = parse_po_graph(infile)
    
    inputs = [(G, m1, m2) for m1 in cotree_modes for m2 in triple_modes]
    
    if parallel_processing:
        with Pool() as p:
            results = p.map(reconstruct, inputs)
    else:
        results = list(map(reconstruct, inputs))
    
    with open(outfile, "w") as f:
        for i in range(len(inputs)):
            
            _, cotree_mode, triple_mode = inputs[i]
            newick, time_needed = results[i]
            f.write(f"# Cotree usage mode: {cotree_mode}, "\
                    f"Max. Consistent Triple Set heuristic: {triple_mode}, "\
                    f"Time: {time_needed}\n")
            f.write(newick + "\n")


if __name__ == "__main__":
    
    import os
    
    directory = "../../testing/test_files_2"
    
    infile = os.path.join(directory, "test.proteinortho-graph")
    outfile = os.path.join(directory, "species_trees_heuristic")
    
    cotree_modes = ["best", "all"]
    triple_modes = ["Mincut", "BPMF", "Greedy"]
    
    reconstruct_trees_and_write(infile, outfile,
                                cotree_modes, triple_modes,
                                parallel_processing=True)

    with open(outfile, "r") as f:
        print(f.read())
    
#    tr = TreeReconstructor(cotree_mode="all")
#    tr.add_ortho_graph(G)
#    spec_tree = tr.build_species_tree(mode="BPMF")
#    
#    print(spec_tree.to_newick())
#    print(tr.newick_with_support())
#    
#    G = parse_best_matches("test_files_2/bmg")
#    print(G.order(), G.size())