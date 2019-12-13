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

import networkx as nx


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
            
            # new color pair begins
            if line.startswith("#") and not line.startswith("# Scores:"):
                colors = line.split()
                if len(colors) == 3:
                    col_a, col_b = colors[1], colors[2]
                
            elif not line.startswith("#"):
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
                
            line = f.readline().strip()
            
    return G


if __name__ == "__main__":
    
#    G = parse_po_graph("test_files_2/test.proteinortho-graph")
    G = parse_po_graph("test_files_2/test.blast-graph")


    from SpeciesTreeFromParalogs import TreeReconstructor
    
    tr = TreeReconstructor()
    tr.add_ortho_graph(G)
    spec_tree = tr.build_species_tree(mode="BPMF")
    
    print(spec_tree.to_newick())