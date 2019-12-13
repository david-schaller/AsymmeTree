# -*- coding: utf-8 -*-

import pickle
import networkx as nx

from best_match_infer.TrueBMG import true_orthology_graph

from tools.GraphTools import performance


def merge_ortho_graphs(gene_trees):
    
    G = nx.Graph()
    
    for fam_id, T in gene_trees:
        
        ortho_graph = true_orthology_graph(T)
        
        for x in ortho_graph.nodes:
            id_x = "fam{}_gene{}".format(fam_id, x)
            G.add_node(id_x)
        
        for x, y in ortho_graph.edges:
            id_x = "{}.faa_fam{}_gene{}".format(ortho_graph.nodes[x]['color'], fam_id, x)
            id_y = "{}.faa_fam{}_gene{}".format(ortho_graph.nodes[y]['color'], fam_id, y)
            
            G.add_edge(id_x, id_y)
    
    return G


from POGraphParser import parse_po_graph

if __name__ == "__main__":
    
    folder = "test_files_2/"
    
    S, gene_trees = pickle.load(open(folder + "scenario.pickle", "rb"))
    print(S.to_newick())
    
    po_graph = parse_po_graph(folder + "test.proteinortho-graph")
    true_ortho_graph = merge_ortho_graphs(gene_trees)
    
    print(po_graph.order(), true_ortho_graph.order())
    
    _, _, tp, tn, fp, fn, accuracy, precision, recall = performance(true_ortho_graph, po_graph)
    print(accuracy, precision, recall)
    
#    print(po_graph.nodes)