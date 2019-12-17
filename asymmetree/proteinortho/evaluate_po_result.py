# -*- coding: utf-8 -*-

import pickle, re
import networkx as nx

from best_match_infer.TrueBMG import true_orthology_graph

from tools.GraphTools import performance


def merge_ortho_graphs(gene_trees):
    
    G = nx.Graph()
    
    for fam_id, T in gene_trees:
        
        ortho_graph = true_orthology_graph(T)
        
        for x in ortho_graph.nodes:
            id_x = "{}.faa_fam{}_gene{}".format(ortho_graph.nodes[x]['color'], fam_id, x)
            G.add_node(id_x)
        
        for x, y in ortho_graph.edges:
            id_x = "{}.faa_fam{}_gene{}".format(ortho_graph.nodes[x]['color'], fam_id, x)
            id_y = "{}.faa_fam{}_gene{}".format(ortho_graph.nodes[y]['color'], fam_id, y)
            
            G.add_edge(id_x, id_y)
    
    return G


def gene_family_subgraphs_from_PO(po_graph):
    
    label = re.compile(r".*_fam([0-9]*)_gene([0-9]+)")
    
    subgraphs = {}
    
    for u, v in po_graph.edges:
        
        label_match_u = label.match(u)
        label_match_v = label.match(v)
        fam_id_u = int(label_match_u.group(1))
        fam_id_v = int(label_match_v.group(1))
        
        if fam_id_u == fam_id_v:
            if fam_id_u not in subgraphs:
                subgraphs[fam_id_u] = nx.Graph()
            subgraphs[fam_id_u].add_edge(u, v)
        else:
            print("Detected match between non-family members {} and {}".format(u,v))
    
    return subgraphs


def gene_family_subgraphs_from_trees(gene_trees):
    
    subgraphs = {}
    
    for fam_id, T in gene_trees:
        
        G = nx.Graph()
        subgraphs[fam_id] = G
        
        ortho_graph = true_orthology_graph(T)
        
        for x in ortho_graph.nodes:
            id_x = "{}.faa_fam{}_gene{}".format(ortho_graph.nodes[x]['color'], fam_id, x)
            G.add_node(id_x)
        
        for x, y in ortho_graph.edges:
            id_x = "{}.faa_fam{}_gene{}".format(ortho_graph.nodes[x]['color'], fam_id, x)
            id_y = "{}.faa_fam{}_gene{}".format(ortho_graph.nodes[y]['color'], fam_id, y)
            
            G.add_edge(id_x, id_y)
        
    
    return subgraphs
        


from POParser import parse_po_graph

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
    
    subgraphs_PO = gene_family_subgraphs_from_PO(po_graph)
    subgraphs_true = gene_family_subgraphs_from_trees(gene_trees)
    print(len(subgraphs_PO), len(subgraphs_true))
    
    sum1, sum2 = 0, 0
    accuracy, precision, recall = [], [], []
    for i in subgraphs_PO.keys():
        order1 = subgraphs_true[i].order()
        order2 = subgraphs_PO[i].order()
        sum1 += order1
        sum2 += order2
        _, _, _, _, _, _, acc_i, prec_i, recall_i = performance(subgraphs_true[i], subgraphs_PO[i])
        accuracy.append(acc_i)
        precision.append(prec_i)
        recall.append(recall_i)
        print("{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}".format(order1, order2, acc_i, prec_i, recall_i))
    print(sum1, sum2)
    
    import matplotlib.pyplot as plt
    plt.boxplot([precision, recall])