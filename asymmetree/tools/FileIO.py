# -*- coding: utf-8 -*-

"""
File I/O.
"""

import itertools

import networkx as nx

__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


def matrix_to_phylip(filename, leaves, matrix):
    """Write the distance matrix in phylip format."""
    with open(filename, "w") as f:
        f.write(str(len(leaves)))
        for i in range(len(leaves)):
            f.write("\n" + str(leaves[i].ID))
            for j in range(len(leaves)):
                f.write("\t" + str(matrix[i,j]))


def species_to_genes(filename, scenario):
    """Write the species and corresponding genes file."""
    with open(filename, "w") as f:
        species = list(scenario.color_dict)
        for i in range(len(species)):
            f.write(str(species[i]))
            for gene in scenario.color_dict[species[i]]:
                f.write("\t" + str(gene.ID))
            if i < len(species)-1:
                f.write("\n")


def species_pairs_outgroups(filename, scenario):
    """Write pairs of species with corresponding outgroups 
    
    ("ok-species" file)."""
    with open(filename, "w") as f:
        begin = True
        for i in range(len(scenario.subtree_list)):
            for speciesX, speciesY in itertools.combinations(scenario.subtree_list[i], 2):
                if begin:
                    f.write(str(speciesX) + "\t" + str(speciesY))
                    begin = False
                else:
                    f.write("\n" + str(speciesX) + "\t" + str(speciesY))
                for j in range(len(scenario.subtree_list)):
                    if i == j:
                        continue
                    for outgroup_species in scenario.subtree_list[j]:
                        f.write("\t" + str(outgroup_species))
   
                 
def species_in_subtree(filename, scenario):
    """Species in the same subtree under the the root
    
    (one line per subtree)."""
    with open(filename, "w") as f:
        begin = True
        for sp_list in scenario.subtree_list:
            if begin:
                f.write("\t".join([str(item) for item in sp_list]))
                begin = False
            else:
                f.write("\n")
                f.write("\t".join([str(item) for item in sp_list]))


def subtree_BMG(filename, scenario):
    """Write the (subgraph of) BMG."""
    with open(filename, "w") as f:
        begin = True
        for u, v in scenario.BMG_subtrees.edges:
            if begin:
                f.write(str(u) + "\t" + str(v))
                begin = False
            else:
                f.write("\n" + str(u) + "\t" + str(v))


def parse_BMG_edges(output, scenario):
    """Parse edge output from external program to a graph
    
    (converts keys to ints if possible)."""
    BMG = nx.DiGraph()
    BMG.add_nodes_from(scenario.BMG.nodes(data=True))
    
    for line in output.split("\n"):
        if line:
            edge = line.split()
            if len(edge) >= 2:
                u = int(edge[0]) if edge[0].isdigit() else edge[0]
                v = int(edge[1]) if edge[1].isdigit() else edge[1]
                if v not in BMG.nodes:
                    print("Not in BMG", v)
                if u not in BMG.nodes:
                    print("Not in BMG", u)
                BMG.add_edge(u, v)
    return BMG


def write_newick(filename, tree):
    """Write Newick tree."""
    with open(filename, "w") as f:
        f.write(tree.to_newick())


if __name__ == "__main__":
    
    # generate new files in folder "test_data"
    
    import asymmetree.simulator.TreeSimulator as ts
    import asymmetree.simulator.TreeImbalancer as tm
    
    from asymmetree.simulator.Scenario import Scenario
    from asymmetree.tools.PhyloTree import PhyloTree
    
    matrix_file = "test_data/matrix.phylip"
    species_file = "test_data/spec_genes.txt"
    outgroup_file = "test_data/outgroups.txt"
    subtrees_file = "test_data/species_subtrees.txt"
    bmg_file = "test_data/bmg.txt"
    tree_file = "test_data/species_tree.txt"
    
    DLH_rates = (1,1,0)
    
    #S = ts.build_species_tree(10, planted=True)
    S = PhyloTree.parse_newick("(((8:0.08603999468839801,(10:0.06055381385164242,(12:0.02750356935270675,(14:0.0071494825768602215,15:0.0071494825768602215)13:0.02035408677584653)11:0.03305024449893567)9:0.025486180836755596)2:0.34305906001624026,(((16:0.036223587639635554,(18:0.032127988304223636,19:0.032127988304223636)17:0.0040955993354119214)6:0.1114990457717915,7:0.14772263341142705)4:0.03331012904283408,5:0.18103276245426114)3:0.24806629225037713)1:0.5709009452953617)0:0.0")
    S.reconstruct_IDs()
    S.reconstruct_timestamps()
    print("------------- S -------------")
    print(S.to_newick())
    
    TGT = ts.build_gene_tree(S, DLH_rates)
    print("done")
    TGT = tm.imbalance_tree(TGT, S, baseline_rate=1,
                                  lognormal_v=0.2,
                                  gamma_param=(0.5, 1.0, 2.2),
                                  weights=(1/3, 1/3, 1/3),
                                  copy_tree=False)
    
    print("------------- TGT -------------")
    print(TGT.to_newick(distance_only=False))
    
    scenario = Scenario(S, TGT, DLH_rates)
    print(scenario.subtree_list)
    D = scenario.get_distance_matrix()
    
    # write the distance matrix
    matrix_to_phylip(matrix_file, scenario.genes, D)
    # write the species and corresponding genes file
    species_to_genes(species_file, scenario)
    # write the "ok-species" files
    species_pairs_outgroups(outgroup_file, scenario)      
    # write species in root subtrees
    species_in_subtree(subtrees_file, scenario)
    # write the BMG
    subtree_BMG(bmg_file, scenario)
    # write species tree to newick format
    write_newick(tree_file, S)
    
    # test reading input
    with open(bmg_file, "r") as f:
        BMG = parse_BMG_edges(f.read(), scenario)
        #print(BMG.edges)