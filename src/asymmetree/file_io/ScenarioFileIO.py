# -*- coding: utf-8 -*-

"""
File I/O related to module 'Scenario'.
"""

import itertools

import networkx as nx


__author__ = 'David Schaller'


def matrix_to_phylip(filename, leaves, matrix):
    """Write the distance matrix in phylip format."""
    
    with open(filename, 'w') as f:
        f.write(str(len(leaves)))
        for i in range(len(leaves)):
            f.write('\n' + str(leaves[i].ID))
            for j in range(len(leaves)):
                f.write('\t' + str(matrix[i,j]))


def species_to_genes(filename, scenario):
    """Write the species and corresponding genes file."""
    
    with open(filename, 'w') as f:
        species = list(scenario.color_dict)
        for i in range(len(species)):
            f.write(str(species[i]))
            for gene in scenario.color_dict[species[i]]:
                f.write('\t' + str(gene.ID))
            if i < len(species)-1:
                f.write('\n')


def species_pairs_outgroups(filename, scenario):
    """Write pairs of species with corresponding outgroups 
    
    ("ok-species" file)."""
    
    with open(filename, 'w') as f:
        begin = True
        for i in range(len(scenario.subtree_list)):
            for speciesX, speciesY in itertools.combinations(scenario.subtree_list[i], 2):
                if begin:
                    f.write('{}\t{}'.format(speciesX, speciesY))
                    begin = False
                else:
                    f.write('\n{}\t{}'.format(speciesX, speciesY))
                for j in range(len(scenario.subtree_list)):
                    if i == j:
                        continue
                    for outgroup_species in scenario.subtree_list[j]:
                        f.write('\t' + str(outgroup_species))
   
                 
def species_in_subtree(filename, scenario):
    """Species in the same subtree under the the root
    
    (one line per subtree)."""
    
    with open(filename, 'w') as f:
        begin = True
        for sp_list in scenario.subtree_list:
            if begin:
                f.write('\t'.join([str(item) for item in sp_list]))
                begin = False
            else:
                f.write('\n')
                f.write('\t'.join([str(item) for item in sp_list]))


def subtree_bmg(filename, scenario):
    """Write the (subgraph of) BMG."""
    
    with open(filename, 'w') as f:
        begin = True
        for u, v in scenario.bmg_subtrees.edges:
            if begin:
                f.write('{}\t{}'.format(u, v))
                begin = False
            else:
                f.write('\n{}\t{}'.format(u, v))


def parse_bmg_edges(output, scenario):
    """Parse edge output from external program to a graph
    
    (converts keys to ints if possible)."""
    
    bmg = nx.DiGraph()
    bmg.add_nodes_from(scenario.bmg.nodes(data=True))
    
    for line in output.split("\n"):
        if line:
            edge = line.split()
            if len(edge) >= 2:
                u = int(edge[0]) if edge[0].isdigit() else edge[0]
                v = int(edge[1]) if edge[1].isdigit() else edge[1]
                if v not in bmg.nodes:
                    print('Not in BMG', v)
                if u not in bmg.nodes:
                    print('Not in BMG', u)
                bmg.add_edge(u, v)
                
    return bmg


def write_newick(filename, tree):
    """Write Newick tree."""
    
    with open(filename, 'w') as f:
        f.write(tree.to_newick())
