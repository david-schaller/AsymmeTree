# -*- coding: utf-8 -*-

import os, pickle, random

import pyvolve as pv

import asymmetree.simulator.TreeSimulator as ts
import asymmetree.simulator.TreeImbalancer as tm


def to_newick_pyvolve(tree, node=None):
    
    if node is None:
        return to_newick_pyvolve(tree, node=tree.root) + ";"
    elif not node.children:
        return "{}_{}:{}".format(node.color, node.label, node.dist)
    else:
        s = ''
        for child in node.children:
            s += to_newick_pyvolve(tree, node=child) + ","
        return "({}):{}".format(s[:-1], node.dist)
    
    
def to_newick_faa(tree, node=None):
    
    if node is None:
        return to_newick_faa(tree, node=tree.root) + ";"
    elif not node.children:
        return "{}.faa".format(node.label)
    else:
        s = ''
        for child in node.children:
            s += to_newick_faa(tree, node=child) + ","
        return "({})".format(s[:-1])


def simulate_gene_families(S, N, seq_length=(200,800),
                           scaling_factor=1.0):
    
    seq_dict = {}
    gene_trees = []
    
    for v in S.preorder():
        if not v.children:
            seq_dict[str(v.label)] = []
            
    for i in range(N):                         
        TGT_simulator = ts.GeneTreeSimulator(S)
        TGT = TGT_simulator.simulate((1.0, 0.5, 0.0))       # only dupl./loss, HGT is disabled
        TGT = tm.imbalance_tree(TGT, S, baseline_rate=1,
                                lognormal_v=0.2,
                                gamma_param=(0.5, 1.0, 2.2),
                                weights=(1, 1, 1))
        OGT = ts.observable_tree(TGT)
        
        for v in OGT.preorder():
            v.dist *= scaling_factor
        
        phylogeny = pv.read_tree(tree=to_newick_pyvolve(OGT))
        model = pv.Model("JTT")
        
        seq_size = random.randint(*seq_length) if isinstance(seq_length, tuple) else seq_length
        partition = pv.Partition(models=model, size=seq_size)
        evolver = pv.Evolver(tree=phylogeny, partitions=partition)
        evolver(seqfile=False, ratefile=False, infofile=False)
        for label, seq in evolver.get_sequences().items():
            species, gene_label = label.split("_", 1)
            seq_dict[species].append( (i, gene_label, seq) )
        
        gene_trees.append( (i, OGT) )
    
    return seq_dict, gene_trees


def write_species_fastas(seq_dict, dirname, pickle_scenario=None):
        
    for species, seq_list in seq_dict.items():
        
        filename = os.path.join(dirname, "{}.faa".format(species))
        with open(filename, "w") as f:
            for fam_id, gene_label, seq in seq_list:
                f.write(">fam{}_gene{}\n".format(fam_id, gene_label))
                pos = 0
                while pos < len(seq):
                    f.write(seq[pos:min(pos+80, len(seq))])
                    pos += 80
                    f.write("\n")
    
    if pickle_scenario:
        pickle_file = os.path.join(dirname, "scenario.pickle")
        pickle.dump(pickle_scenario, open(pickle_file , "wb" ))