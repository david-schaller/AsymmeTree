# -*- coding: utf-8 -*-

import os, pickle, random

import pyvolve as pv

import simulator.TreeSimulator as ts
import simulator.TreeImbalancer as tm
from simulator.Tree import Tree


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
        TGT = ts.build_gene_tree(S, (1,0.5,0))                    # only dupl./loss, HGT is disabled
        TGT = tm.imbalance_tree(TGT, S, baseline_rate=1,
                                lognormal_v=0.2,
                                gamma_param=(0.5, 1.0, 2.2),
                                weights=(1, 1, 1),
                                copy_tree=False)
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


for i in range(100):
    
    print("Scenario", i+1, "/ 100")
    
    directory = os.path.join(os.path.expanduser("~"), "Documents", "proteomes", "scenario_{}".format(i))
    
    if os.path.exists(directory) and not os.path.isdir(directory):
        raise NotADirectoryError("Path {} exists and is not a directory!".format(directory))
    elif not os.path.exists(directory):
        os.makedirs(directory)
    
#    S = Tree.parse_newick("(((((16:0.38786287055727103,(18:0.2071277923445058,19:0.2071277923445058)17:0.18073507821276524)12:0.10075553853805931,(14:0.13224182895383052,15:0.13224182895383052)13:0.3563765801414998)4:0.07517286794939665,(6:0.5373882998574596,(8:0.4434182448023457,(10:0.04929450217312242,11:0.04929450217312242)9:0.3941237426292233)7:0.0939700550551139)5:0.02640297718726732)2:0.2512472266526016,3:0.8150385036973286)1:0.18496149630267142)0:0.0;")
#    S.reconstruct_IDs()
#    S.reconstruct_timestamps()
    
    S = ts.build_species_tree(10)
    
    seqs, gene_trees = simulate_gene_families(S, 1000, seq_length=(200,800),
                                              scaling_factor=0.5)
    print(seqs)
    
    write_species_fastas(seqs, directory, pickle_scenario=(S, gene_trees))
    
    with open(os.path.join(directory, "true_species_tree"), "w") as f:
        f.write(to_newick_faa(S))