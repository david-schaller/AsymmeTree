# -*- coding: utf-8 -*-

"""
This simulation script accompanies the following publication:

D. Schaller, M. Hellmuth, P. F. Stadler (2023)
Orientation of Fitch Graphs and Reconciliation-Free Inference of Horizontal 
Gene Transfer in Gene Trees.
To appear in SIAM Journal on Discrete Mathematics.
https://arxiv.org/abs/2112.00403
"""


import os, itertools
from concurrent import futures

import numpy as np

from tralda.datastructures import Tree, LCA
from tralda.cograph.Cograph import to_cotree
import tralda.tools.GraphTools as gt

import asymmetree.treeevolve as te
import asymmetree.analysis.BestMatches as bm
import asymmetree.analysis.HGT as hgt
import asymmetree.tools.PhyloTreeTools as ptt


__author__ = 'David Schaller'


# --------------------- parameters ---------------------

result_dir = 'fitch_graph_orientation_results'
outfile = 'simulation.csv'
tree_dir = 'trees'
parallel_processes = os.cpu_count()

# species tree
S_min, S_max = 10, 100   # range of number of species

# frequency of duplication/loss/HGT
rates = [(0.25, 0.25, 0.25), (0.5, 0.5, 0.5),
         (0.5, 0.5, 1.0), (0.5, 0.5, 1.5),
         (1.0, 1.0, 0.5), (1.0, 1.0, 1.0),
         (1.5, 1.5, 0.5), (1.5, 1.5, 1.0),]

# number of simulatios per combination of rates
repeats = 5000


# --------------------- simulation ---------------------

def write_header(filename):
    
    items = ['sim_ID', 'repeat', 'dupl_rate', 'loss_rate', 'hgt_rate',
             'extant', 'D_count', 'L_count', 'H_count',
             'no_sets', 'dfitch_order', 'dfitch_size', 'ufitch_size',
             'HGT_ess', 'HGT_forb', 'HGT_ambp', 'HGT_amba',]
    
    for tree in ('full', 'cotree', 'lrt', 'brt'):
        items.append('{}_res'.format(tree))
        for comp in ('comp', 'rcomp'):
            for stat in ('pairs', 'ess', 'forb', 'ambp', 'amba', 'wrong'):
                for mode in ('', '_tot'):
                    items.append('{}_{}_{}{}'.format(tree, comp, stat, mode))
    
    with open(filename, 'w') as f:
        f.write(','.join(items))
        

def write_line(filename, line):
    
    with open(filename, 'a') as f:
        f.write('\n')
        f.write( ','.join(str(item) for item in line) )
        
        
def edge_classification(T, P):
    
    vertex_coloring, _ = hgt.is_compatible(T, P)
    if not vertex_coloring:
        return [False, False, False, False]
    
    essential, forbidden, ambiguous_p, ambiguous_a = 0, 0, 0, 0
    
    for u, v in T.edges():
        
        if u in vertex_coloring and v in vertex_coloring:
            if vertex_coloring[u] == vertex_coloring[v]:
                forbidden += 1
            else:
                essential += 1
        elif v.transferred:
            ambiguous_p += 1
        else:
            ambiguous_a += 1
    
    return essential, forbidden, ambiguous_p, ambiguous_a
        

def fitch_stats(P, dfitch, matrix):
    
    pairs, pairs_tot =           0, 0
    essential, essential_tot =   0, 0
    forbidden, forbidden_tot =   0, 0
    ambiguousp, ambiguousp_tot = 0, 0
    ambiguousa, ambiguousa_tot = 0, 0
    wrong, wrong_tot =           0, 0
    
    for i in range(len(P)):
        for j in range(len(P)):
            if i == j:
                continue
            
            # representatives
            x = next(iter(P[i]))
            y = next(iter(P[j]))
            
            # total number of pairs
            tot_pairs = len(P[i]) * len(P[j])
            pairs += 1
            pairs_tot += tot_pairs
            
            if matrix[i][j] == 'essential' and dfitch.has_edge(x, y):
                essential += 1
                essential_tot += tot_pairs
            elif matrix[i][j] == 'forbidden' and not dfitch.has_edge(x, y):
                forbidden += 1
                forbidden_tot += tot_pairs
            elif matrix[i][j] == 'ambiguous' and dfitch.has_edge(x, y):
                ambiguousp += 1
                ambiguousp_tot += tot_pairs
            elif matrix[i][j] == 'ambiguous' and not dfitch.has_edge(x, y):
                ambiguousa += 1
                ambiguousa_tot += tot_pairs
            else:
                wrong += 1
                wrong_tot += tot_pairs
    
    return (pairs, pairs_tot,
            essential, essential_tot,
            forbidden, forbidden_tot,
            ambiguousp, ambiguousp_tot,
            ambiguousa, ambiguousa_tot,
            wrong, wrong_tot)
    

def simulation(arg):
    
    sim_ID, repeat, (dupl_rate, loss_rate, hgt_rate) = arg
    
    print(f'Started task sim_ID {sim_ID}, '\
          f'repeat {repeat+1}/{repeats}, '\
          f'dupl. rate: {dupl_rate}, '\
          f'loss rate: {loss_rate}, '\
          f'hgt rate {hgt_rate}')
    
    # fixes issue with numpy.random and multiprocessing
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    
    # species tree
    S = te.species_tree_n_age(np.random.randint(S_min, S_max+1), 1.0,
                              model='yule',
                              innovation=True)
    
    # this step is only necessary to ensure consistency with the results in 
    # the paper since they were generated using an older version of 
    # AsymmeTree
    ptt.random_ultrametric_timing(S, inplace=True, adjust_distances=True)
    
    S.serialize(os.path.join(result_dir, tree_dir, f'S_{sim_ID}.pickle'))
    
    # true gene tree
    TGT = te.dated_gene_tree(S,
                             dupl_rate=dupl_rate,
                             loss_rate=loss_rate,
                             hgt_rate=hgt_rate,
                             prohibit_extinction='per_family')
    TGT.serialize(os.path.join(result_dir, tree_dir, f'TGT_{sim_ID}.pickle'))
    
    S = Tree.load(os.path.join(result_dir, tree_dir, f'S_{sim_ID}.pickle'))
    
    TGT = Tree.load(os.path.join(result_dir, tree_dir, f'TGT_{sim_ID}.pickle'))
    
    counts_TGT = ptt.count_node_types(TGT)      # use loss count in TGT
    
    # observable gene tree
    OGT = te.prune_losses(TGT)
    counts_OGT = ptt.count_node_types(OGT)        # use dupl./HGT counts in OGT
    
    dfitch, ufitch = hgt.fitch(OGT, hgt.true_transfer_edges(OGT),
                               supply_undirected=True)
    P = gt.independent_sets(ufitch)
    
    cotree = to_cotree(bm.orthology_from_tree(OGT))
    lrt = bm.lrt_from_tree(OGT)
    brt = bm.binary_refinable_tree(bm.bmg_from_tree(OGT))
    
    data_line = [sim_ID, repeat, dupl_rate, loss_rate, hgt_rate,
                 counts_TGT['extant'], counts_OGT['D'],
                 counts_TGT['L'], counts_OGT['H'],
                 len(P), dfitch.order(), dfitch.size(), ufitch.size(),
                 
                 # HGT edge classification for full tree
                 *edge_classification(OGT, P),
                ]
    
    for T in (OGT, cotree, lrt, brt):
        
        data_line.append(sum(1 for _ in T.inner_nodes()) - 1)
        
        lca_T = LCA(T)
        
        for matrix in [hgt.fitch_orientation(T, P, lca=lca_T),
                       hgt.fitch_orientation_for_refinements(T, P, lca=lca_T)]:
            if matrix:
                data_line.extend(fitch_stats(P, dfitch, matrix))
            else:
                data_line.extend(12 * [False])
    
    print(f'Finished task sim_ID {sim_ID}, '\
          f'repeat {repeat+1}/{repeats}, '\
          f'dupl. rate: {dupl_rate}, '\
          f'loss rate: {loss_rate}, '\
          f'hgt rate {hgt_rate}')
    
    return data_line
                

# --------------------- main ---------------------
            
if __name__ == '__main__':

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    if not os.path.exists( os.path.join(result_dir, tree_dir) ):
        os.makedirs(os.path.join(result_dir, tree_dir))
    
    outfile = os.path.join(result_dir, outfile)
    write_header(outfile)
    
    tasks = [(i, *item) for i, item in
             enumerate(itertools.product(range(repeats), rates))]
    
    with futures.ProcessPoolExecutor(max_workers=parallel_processes) as e:
        
        for result in e.map(simulation, tasks):
            write_line(outfile, result)
