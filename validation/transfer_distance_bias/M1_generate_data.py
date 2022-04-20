#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, itertools, multiprocessing, threading

import numpy as np

from tralda.datastructures import LCA

import asymmetree.treeevolve as te
import asymmetree.tools.PhyloTreeTools as ptt


__author__ = 'David Schaller'


# --------------------- parameters ---------------------

result_dir = 'results'
outfile = 'simulation.csv'
outfile2 = 'distances.csv'
# tree_dir = 'trees1000'
threads = os.cpu_count()
timeout_sec = 10000

if not os.path.exists(result_dir):
    os.makedirs(result_dir)
    
# species tree
S_min, S_max = 10, 100   # number of species
contraction_probability = 0.0
    
rates = [(0.25, 0.25, 0.25), (0.5, 0.5, 0.5),
         (0.5, 0.5, 1.0), (0.5, 0.5, 1.5),
         (1.0, 1.0, 0.5), (1.0, 1.0, 1.0),
         (1.5, 1.5, 0.5), (1.5, 1.5, 1.0),]

bias_modes = (False, 'inverse', 'exponential')

repeats = 100

# --------------------- simulation ---------------------

def write_header(filename1, filename2):
    
    items = ['sim_ID', 'repeat', 'bias_mode', 
             'dupl_rate', 'loss_rate', 'hgt_rate',
             'extant', 'D_count', 'L_count', 'H_count',]
    
    with open(filename1, 'w') as f:
        f.write(','.join(items))
    
    items = ['sim_ID', 'bias_mode', 'type', 'distance',]
    
    with open(filename2, 'w') as f:
        f.write(','.join(items))
        

def write_line(filename1, filename2, line, distances_repl, distances_add):
    
    with open(filename1, 'a') as f:
        f.write('\n')
        f.write( ','.join(str(item) for item in line) )
    
    with open(filename2, 'a') as f:
        for d in distances_repl:
            f.write(f'\n{line[0]},{line[2]},r,{d}')
        for d in distances_add:
            f.write(f'\n{line[0]},{line[2]},a,{d}')
    

def simulation(sim_ID, repeat, rate_combination, bias_mode):
    
    dupl_rate, loss_rate, hgt_rate = rate_combination
    
    print(f'sim_ID {sim_ID}, '\
          f'repeat {repeat+1}/{repeats}, '\
          f'bias mode {bias_mode}, '\
          f'dupl. rate: {dupl_rate}, '\
          f'loss rate: {loss_rate}, '\
          f'hgt rate {hgt_rate}')
    
    # fixes issue with numpy.random and multiprocessing
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    
    # species tree
    S = te.simulate_species_tree(np.random.randint(S_min, S_max+1),
                                 model='innovations',
                                 contraction_probability=contraction_probability,)
    # S.serialize(os.path.join(result_dir, tree_dir, f'S_{sim_ID}.pickle'))
    
    # gene tree
    T = te.simulate_dated_gene_tree(S,
                                    dupl_rate=dupl_rate,
                                    loss_rate=loss_rate,
                                    hgt_rate=hgt_rate,
                                    prohibit_extinction='per_family',
                                    replace_prob=0.5,
                                    transfer_distance_bias=bias_mode)
    # T.serialize(os.path.join(result_dir, tree_dir, f'T_{sim_ID}.pickle'))
    
    # S = Tree.load(os.path.join(result_dir, tree_dir,
    #                             f'S_{sim_ID}.pickle'))
    
    # T = Tree.load(os.path.join(result_dir, tree_dir,
    #                            f'T_{sim_ID}.pickle'))
    
    counts_T = ptt.count_node_types(T)
    
    S_label_to_node = {v.label: v for v in S.preorder()}
    T_label_to_node = {v.label: v for v in T.preorder()}
    lca_S = LCA(S)
    lca_T = LCA(T)
    
    distances_repl = []
    distances_add = []
    
    for v in T.preorder():
        if v.event == 'H':
            
            # replacing HGT
            if hasattr(v, 'replaced_gene'):
                v_repl = T_label_to_node[v.replaced_gene]
                
                distances_repl.append( lca_T(v, v_repl).tstamp - v.tstamp )
            
            # additive HGT
            else:
                donor_edge = S_label_to_node[v.color[1]]
                for c in v.children:
                    if c.transferred:
                        color = c.color
                        break
                if isinstance(color, tuple):
                    color = color[1]
                recip_edge = S_label_to_node[color]
                distances_add.append( lca_S(donor_edge, recip_edge).tstamp - v.tstamp )
                
    
    data_line = [sim_ID, repeat, bias_mode,
                 dupl_rate, loss_rate, hgt_rate,
                 counts_T['extant'], counts_T['D'],
                 counts_T['L'], counts_T['H'],
                ]
    
    return data_line, distances_repl, distances_add


def mp_with_timeout(tasks, file_lock, task_lock, outfile):
    
    while True:
        
        with task_lock:
            if tasks:
                task = tasks.pop()
            else:
                task = None
        
        if task is None:
            break
            
        with multiprocessing.Pool(processes=1) as pool:
            async_result = pool.apply_async(simulation, args=(*task,))
            pool.close()
            
            try:
                line, distances_repl, distances_add = async_result.get(timeout=timeout_sec)
                
                with file_lock:
                    write_line(outfile, outfile2,
                               line, distances_repl, distances_add)
                
            except multiprocessing.TimeoutError:
                pool.terminate()
                print(f'{task} took too long ...')
                

# --------------------- main ---------------------
            
if __name__ == '__main__':

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    # if not os.path.exists( os.path.join(result_dir, tree_dir) ):
    #     os.makedirs(os.path.join(result_dir, tree_dir))
    
    outfile = os.path.join(result_dir, outfile)
    outfile2 = os.path.join(result_dir, outfile2)
    write_header(outfile, outfile2)
        
    q = multiprocessing.JoinableQueue()
    file_lock = threading.Lock()
    task_lock = threading.Lock()
    
    tasks = [(i, *item) for i, item in
             enumerate(itertools.product([i for i in range(repeats)],
                                         rates, bias_modes))]
        
    
    tasks = list(reversed(tasks))
    print(tasks)
        
    threads = [threading.Thread(target=mp_with_timeout,
                                args=(tasks, file_lock, task_lock, outfile))
                 for x in range(threads)]
    
    for t in threads:
        t.start()
    
    for t in threads:
        t.join()