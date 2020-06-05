# -*- coding: utf-8 -*-

import time
from multiprocessing import Pool

from asymmetree.paraphylo.SpeciesTreeFromParalogs import TreeReconstructor
from asymmetree.file_io.ProteinOrthoParser import parse_po_graph


__author__ = 'David Schaller'


def _reconstruct(params):
    
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
            results = p.map(_reconstruct, inputs)
    else:
        results = list(map(_reconstruct, inputs))
    
    with open(outfile, 'w') as f:
        for i in range(len(inputs)):
            
            _, cotree_mode, triple_mode = inputs[i]
            newick, time_needed = results[i]
            f.write('# Cotree usage mode: {}, Max. Consistent Triple Set '\
                    'heuristic: {}, Time: {}\n'.format(cotree_mode,
                                                       triple_mode,
                                                       time_needed))
            f.write(newick + '\n')
            
            
def reconstruct_from_proteinortho(filename, triple_mode='BPMF'):
    
    G = parse_po_graph(filename)
    tr = TreeReconstructor()
    tr.add_ortho_graph(G)
    S = tr.build_species_tree(mode=triple_mode)
    support_newick = tr.newick_with_support()
    
    return S, support_newick