# -*- coding: utf-8 -*-

import os, glob, re, itertools

from time import time

# import asymmetree.treeevolve as te
import asymmetree.seqevolve as se
from asymmetree.tools.PhyloTreeTools import parse_newick

import pyvolve

tree_directory = 'testfiles_trees'
seq_directory = 'testfiles_seqs'
result_directory = 'results'
resultfile = os.path.join(result_directory, 'times.csv')

# parameters for among-site heterogeneity
het_alpha = 1.0
het_classes = 5

for direct in (seq_directory, result_directory):
    if not os.path.exists(direct):
        os.mkdir(direct)

files = glob.glob(os.path.join(tree_directory, 'tree*'))
regex = re.compile(r".*tree\_([0-9]+)\_([0-9]+)\_([0-9]+)\.indelible")

tree_sizes = set()
lengths = set()
repeats = set()

for file in files:
    match = regex.match(file)
    if match:
        repeats.add(int(match.group(1)))
        tree_sizes.add(int(match.group(2)))
        lengths.add(int(match.group(3)))

repeats = sorted(repeats)
tree_sizes = sorted(tree_sizes)
lengths = sorted(lengths)
print(repeats)
print(tree_sizes)
print(lengths)

with open(resultfile, 'w') as f:
    f.write('repeat,N,l,tool,mode,time')

for i, N, l in itertools.product(repeats[:2], tree_sizes, lengths):
    print(i, N, l)
    tree_file = os.path.join(tree_directory, f'tree_{i}_{N}')
    
    
    tool = 'asymmetree'
    
    for j in range(1,3):
        outfile = os.path.join(seq_directory, f'{tool}{j}_{i}_{N}_{l}')
        
        start_time = time()
        with open(tree_file, 'r') as f:
            newick = f.readline().strip()
        
        tree = parse_newick(newick)
        
        subst_model = se.SubstModel('a', 'WAG')
        if j == 2:
            het_model = se.HetModel(het_alpha, classes=het_classes)
        else:
            het_model = None
        evolver = se.Evolver(subst_model, het_model=het_model)
        evolver.evolve_along_tree(tree, start_length=l)
        evolver.write_sequences(outfile, include_inner=False)
        end_time = time() - start_time
        
        with open(resultfile, 'a') as f:
            f.write(f'\n{i},{N},{l},{tool},{j},{end_time}')
    
    
    tool = 'pyvolve'
    
    for j in range(1,3):
        outfile = os.path.join(seq_directory, f'{tool}{j}_{i}_{N}_{l}')
        
        start_time = time()
        tree = pyvolve.read_tree(file=tree_file)
        if j == 2:
            aa_model = pyvolve.Model('WAG', alpha=het_alpha, 
                                            num_categories=het_classes)
        else:
            aa_model = pyvolve.Model('WAG')
        partition = pyvolve.Partition(models=aa_model, size=l)
        evolver = pyvolve.Evolver(partitions=partition, tree=tree) 
        evolver(seqfile=outfile)
        
        end_time = time() - start_time
        
        with open(resultfile, 'a') as f:
            f.write(f'\n{i},{N},{l},{tool},{j},{end_time}')