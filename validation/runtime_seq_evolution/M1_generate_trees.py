# -*- coding: utf-8 -*-

import os

import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick

directory = 'testfiles_trees'
tree_sizes = [10, 20, 50] #, 100, 200]
lengths = [50, 100, 200] #, 500, 1000, 2000]
repeats = 100

if not os.path.exists(directory):
    os.mkdir(directory)


for i in range(repeats):
    
    print(i)
    
    for N in tree_sizes:
        
        S = te.simulate_species_tree(N, model='innovation',
                                     planted=False)
        
        for v in S.preorder():
            v.label = f't{v.label}'
        
        newick = to_newick(S, color=False, distance=True, label_inner=True)
        
        with open(os.path.join(directory, f'tree_{i}_{N}'), 'w') as f:
            f.write(newick)
            f.write('\n')
        
        for l in lengths:
            
            outfile1 = os.path.join('indelible1', f'seqs_{i}_{N}_{l}')
            
            # configuration files for INDELible
            with open(os.path.join(directory, f'tree_{i}_{N}_{l}.indelible1'), 'w') as f:
                f.write( '[TYPE] AMINOACID 1\n'\
                         '[MODEL] m1\n'\
                         '  [submodel] WAG\n')
                f.write(f'[TREE] t1 {newick}\n')
                f.write(f'[PARTITIONS] p1 [t1 m1 {l}]\n')
                f.write(f'[EVOLVE] p1 1 {outfile1}\n')
            
            outfile2 = os.path.join('indelible2', f'seqs_{i}_{N}_{l}')
            
            # configuration files for INDELible
            with open(os.path.join(directory, f'tree_{i}_{N}_{l}.indelible2'), 'w') as f:
                f.write( '[TYPE] AMINOACID 2\n'\
                         '[MODEL] m1\n'\
                         '  [submodel] WAG\n'\
                         '  [rates] 0 1 5\n')
                f.write(f'[TREE] t1 {newick}\n')
                f.write(f'[PARTITIONS] p1 [t1 m1 {l}]\n')
                f.write(f'[EVOLVE] p1 1 {outfile2}\n')
            