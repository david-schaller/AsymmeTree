# -*- coding: utf-8 -*-

import os, glob

from tralda.datastructures.Tree import Tree, TreeNode

from asymmetree.seqevolve import Evolver, SubstModel
import asymmetree.tools.DistanceCalculation as dc
from asymmetree.file_io.SeqFileIO import write_alignment

from validation import aux_functions

lengths = [100, 500, 1000, 5000, 20000]
distances = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]

subst_models = [SubstModel('n', 'JC69'),
                SubstModel('n', 'K80', kappa=2.0),
                SubstModel('a', 'WAG'),
                SubstModel('a', 'JTT'),]

repeats = (0, 100)

directory = 'testfiles_sequences'

if not os.path.exists(directory):
    os.mkdir(directory)
   

def simulate():
    
    for i in range(*repeats):
        
        for d in distances:
            
            root = TreeNode(label=0, dist=0.0)
            child = TreeNode(label=1, dist=d)
            root.add_child(child)
            T = Tree(root)
            
            for subst_model in subst_models:
                
                evolver = Evolver(subst_model, jump_chain=False)
                
                for l in lengths:
                    
                    outfile = "{}/{}_{}_{}_{}.phylip".format(directory,
                                                             subst_model.model_name,
                                                             d, l, i)
                    print("Processing '{}'".format(outfile))
                    
                    evolver.evolve_along_tree(T, start_length=l)
                    if subst_model.model_name in ('WAG', 'JTT'):
                        
                        # treepuzzle accepts only >= 3 sequences
                        alignment = evolver.true_alignment(include_inner=True)
                        for k, v in alignment.items():
                            seq = v
                            break
                        alignment[TreeNode(label=2)] = seq
                        alignment[TreeNode(label=3)] = seq
                        
                        write_alignment(outfile, alignment,
                                        alignment_format='phylip')
                        
                    else:
                        evolver.true_alignment(include_inner=True,
                                               write_to=outfile)


def calculate_distances(outfile):
    
    table = []
    
    for file in glob.glob(directory + '/*.phylip'):
        
        model, d, l, i = os.path.basename(file).rsplit('.', 1)[0].split('_')
        d, l, i = float(d), int(l), int(i)
        print(model, d, l, i)
        
        if model == 'JC69':
            
            labels, seqs = aux_functions.parse_phylip_alignment(file)
            d_hat = dc.JC69_distance(seqs['0'], seqs['1'])
            
        elif model == 'K80':
            
            labels, seqs = aux_functions.parse_phylip_alignment(file)
            d_hat, kappa = dc.K80_distance(seqs['0'], seqs['1'])
            
        elif model in ('WAG', 'JTT'):
            
            aux_functions.calc_distances(file, model=model)
            D_unsorted, D_labels = aux_functions.parse_distance_matrix(file + '.dist')
            D = aux_functions.rearrange_matrix(D_unsorted, D_labels, ['0', '1', '2', '3'])
            d_hat = D[0, 1]
            os.remove(file + '.dist')
            os.remove(file + '.puzzle')
    
        table.append([model, d, l, i, d_hat])
        
    with open(outfile, 'w') as f:
    
        f.write('model,d,length,i,d_hat')
        
        for row in table:
            f.write('\n' + ','.join(str(value) for value in row))
            
            
        

if __name__ == '__main__':
    
#    simulate()
    
    calculate_distances('results/all_distances.csv')