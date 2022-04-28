# -*- coding: utf-8 -*-

import os, decimal

import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick

tree_directory = 'testfiles_trees'
seq_directory = 'testfiles_seqs'

tree_sizes = [10, 20, 50] #, 100, 200]
lengths = [50, 100, 200, 500, 1000, 2000]
repeats = 100

if not os.path.exists(tree_directory):
    os.mkdir(tree_directory)



ctx = decimal.Context()
ctx.prec = 20

def float_to_str(f):
    """Convert the given float to a string, without resorting to scientific 
    notation.
    """
    d1 = ctx.create_decimal(repr(f))
    return format(d1, 'f')


for i in range(repeats):
    
    print(i)
    
    for N in tree_sizes:
        
        S = te.simulate_species_tree(N, model='innovation',
                                     planted=False)
        
        for v in S.preorder():
            v.label = f't{v.label}'
            v.dist = float_to_str(v.dist)
        
        newick = to_newick(S, color=False, distance=True, label_inner=True)
        
        delattr(S.root, 'dist')
        newick2 = to_newick(S, color=False, distance=True, label_inner=False)
        
        with open(os.path.join(tree_directory, f'tree_{i}_{N}'), 'w') as f:
            f.write(newick)
            f.write('\n')
        
        with open(os.path.join(tree_directory, f'tree_{i}_{N}_nolabel'), 'w') as f:
            f.write(newick2)
            f.write('\n')
        
        for l in lengths:
            
            for j in range(1,3):
            
                outfile = os.path.join(seq_directory, f'indelible{j}_{i}_{N}_{l}')
                
                # configuration files for INDELible
                with open(os.path.join(tree_directory,
                                       f'indelible{j}_{i}_{N}_{l}.txt'), 'w') as f:
                    f.write(f'[TYPE] AMINOACID {j}\n'\
                             '[MODEL] m1\n'\
                             '  [submodel] WAG\n')
                    if j == 2:
                        f.write('  [rates] 0 1 5\n')
                    f.write(f'[TREE] t1 {newick2}\n')
                    f.write(f'[PARTITIONS] p1 [t1 m1 {l}]\n')
                    f.write(f'[EVOLVE] p1 1 {outfile}\n')
                
                wdir = os.path.join(seq_directory, f'alf{j}_{i}_{N}_{l}')
                tree_file = os.path.join(tree_directory, f'tree_{i}_{N}_nolabel')
                rate_var_model = 'Gamma, 5, 0, 1' if j == 2 else ''
                
                # configuration files for INDELible
                with open(os.path.join(tree_directory,
                                       f'alf{j}_{i}_{N}_{l}.drw'), 'w') as f:
                    f.write(f"# name of simulation\n"
                             "mname := uuid;\n\n"

                             "# directories for file storage\n"
                            f"wdir := '{wdir}';\n"
                             "dbdir := 'DB/';\n"
                             "dbAncdir := 'DBancestral/';\n\n"

                             "# time scale for simulation (PAM is default)\n"
                             "unitIsPam := false:\n\n"

                             "# parameters concerning the root genome\n"
                             "realseed := false;\n"
                             "protStart := 1;\n"
                            f"minGeneLength := {l};\n"
                             "gammaLengthDist := [1, 1];\n"
                             "blocksize := 1:\n\n"

                             "# parameters concerning the species tree\n"
                             "treeType := 'Custom';\n"
                            f"treeFile := '{tree_file}';\n\n"

                             "# parameters concerning the substitution models\n"
                             "substModels := [SubstitutionModel('WAG')];\n"
                             "indelModels := [IndelModel(0)];\n"
                            f"rateVarModels := [RateVarModel({rate_var_model})];\n"
                             "modelAssignments := [1]:\n"
                             "modelSwitchS := [[1]]:\n"
                             "modelSwitchD := [[1]]:\n\n"
                             
                             "simOutput := {'Fasta', NULL}:\n"
                            )

            