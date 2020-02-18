# -*- coding: utf-8 -*-

import os, glob

import asymmetree.proteinortho.POParser as parser
    
global_directory = os.path.expanduser("~") + "/Documents/PhD_local_files/ParaphyloVsHeuristic"

scenario_dir = os.path.join(global_directory, "proteomes")
result_dir = os.path.join(global_directory, "species_trees_heuristic")

cotree_modes = ["best", "all"]
triple_modes = ["Mincut", "BPMF", "Greedy"]

if os.path.exists(result_dir) and not os.path.isdir(result_dir):
    raise NotADirectoryError(f"Path '{result_dir}' exists and is not a directory!")
elif not os.path.exists(result_dir):
    os.makedirs(result_dir)
elif os.path.exists(result_dir):
    print(result_dir)
    print(scenario_dir)
    
directories = sorted(glob.glob(scenario_dir + "/*"))
counter = 0
for directory in directories:
    
    counter += 1
    print(f"Processing directory {counter}/{len(directories)}: {directory}")
    
    i = int(directory.rsplit("_", 1)[1])
    infile = os.path.join(directory, "species.proteinortho-graph")
    outfile = os.path.join(result_dir, f"species_tree_{i}")
    
    parser.reconstruct_trees_and_write(infile, outfile,
                                       cotree_modes, triple_modes,
                                       parallel_processing=True)
    
    with open(outfile, "r") as f:
        print(f.read())