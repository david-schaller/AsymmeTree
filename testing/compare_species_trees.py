# -*- coding: utf-8 -*-

import os, glob, re

def parse_paraphylo_output(filename):
    
    build_tree, lrt_tree, time = "", "", 0.0
    
    with open(filename, "r") as f:
        
        line = f.readline()
        while line:
            if line.startswith("#BUILD tree"):
                build_tree = f.readline().strip()
            elif line.startswith("#least resolved tree"):
                lrt_tree = f.readline().strip()
            elif line.startswith("#total time"):
                time = float(f.readline().strip()) / 1000.0
            line = f.readline()
    
    return build_tree, lrt_tree, time


def parse_heuristic(filename):
    
    trees = []
    
    info_line = re.compile(r"# Cotree usage mode: ([A-Za-z]+), Max\. Consistent Triple Set heuristic: ([A-Za-z]+), Time: ([0-9]+\.[0-9]+)")
    
    with open(filename, "r") as f:
        
        line = f.readline()
        while line:
            match = info_line.match(line)
            if match:
                tree = f.readline().strip()
                cotree_mode = match.group(1)
                triple_mode = match.group(2)
                time = float(match.group(3))
                trees.append((tree, cotree_mode, triple_mode, time))
            line = f.readline()
    
    return trees


def remove_planted_root(newick):
    
    return newick[1:-2] + ";"

    
global_directory = os.path.expanduser("~") + "/Documents/PhD_local_files/ParaphyloVsHeuristic"

scenario_dir = os.path.join(global_directory, "proteomes")
result_dir = os.path.join(global_directory, "species_trees_heuristic")
 
result = {}
directories = sorted(glob.glob(scenario_dir + "/*"))
counter = 0
for directory in directories:
    
    i = int(directory.rsplit("_", 1)[1])
    counter += 1
    print(f"Processing directory {counter}/{len(directories)}: {directory}")
    
    true_tree_file = os.path.join(directory, "true_species_tree")
    with open(true_tree_file, "r") as f:
        true_tree = f.readline().strip()
        true_tree = remove_planted_root(true_tree)
        print(true_tree)
    
    paraphylo_file = os.path.join(directory, "species")
    build_tree, lrt_tree, paraphylo_time = parse_paraphylo_output(paraphylo_file)
    print(build_tree)
    
    heuristic_file = os.path.join(result_dir, f"species_tree_{i}")
    heuristic_result = parse_heuristic(heuristic_file)
    heuristic_trees = [item[0] for item in heuristic_result]
    print(heuristic_trees)
    
    result[i] = (true_tree, build_tree, lrt_tree, *heuristic_trees)
    
result_file = os.path.join(global_directory, "all_trees2")

with open(result_file, "w") as f:
    
    for key in sorted(result.keys()):
        for tree in result[key]:
            f.write(tree + "\n")
     
