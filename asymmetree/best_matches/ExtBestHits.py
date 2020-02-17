# -*- coding: utf-8 -*-

"""
Extended Best Hits Method.

Implementation of the Extended Best Hits method for best match inference.

Methods in this module:
    - ebh_qinfer
    - ebh_from_scenario
    - ebh
"""

import os, subprocess, time
import networkx as nx

from asymmetree.tools import FileIO
from asymmetree.best_matches import TrueBMG


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"



# --------------------------------------------------------------------------
#                           PYTHON IMPLEMENTATION
#
# --------------------------------------------------------------------------
    
def ebh(leaves, D, epsilon=0.000_000_01):
    """Compute BMG and RBMG from a distances matrix D.
    
    Keyword arguments:
        epsilon -- epsilon for relative BM threshold: (x,y) in BMG if
                   D(x,y) <= (1+epsilon) * min d(x,y'),
                   default=10E-8 (for limited float precision).
    """
    BMG, RBMG = nx.DiGraph(), nx.Graph()
    colors = set()
    relative_threshold = 1 + epsilon
    
    for v in leaves:
        BMG.add_node(v.ID, label=v.label, color=v.color)
        RBMG.add_node(v.ID, label=v.label, color=v.color)
        colors.add(v.color)
    
    # ---- build BMG ----
    for u in range(len(leaves)):
        minima = {color: float('inf') for color in colors}
        for v in range(len(leaves)):
            if D[u,v] < minima[leaves[v].color]:
                minima[leaves[v].color] = D[u,v]
        for v in range(len(leaves)):
            if (leaves[u].color != leaves[v].color and
                D[u,v] <= relative_threshold * minima[leaves[v].color]):
                BMG.add_edge(leaves[u].ID, leaves[v].ID, distance = D[u,v])
    
    # ---- build RBMG as symmetric part of the BMG ----
    for x, neighbors in BMG.adjacency():
        for y in neighbors:
            if BMG.has_edge(y,x):
                RBMG.add_edge(x,y)
    
    return BMG, RBMG


# --------------------------------------------------------------------------
#                            EXTERNAL C++ PROGRAM
#
# --------------------------------------------------------------------------
 
def ebh_qinfer(scenario,
               matrix_filename, species_filename,
               epsilon=0.5,
               benchmark_file=None,
               binary_path=None):
    """Compute BMG and RBMG from a distances matrix D using 'qinfer'.
    
    Keyword arguments:
        epsilon -- epsilon for relative BM threshold: (x,y) in BMG if
                   D(x,y) <= (1+epsilon) * min d(x,y'),
                   default=10E-8 (for limited float precision).
        benchmark_file -- activate benchmarking in 'qinfer' and
                          specify the filename
        binary_path -- path to 'qinfer' binary (if not available
                       within path)
    """
    
    if not binary_path:
        qinfer_command = "qinfer"
    elif os.path.exists(binary_path):
        qinfer_command = binary_path
    else:
        raise FileNotFoundError(f"Path to qinfer binary file '{binary_path}' does not exist!")
    
    output = -1
    command = [qinfer_command, matrix_filename, species_filename,
               "--disable-quartet", "--epsilon=" + str(epsilon)]

    if benchmark_file is not None:
        command.append( "--benchmark=" + benchmark_file )
    
    # call 'qinfer' and measure execution time
    start = time.time()
    
    try:
        output = subprocess.run(command, stdout=subprocess.PIPE)
    except:
        raise FileNotFoundError("Calling qinfer failed!")
    
    exec_time = time.time() - start
    
    if output == -1:
        raise Exception("No output from qinfer!")
    
    BMG = FileIO.parse_BMG_edges(output.stdout.decode(), scenario)
    RBMG = TrueBMG.RBMG_from_BMG(BMG)
    
    return BMG, RBMG, exec_time


def ebh_from_scenario(scenario, epsilon=0.5):
    """Compute BMG and RBMG from a scenario using 'qinfer'.
    
    Keyword arguments:
        epsilon -- epsilon for relative BM threshold: (x,y) in BMG if
                   D(x,y) <= (1+epsilon) * min d(x,y'),
                   default=10E-8 (for limited float precision).
    """

    matrix_filename = "temp.phylip"
    species_filename = "temp_species.txt"
    
    matrix = scenario.get_distance_matrix()
    FileIO.matrix_to_phylip(matrix_filename, scenario.genes, matrix)
    FileIO.species_to_genes(species_filename, scenario)
    
    BMG, RBMG, exec_time = ebh_qinfer(scenario,
                                      matrix_filename, species_filename,
                                      epsilon=epsilon)
    
    os.remove(matrix_filename)
    os.remove(species_filename)
    
    return BMG, RBMG, exec_time