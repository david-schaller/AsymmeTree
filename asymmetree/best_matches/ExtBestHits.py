# -*- coding: utf-8 -*-

"""
Extended Best Hits Method.

Implementation of the Extended Best Hits method for best match inference.
"""

import os, subprocess, time
import networkx as nx

from asymmetree.file_io.ScenarioFileIO import parse_bmg_edges, matrix_to_phylip, species_to_genes
from asymmetree.tools.GraphTools import symmetric_part


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                           PYTHON IMPLEMENTATION
#
# --------------------------------------------------------------------------
    
def ebh(leaves, D, epsilon=1e-8):
    """Compute BMG and RBMG from a distances matrix D.
    
    Keyword arguments:
        epsilon -- epsilon for relative BM threshold: (x,y) in BMG if
                   D(x,y) <= (1+epsilon) * min d(x,y'),
                   default=1e-8 (for limited float precision).
    """
    
    bmg, rbmg = nx.DiGraph(), nx.Graph()
    colors = set()
    relative_threshold = 1 + epsilon
    
    for v in leaves:
        bmg.add_node(v.ID, label=v.label, color=v.color)
        rbmg.add_node(v.ID, label=v.label, color=v.color)
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
                bmg.add_edge(leaves[u].ID, leaves[v].ID, distance = D[u,v])
    
    return bmg, symmetric_part(bmg)


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
        raise FileNotFoundError("path to qinfer binary file '{}' does not exist".format(binary_path))
    
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
        raise FileNotFoundError("calling qinfer failed")
    
    exec_time = time.time() - start
    
    if output == -1:
        raise Exception("no output from qinfer")
    
    bmg = parse_bmg_edges(output.stdout.decode(), scenario)
    
    return bmg, symmetric_part(bmg), exec_time


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
    matrix_to_phylip(matrix_filename, scenario.genes, matrix)
    species_to_genes(species_filename, scenario)
    
    bmg, rbmg, exec_time = ebh_qinfer(scenario,
                                      matrix_filename, species_filename,
                                      epsilon=epsilon)
    
    os.remove(matrix_filename)
    os.remove(species_filename)
    
    return bmg, rbmg, exec_time