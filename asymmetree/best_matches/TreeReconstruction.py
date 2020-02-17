# -*- coding: utf-8 -*-

import os, subprocess, itertools, time

from asymmetree.tools import FileIO
from asymmetree.tools.PhyloTree import PhyloTree, PhyloTreeNode


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


# --------------------------------------------------------------------------
#                            TREE RECONSTRUCTION
#
# --------------------------------------------------------------------------

def neighbor_joining(leaves, leaf_index, matrix_filename,
                     return_calltime=False,
                     binary_path=None):
    """
    Keyword argument:
        return_calltime -- return the time needed for RapidNJ call
        binary_path -- path to 'qinfer' binary (if not available
                       within path)
    """
    
    if not binary_path:
        nj_command = "rapidnj"
    elif os.path.exists(binary_path):
         nj_command = binary_path
    else:
        raise FileNotFoundError(f"Path to RapidNJ binary file '{binary_path}' does not exist!")
            
    start_time = time.time()
    
    try:
        output = subprocess.run([nj_command, matrix_filename, "-i", "pd"],
                                stdout=subprocess.PIPE)
    except:
        raise FileNotFoundError("Calling RapidNJ failed!")
        
    calltime = time.time() - start_time
    
    newick = output.stdout.decode()
    
    tree = PhyloTree.parse_newick(newick)
    
    # restore (leaf) indeces and colors
    leaf_dict = {leaf.ID: leaf for leaf in leaves}
    index_counter = 0
    for v in tree.preorder():
        if not v.children:
            v.ID = int(v.label)
            v.color = leaf_dict[v.ID].color
        else:
            while index_counter in leaf_dict:
                index_counter += 1
            v.ID = index_counter
            index_counter += 1
    
    if return_calltime:
        return tree, calltime
    else:
        return tree
    

def reroot(tree, node):
    
    edges = []                              # edges (u, v, distance) to be changed
    pos = node
    while pos.parent:
        edges.append( (pos.parent, pos, pos.dist) )    
        pos = pos.parent
    old_root = pos
    
    for u, v, dist in reversed(edges):      # change direction of edges
        u.remove_child(v)
        v.add_child(u)
        u.dist = dist
    
    node.detach()
    node.dist = 0.0
    tree.root = node
    
    if len(old_root.children) <= 1:         # delete old root if its out-degree is 1
        tree.delete_and_reconnect(old_root)


def midpoint_rooting(tree):
    
    # identify the two most distant leaves and their l.c.a.
    distance_dict, _ = tree.distances_from_root()
    tree.supply_leaves()
    max_dist, leaf1, leaf2, lca = float('-inf'), None, None, None
    max_ID = 0
    
    for v in tree.preorder():
        if v.children:
            for c1, c2 in itertools.combinations(v.children, 2):
                for x in c1.leaves:
                    x_dist = distance_dict[x] - distance_dict[v]
                    for y in c2.leaves:
                        y_dist = distance_dict[y] - distance_dict[v]
                        if x_dist + y_dist > max_dist:
                            max_dist, leaf1, leaf2, lca = x_dist + y_dist, x, y, v
        max_ID = max([max_ID, v.ID])
    
    # identify the edge for the new root
    edge, cut_at = None, 0
    
    pos, remaining = leaf1, max_dist/2
    while (not edge) and (pos is not lca):
        if remaining < pos.dist:
            edge = (pos.parent, pos)
            cut_at = pos.dist - remaining
        remaining -= pos.dist
        pos = pos.parent
    
    pos, remaining = leaf2, max_dist/2
    while (not edge) and (pos is not lca):
        if remaining < pos.dist:
            edge = (pos.parent, pos)
            cut_at = pos.dist - remaining
        remaining -= pos.dist
        pos = pos.parent
    
    # rerooting
    if edge:
        u, v = edge
        u.remove_child(v)
        new_root = PhyloTreeNode(max_ID+1, dist=cut_at)
        u.add_child(new_root)
        new_root.add_child(v)
        v.dist -= new_root.dist
        reroot(tree, new_root)
    else:
        reroot(tree, lca)
        

def nj_from_numpy_matrix(leaves, leaf_index, matrix,
                         filename="temp.phylip"):

    FileIO.matrix_to_phylip(filename, leaves, matrix)
    tree = neighbor_joining(leaves, leaf_index, filename)
    os.remove(filename)
    
    return tree