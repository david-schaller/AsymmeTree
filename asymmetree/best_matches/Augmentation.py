# -*- coding: utf-8 -*-


import networkx as nx

from asymmetree import TreeNode, PhyloTreeNode


__author__ = 'David Schaller'


def augment_and_label(tree, inplace=False):
    """Augment tree and add event labeling based on color intersections."""
    
    if not inplace:
        tree = tree.copy()
    
    # precompute since tree will be modified
    inner_vertices = [u for u in tree.inner_vertices()]
    tree.supply_leaves()
    current_ID = tree.get_max_ID() + 1
    
    for u in inner_vertices:
        
        # should not appear except if tree is planted, otherwise ignore
        if len(u.children) == 1:
            continue
        
        C = _color_intersection_components(u)
        
        if len(C) == 1:
            u.label = 'D'
            
        else:
            u.label = 'S'
            
            for cc in C:
                if len(cc) > 1:
                    w = PhyloTreeNode(current_ID, label='D')
                    current_ID += 1
                    u.add_child(w)
                    for v in cc:
                        u.remove_child(v)
                        w.add_child(v)
        
    return tree


def _color_intersection_components(u):
    
    aux_graph = nx.Graph()
    
    for v in u.children:
        aux_graph.add_node(v)
        
        for x in v.leaves:
            aux_graph.add_edge(v, x.color)
            
    result = []
    
    for cc in nx.connected_components(aux_graph):
        
        # discard the colors in the components
        result.append( [v for v in cc if isinstance(v, TreeNode)] )
        
    return result


if __name__ == '__main__':
    
    from asymmetree import PhyloTree
    from asymmetree.best_matches.LRTConstructor import lrt_from_observable_tree
    
    T = PhyloTree.random_colored_tree(30, 5)
    print('--- T ---\n', T.to_newick())
    
    lrt = lrt_from_observable_tree(T)
    print('--- LRT ---\n', lrt.to_newick())
    
    aug_tree = augment_and_label(lrt)
    print('--- Augmented tree ---\n', aug_tree.to_newick())