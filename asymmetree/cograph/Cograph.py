# -*- coding: utf-8 -*-

import itertools, random
import networkx as nx

from asymmetree.tools import Tree


__author__ = 'David Schaller'
    

class CotreeNode(Tree.TreeNode):
    
    __slots__ = ['aux_counter']
    
    
    def __init__(self, ID, label=None, parent=None,
                 aux_counter=0):
        
        super().__init__(ID, label=label)
        
        self.aux_counter = aux_counter      # auxiliary counting variable
        
    
    def __str__(self):
        
        if not self.children:
            return str(self.ID)
        
        elif self.label == 'series':
            return '<1>'
        
        elif self.label == 'parallel':
            return '<0>'
        
        else:
            return '<>'
        

class Cotree(Tree.Tree):
    
    # corresponding node type
    node_type = CotreeNode
    
    
    def __init__(self, root):
        super().__init__(root)
        
    
    def to_cograph(self):
        
        self.supply_leaves()
        G = nx.Graph()
        
        for v in self.root.leaves:
            G.add_node(v.ID)
        
        for u in self.preorder():
            if u.label == 'series':
                for v1, v2 in itertools.combinations(u.children, 2):
                    for l1, l2 in itertools.product(v1.leaves, v2.leaves):
                        G.add_edge(l1.ID, l2.ID)
        
        return G
        
        
    @staticmethod
    def cotree(G):
        """Checks if a graph is a cograph and returns its cotree.
        
        Simple O(n^3) implementation.
        """
        
        def build_cotree(G, label=None):
            
            v = CotreeNode(None)
            v.label = label
            child_nodes = []
            
            if G.order() == 1:
                v.label = 'leaf'
                for ID in G.nodes():
                    v.ID = ID
                    
                return v
            
            child_nodes = []
            
            if v.label is None:
                ccs = [cc for cc in nx.connected_components(G)]
                
                if len(ccs) > 1:
                    v.label = 'parallel'
                else:
                    G = nx.complement(G)
                    ccs = [cc for cc in nx.connected_components(G)]
                    if len(ccs) > 1:
                        v.label = 'series'
                    else:
                        return False
            
            else:
                G = nx.complement(G)
                ccs = [cc for cc in nx.connected_components(G)]
                
                if len(ccs) == 1:
                    return False
            
            child_label = 'series' if v.label == 'parallel' else 'parallel'
            for cc in ccs:
                G_i = G.subgraph(cc).copy()
                v_i = build_cotree(G_i, label=child_label)
                if v_i:
                    child_nodes.append(v_i)
                else:
                    return False
            
            for child in child_nodes:
                v.children.append(child)
                child.parent = v
            
            return v
        
        root = build_cotree(G)
        
        if root:
            return Cotree(root)
        else:
            return False
       
    
    @staticmethod
    def random_cotree(N_leaves, force_series_root=False):
        """Creates a random cotree."""
        
        root = CotreeNode(ID=0)
        cotree = Cotree(root)
        node_list = [root]
        nr, leaf_count = 1, 1
        
        while leaf_count < N_leaves:
            node = random.choice(node_list)
            
            # avoid nodes with outdegree 1
            if not node.children:
                new_child1 = CotreeNode(nr)
                new_child2 = CotreeNode(nr+1)
                node.add_child(new_child1)
                node.add_child(new_child2)
                node_list.extend(node.children)
                nr += 2
                leaf_count += 1
                
            # add only one child if there are already children
            elif node.children:
                new_child = CotreeNode(nr)
                node.add_child(new_child)
                node_list.append(new_child)
                nr += 1
                leaf_count += 1
        
        # assign labels ('series', 'parallel', 'leaf')
        for v in cotree.preorder():
            if not v.children:
                v.label == 'leaf'
            elif v.parent is None:
                if force_series_root:
                    v.label = 'series'
                else:
                    v.label = 'series' if random.random() < 0.5 else 'parallel'
            else:
                v.label = 'series' if v.parent.label == 'parallel' else 'parallel'
                
        return cotree