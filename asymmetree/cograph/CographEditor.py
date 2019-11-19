# -*- coding: utf-8 -*-

"""
Heuristic for cograph editing.

Implementation the O(n^2) algorithm in:
    Christophe Crespelle.
    Linear-Time Minimal Cograph Editing.
    Preprint published 2019.
"""

from collections import deque

import tools.DoublyLinkedList as dll

try:
    from .Cograph import Cotree, CotreeNode
except ModuleNotFoundError:
    from cograph.Cograph import Cotree, CotreeNode


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


class CENode(CotreeNode):
    """Treenode for cograph editing."""
    
    def __init__(self, ID, label=None, parent=None, leaf_number=0):
        super().__init__(ID, label=label, parent=parent)
#        self.children = dll.DLList()
#        self.parent_dll_element = None      # reference to doubly-linked list element
#                                            # in the parents' children
        self.leaf_number = leaf_number


class CographEditor:
    
    def __init__(self, G):
        
        self.G = G
        self.V = list(G.adj_list.keys())
        
        self.T = Cotree(None)
        self.already_in_T = set()
        self.leaf_map = {}
        
        self.total_cost = 0
        
#        self.node_counter = 0
        
#        self.marked = set()
#        self.m_u_children = {}              # lists of marked and unmarked children
#        self.mark_counter = 0
#        self.unmark_counter = 0
#        
#        self.error_message = ""
#    
    
    def cograph_edit(self):
        
        # build starting tree for the first two vertices
        self.start_tree()
        
        if len(self.V) <= 2:
            return self.T
        
        # incrementally insert vertices while updating total cost
        for x in self.V[2:]:
            self._insert(x)
            self.already_in_T.add(x)
            
        return self.T
        
    
    def _start_tree(self):
        
        if len(self.V) == 0:
            print("Empty graph in Cograph Recognition!")
            return self.T
        
        elif len(self.V) == 1:
            self.T.root = CENode(self.V[0], label="leaf", leaf_number=1)
            return self.T
        
        v1, v2 = self.V[0], self.V[1]
        self.already_in_T.update([v1, v2])
        
        R = CENode(None, label="series", leaf_number=2)
        self.T.root = R
        
        if self.G.has_edge(v1, v2):
            v1_node = CENode(v1, label="leaf", parent=R, leaf_number=1)
            v2_node = CENode(v2, label="leaf", parent=R, leaf_number=1)
            v1_node.parent_dll_element = R.children.append(v1_node)
            v2_node.parent_dll_element = R.children.append(v2_node)
#            self.node_counter = 3
        else:
            N = CENode(None, label="parallel", parent=R, leaf_number=2)
            N.parent_dll_element = R.children.append(N)
            v1_node = CENode(v1, label="leaf", parent=N, leaf_number=1)
            v2_node = CENode(v2, label="leaf", parent=N, leaf_number=1)
            v1_node.parent_dll_element = N.children.append(v1_node)
            v2_node.parent_dll_element = N.children.append(v2_node)
#            self.node_counter = 4
            
        self.leaf_map[v1] = v1_node
        self.leaf_map[v2] = v2_node
    
    
    def _get_nh_node_info(self, v, C_):
        pass
    
        
    def _insert(self, x):
        
        # information for non-hollow nodes
        C_nh = {}                   # node u: list of non-hollow children
        C_nh_number = {}            # node u: number of non-hollow children
        
        Nx_number = {}              # node u: number of neighbors of x in V(u)
        completion_forced = {}      # node u: completion-forced or not
        mixed = {}                  # node u: mixed or not (i.e. full)
        deletion_forced = {}        # node u: deletion-forced or not
        
        # ---- FIRST STEP ----
        
        # first traversal (C_nh, C_nh_number)
        for v in self.G.neighbors(x):
            if v not in self.already_in_T:
                continue
            
            v_node = self.leaf_map[v]
            C_nh[v_node]                = []
            C_nh_number[v_node]         = 0
            Nx_number[v_node]           = 1
            completion_forced[v_node]   = True
            mixed[v_node]               = False
            deletion_forced[v_node]     = False
            
            current = self.leaf_map[v]
            while current.parent:
                if current.parent in C_nh:
                    C_nh[current.parent].append(current)
                    C_nh_number[current.parent] += 1
                    break                               # rest of the path is already done
                else:
                    C_nh[current.parent] = [current]
                    C_nh_number[current.parent] = 1
                current = current.parent                # continue path to root
                
        # second traversal (C_nh, C_nh_number)
        pass
    
    
if __name__ == "__main__":
    
    from cograph.Cograph import SimpleGraph
    cotree = Cotree.random_cotree(5)
    print(cotree.to_newick())
    cograph = cotree.to_cograph()