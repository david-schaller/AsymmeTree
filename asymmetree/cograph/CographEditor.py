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
        
    
    def cograph_edit(self):
        
        # build starting tree for the first two vertices
        self._start_tree()
        print(self.T.to_newick())
        
        if len(self.V) <= 2:
            return self.T
        
        # incrementally insert vertices while updating total cost
        for x in self.V[2:]:
            cost, x_node = self._insert(x)
            
            # update the number of leaves on the path to the root
            current = x_node.parent
            while current:
                current.leaf_number += 1
                current = current.parent
            
            self.leaf_map[x] = x_node
            self.already_in_T.add(x)
            self.total_cost += cost
            print(x, self.T.to_newick())
            
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
            R.children.append(v1_node)
            R.children.append(v2_node)
        else:
            N = CENode(None, label="parallel", parent=R, leaf_number=2)
            R.children.append(N)
            v1_node = CENode(v1, label="leaf", parent=N, leaf_number=1)
            v2_node = CENode(v2, label="leaf", parent=N, leaf_number=1)
            N.children.append(v1_node)
            N.children.append(v2_node)
            
        self.leaf_map[v1] = v1_node
        self.leaf_map[v2] = v2_node
    
        
    def _insert(self, x):
        
        # information for non-hollow nodes
        C_nh = {}                   # node u: list of non-hollow children
        C_nh_number = {}            # node u: number of non-hollow children
        
        Nx_number = {}              # node u: number of neighbors of x in V(u)
        completion_forced = {}      # node u: completion-forced or not
        full = {}                   # node u: full or not (i.e. mixed)
        deletion_forced = {}        # node u: deletion-forced or not
        
        Nx_counter = 0
        
        # ---- FIRST STEP ----
        
        # first traversal (C_nh)
        for v in self.G.neighbors(x):
            if v not in self.already_in_T:
                continue
            
            Nx_counter += 1
            
            v = self.leaf_map[v]
            C_nh[v]                = []
            C_nh_number[v]         = 0
            Nx_number[v]           = 1
            completion_forced[v]   = True
            full[v]                = True
            deletion_forced[v]     = False
            
            current = v
            while current.parent:
                if current.parent in C_nh:
                    C_nh[current.parent].append(current)
                    break                               # rest of the path is already done
                else:
                    C_nh[current.parent] = [current]
                current = current.parent                # continue path to root
                
        # all nodes in T are adjacent to x
        if self.T.root.leaf_number == Nx_counter:
            R = self.T.root
            x_node = CENode(x, label="leaf", parent=R, leaf_number=1)
            R.children.append(x_node)
            self.leaf_map[x] = x_node
            return 0, x_node                            # cost is 0
        # no nodes in T are adjacent to x
        elif Nx_counter == 0:
            # d(R)=1
            if len(self.T.root.children) == 1:
                N = self.T.root.children[0]
                x_node = CENode(x, label="leaf", parent=N, leaf_number=1)
                N.children.append(x_node)
            else:
                R_old = self.T.root
                R_new = CENode(None, label="series", leaf_number=R_old.leaf_number)
                N = CENode(None, label="parallel", parent=R_new, leaf_number=R_old.leaf_number)
                R_new.children.append(N)
                R_old.parent = N
                N.children.append(R_old)
                self.T.root = R_new
                
                x_node = CENode(x, label="leaf", parent=N, leaf_number=1)
                N.children.append(x_node)
            return 0, x_node                            # cost is 0
                
        # second traversal, postorder
        # (C_nh_number, Nx_number, completion_forced, mixed, deletion_forced)
        stack = [self.T.root]
        while stack:
            u = stack.pop()
            if not u.children:
                continue
            elif u not in C_nh_number:
                stack.append(u)
                C_nh_number[u] = len(C_nh[u])
                for nh_child in C_nh[u]:
                    stack.append(nh_child)
            else:
                completion_forced_counter = 0
                deletion_forced_counter = 0
                full_counter = 0
                Nx_number[u] = 0
                for nh_child in C_nh[u]:
                    Nx_number[u] += Nx_number[nh_child]
                    if completion_forced[nh_child]:
                        completion_forced_counter += 1
                    if deletion_forced[nh_child]:
                        deletion_forced_counter += 1
                    if full[nh_child]:
                        full_counter += 1
                full[u] = True if Nx_number[u] == u.leaf_number else False
                
                # u completion-forced?
                if (full[u] or
                    (u.label == "parallel" and len(u.children) == C_nh_number[u]) or
                    (u.label == "series" and len(u.children) == completion_forced_counter)):
                    completion_forced[u] = True
                else:
                    completion_forced[u] = False
                    
                # u deletion-forced?
                if ((u.label == "series" and full_counter == 0) or
                    (u.label == "parallel" and C_nh_number[u] == deletion_forced_counter)):
                    deletion_forced[u] = True
                else:
                    deletion_forced[u] = False
                    
        # ---- SECOND STEP ----
        cost_above = {}                     # node u: cost above u
        C_mixed = {}                        # node u: mixed children of u
        C_full = {}                         # node u: full children of u
        
        stack = [self.T.root]
        while stack:        
            u = stack.pop()
            C_mixed[u] = []
            C_full[u] = []
            for nh_child in C_nh[u]:
                if not full[nh_child]:      # only mixed nodes are considered
                    stack.append(nh_child)
                    C_mixed[u].append(nh_child)
                else:
                    C_full[u].append(nh_child)
            v = u.parent
            if not v:                       # root has no parent
                cost_above[u] = 0
            elif u.parent.label == "parallel":
                cost_above[u] = cost_above[v] + Nx_number[v] - Nx_number[u]
            else:
                cost_above[u] = cost_above[v] + ((v.leaf_number - Nx_number[v]) - 
                                                 (u.leaf_number - Nx_number[u]))
        
        mincost = {}
        
        for u in cost_above.keys():         # contains exactly the mixed nodes
            
            # apply Lemma 7 directly if u has exactly 2 children
            if len(u.children) == 2 and u.label == "parallel":
                for i in range(2):
                    v, v2 = u.children[i], u.children[-1-i]
                    if v in completion_forced and completion_forced[v]:
                        Nx_v2 = Nx_number[v2] if (v2 in Nx_number) else 0
                        cost = cost_above[u] + (v.leaf_number - Nx_number[v]) + Nx_v2
                        if (u not in mincost) or (cost < mincost[u][0]):
                            mincost[u] = (cost, [v])
            elif len(u.children) == 2 and u.label == "series":
                for i in range(2):
                    v, v2 = u.children[i], u.children[-1-i]
                    if v not in deletion_forced or deletion_forced[v]:
                        Nx_v = Nx_number[v] if (v in Nx_number) else 0
                        Nx_v2 = Nx_number[v2] if (v2 in Nx_number) else 0
                        cost = cost_above[u] + Nx_v + (v2.leaf_number - Nx_v2)
                        if (u not in mincost) or (cost < mincost[u][0]):
                            mincost[u] = (cost, [v2])
            
            # u has more than 2 children
            else:
                red = set()
                blue = set()
                
                if u.label == "parallel":                                   # u is a parallel node
                    # first step
                    for v in C_mixed[u]:
                        if Nx_number[v] >= v.leaf_number - Nx_number[v]:
                            red.add(v)
                        else:
                            blue.add(v)
                    
                    # second step
                    if C_nh_number[u] == len(u.children) and not blue:
                        current_min, current_min_node = float("inf"), None
                        for v in red:
                            diff = 2 * Nx_number[v] - v.leaf_number         # = Nx_number[v] - (v.leaf_number - Nx_number[v])
                            if diff < current_min:
                                current_min, current_min_node = diff, v
                        red.remove(current_min_node)
                        blue.add(current_min_node)
                        
                    # third step
                    nb_filled = len(C_full[u]) + len(red)
                    if nb_filled < 2:
                        for i in range(2 - nb_filled):
                            current_min, current_min_node = float("inf"), None
                            for v in blue:
                                diff = v.leaf_number - 2 * Nx_number[v]     # = (v.leaf_number - Nx_number[v]) - Nx_number[v]
                                if diff < current_min:
                                    current_min, current_min_node = diff, v
                            blue.remove(current_min_node)
                            red.add(current_min_node)
                else:                                                       # u is a series node
                    # first step
                    for v in C_mixed[u]:
                        if v.leaf_number - Nx_number[v] >= Nx_number[v]:
                            blue.add(v)
                        else:
                            red.add(v)
                    
                    # second step
                    if not C_full[u] and not red:
                        current_min, current_min_node = float("inf"), None
                        for v in blue:
                            diff = v.leaf_number - 2 * Nx_number[v]         # = (v.leaf_number - Nx_number[v]) - Nx_number[v]
                            if diff < current_min:
                                current_min, current_min_node = diff, v
                        blue.remove(current_min_node)
                        red.add(current_min_node)
                        
                    # third step
                    nb_hollowed = len(u.children) - C_nh_number[u] + len(blue)    # number of hollow or blue children
                    if nb_hollowed < 2:
                        for i in range(2 - nb_hollowed):
                            current_min, current_min_node = float("inf"), None
                            for v in red:
                                diff = 2 * Nx_number[v] - v.leaf_number      # = Nx_number[v] - (v.leaf_number - Nx_number[v])
                                if diff < current_min:
                                    current_min, current_min_node = diff, v
                            red.remove(current_min_node)
                            blue.add(current_min_node)
                    
                cost = cost_above[u]
                for v in red:
                    cost += (v.leaf_number - Nx_number[v])
                for v in blue:
                    cost += Nx_number[v]
                    
                mincost[u] = (cost, C_full[u] + list(red[u]))
        
        # determine a minimum cost settling
        insertion_mincost, settling_node = float("inf"), None
        for u in mincost.keys():
            if mincost[u][0] < insertion_mincost:
                insertion_mincost, settling_node = mincost[u][0], u
        
        # insert x into tree
        u = settling_node           # node under which x is inserted
        print("settling node:", u)
        filled = mincost[u][1]      # children of u to be filled (marked nodes in LinearCographDetector)
        print("filled", filled[0].ID)
        
        # parallel node where one child is to be filled
        if u.label == "parallel" and len(filled) == 1:
            w = filled[0]
            if w.label == "leaf":
                new_node = CENode(None, label="series", parent=u, leaf_number=2)  # leaves w and x
                u.children.remove(w)
                u.children.append(new_node)
                w.parent = new_node
                new_node.children.append(w)
                
                x_node = CENode(x, label="leaf", parent=new_node, leaf_number=1)
                new_node.children.append(x_node)
            else:
                x_node = CENode(x, label="leaf", parent=w, leaf_number=1)
                w.children.append(x_node)
        
        # series node and only one child is not to be filled
        elif (u.label == "series" and 
              len(u.children) - len(filled) == 1):
            set_A = set(filled)       # auxiliary set
            w = None
            for child in u.children:
                if child not in set_A:
                    w = child
                    break
            if w.label == "leaf":
                new_node = CENode(None, label="parallel", parent=u, leaf_number=2)  # leaves w and x
                u.children.remove(w)
                u.children.append(new_node)
                w.parent = new_node
                new_node.children.append(w)
                
                x_node = CENode(x, label="leaf", parent=new_node, leaf_number=1)
                new_node.children.append(x_node)
            else:
                x_node = CENode(x, label="leaf", parent=w, leaf_number=1)
                w.children.append(x_node)
        
        else:
            y = CENode(None, label=u.label, leaf_number=len(filled))
            for a in filled:
                u.children.remove(a)
                a.parent = y
                y.children.append(a)
                
            if u.label == "parallel":
                new_node = CENode(None, label="series", parent=u, leaf_number=len(filled))
                u.children.append(new_node)
                
                y.parent = new_node
                new_node.children.append(y)
                x_node = CENode(x, label="leaf", parent=new_node, leaf_number=1)
                new_node.children.append(x_node)
            else:
                par = u.parent
                if par is not None:             # u was the root of T
                    par.children.remove(u)
                    par.children.append(y)
                else:
                    self.T.root = y             # y becomes the new root
                y.parent = par
                
                new_node = CENode(None, label="parallel", parent=y)
                y.children.append(new_node)
                u.parent = new_node
                new_node.children.append(u)
                x_node = CENode(x, label="leaf", parent=new_node, leaf_number=1)
                new_node.children.append(x_node)
                
                # update the leaf numbers
                y.leaf_number = u.leaf_number
                u.leaf_number -= len(filled)
                new_node.leaf_number = u.leaf_number

        return insertion_mincost, x_node
                
    
    
if __name__ == "__main__":
    
#    cotree = Cotree.random_cotree(5)
#    print(cotree.to_newick())
#    cograph = cotree.to_cograph()
#    print(cograph.adj_list)
    
    from cograph.Cograph import SimpleGraph
    print("(1,2,(4,5,6)<0>)<1>;")
    cograph = SimpleGraph(initial={1: {2, 4, 5, 6}, 2: {1, 4, 5, 6}, 4: {1, 2}, 5: {1, 2}, 6: {1, 2}})
    print(cograph.adj_list)
    
    CE = CographEditor(cograph)
    new_cotree = CE.cograph_edit()
    print("done")
    if new_cotree:
        print(new_cotree.to_newick())
        new_cograph = new_cotree.to_cograph()
        print(new_cograph.adj_list)
        print(CE.total_cost)
        print(cograph.graphs_equal(new_cograph))
    else:
        print("Not a cograph!")