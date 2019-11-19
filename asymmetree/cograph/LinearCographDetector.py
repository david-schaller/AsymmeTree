# -*- coding: utf-8 -*-

"""
Linear cograph detector.

Implementation of:
    D. G. Corneil, Y. Perl, and L. K. Stewart.
    A Linear Recognition Algorithm for Cographs.
    SIAM J. Comput., 14(4), 926–934 (1985).
    DOI: 10.1137/0214065
"""

from collections import deque

import tools.DoublyLinkedList as dll

try:
    from .Cograph import Cotree, CotreeNode
except ModuleNotFoundError:
    from cograph.Cograph import Cotree, CotreeNode


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


class LCDNode(CotreeNode):
    
    def __init__(self, ID, label=None, parent=None):
        super().__init__(ID, label=label, parent=parent)
        self.children = dll.DLList()
        self.parent_dll_element = None      # reference to doubly-linked list element
                                            # in the parents' children
        self.md = 0                         # number of marked and unmarked children


class LCD:
    
    def __init__(self, G):
        
        self.G = G
        self.V = list(G.adj_list.keys())
        
        self.T = Cotree(None)
        self.already_in_T = set()
        self.leaf_map = {}
        self.node_counter = 0
        
        self.marked = set()
        self.m_u_children = {}              # lists of marked and unmarked children
        self.mark_counter = 0
        self.unmark_counter = 0
        
        self.error_message = ""
    
    
    def cograph_recognition(self):
        
        if len(self.V) == 0:
            print("Empty graph in Cograph Recognition!")
            return self.T
        
        elif len(self.V) == 1:
            self.T.root = LCDNode(self.V[0], label="leaf")
            return self.T
        
        v1, v2 = self.V[0], self.V[1]
        self.already_in_T.update([v1, v2])
        
        R = LCDNode(None, label="series")
        self.T.root = R
        
        if self.G.has_edge(v1, v2):
            v1_node = LCDNode(v1, label="leaf", parent=R)
            v2_node = LCDNode(v2, label="leaf", parent=R)
            v1_node.parent_dll_element = R.children.append(v1_node)
            v2_node.parent_dll_element = R.children.append(v2_node)
            self.node_counter = 3
        else:
            N = LCDNode(None, label="parallel", parent=R)
            N.parent_dll_element = R.children.append(N)
            v1_node = LCDNode(v1, label="leaf", parent=N)
            v2_node = LCDNode(v2, label="leaf", parent=N)
            v1_node.parent_dll_element = N.children.append(v1_node)
            v2_node.parent_dll_element = N.children.append(v2_node)
            self.node_counter = 4
            
        self.leaf_map[v1] = v1_node
        self.leaf_map[v2] = v2_node
        
        if len(self.V) == 2:
            return self.T
        
        for x in self.V[2:]:
            
#            print(str(x)+" before "+self.T.to_newick())
            
            # initialization (necessary?)
            self.marked.clear()
            self.m_u_children.clear()
            self.mark_counter = 0
            self.unmark_counter = 0
            self.already_in_T.add(x)        # add x for subsequent iterations
            
            # call procedure _MARK(x)
            self._MARK(x)
#            print(self.marked)
#            print(self.node_counter)
#            print(self.mark_counter)
#            print(self.unmark_counter)
            
            # all nodes in T were marked and unmarked
            if self.node_counter == self.unmark_counter:
                R = self.T.root
                x_node = LCDNode(x, label="leaf", parent=R)
                x_node.parent_dll_element = R.children.append(x_node)
                self.node_counter += 1
                self.leaf_map[x] = x_node
#                print(str(x)+" after "+self.T.to_newick())
                continue
            # no nodes in T were marked and unmarked
            elif self.mark_counter == 0:
                # d(R)=1
                if len(self.T.root.children) == 1:
                    N = self.T.root.children[0]
                    x_node = LCDNode(x, label="leaf", parent=N)
                    x_node.parent_dll_element = N.children.append(x_node)
                    self.node_counter += 1
                else:
                    R_old = self.T.root
                    R_new = LCDNode(None, label="series")
                    N = LCDNode(None, label="parallel", parent=R_new)
                    N.parent_dll_element = R_new.children.append(N)
                    R_old.parent = N
                    R_old.parent_dll_element = N.children.append(R_old)
                    self.T.root = R_new
                    
                    x_node = LCDNode(x, label="leaf", parent=N)
                    x_node.parent_dll_element = N.children.append(x_node)
                    self.node_counter += 3
                self.leaf_map[x] = x_node
#                print(str(x)+" after "+self.T.to_newick())
                continue
            
            u = self._find_lowest()
            if not u:
                return False
#            print(u)
            
            # label(u)=0 and |A|=1
            if u.label == "parallel" and len(self.m_u_children[u]) == 1:
                w = self.m_u_children[u][0]
                if w.label == "leaf":
                    new_node = LCDNode(None, label="series", parent=u)
                    u.children.remove_element(w.parent_dll_element)
                    new_node.parent_dll_element = u.children.append(new_node)
                    w.parent = new_node
                    w.parent_dll_element = new_node.children.append(w)
                    
                    x_node = LCDNode(x, label="leaf", parent=new_node)
                    x_node.parent_dll_element = new_node.children.append(x_node)
                    self.node_counter += 2
                else:
                    x_node = LCDNode(x, label="leaf", parent=w)
                    x_node.parent_dll_element = w.children.append(x_node)
                    self.node_counter += 1 
            
            # label(u)=1 and |B|=1
            elif (u.label == "series" and 
                  len(u.children) - len(self.m_u_children[u]) == 1):
                set_A = set(self.m_u_children[u])       # auxiliary set bounded by O(deg(x))
                w = None
                for child in u.children:
                    if child not in set_A:
                        w = child
                        break
                if w.label == "leaf":
                    new_node = LCDNode(None, label="parallel", parent=u)
                    u.children.remove_element(w.parent_dll_element)
                    new_node.parent_dll_element = u.children.append(new_node)
                    w.parent = new_node
                    w.parent_dll_element = new_node.children.append(w)
                    
                    x_node = LCDNode(x, label="leaf", parent=new_node)
                    x_node.parent_dll_element = new_node.children.append(x_node)
                    self.node_counter += 2
                else:
                    x_node = LCDNode(x, label="leaf", parent=w)
                    x_node.parent_dll_element = w.children.append(x_node)
                    self.node_counter += 1
            
            else:
                y = LCDNode(None, label=u.label)
                for a in self.m_u_children[u]:
                    u.children.remove_element(a.parent_dll_element)
                    a.parent = y
                    a.parent_dll_element = y.children.append(a)
                    
                if u.label == "parallel":
                    new_node = LCDNode(None, label="series", parent=u)
                    new_node.parent_dll_element = u.children.append(new_node)
                    
                    y.parent = new_node
                    y.parent_dll_element = new_node.children.append(y)
                    x_node = LCDNode(x, label="leaf", parent=new_node)
                    x_node.parent_dll_element = new_node.children.append(x_node)
                else:
                    par = u.parent
                    if par is not None:             # u was the root of T
                        par.children.remove_element(u.parent_dll_element)
                        y.parent_dll_element = par.children.append(y)
                    else:
                        self.T.root = y             # y becomes the new root
                    y.parent = par
                    
                    new_node = LCDNode(None, label="parallel", parent=y)
                    new_node.parent_dll_element = y.children.append(new_node)
                    u.parent = new_node
                    u.parent_dll_element = new_node.children.append(u)
                    x_node = LCDNode(x, label="leaf", parent=new_node)
                    x_node.parent_dll_element = new_node.children.append(x_node)
                self.node_counter += 3
                
            self.leaf_map[x] = x_node
#            print(str(x)+" after "+self.T.to_newick())
        
        return self.T
    
    
    def _MARK(self, x):
        
        for v in self.G.neighbors(x):
            if v in self.already_in_T:
                self.marked.add(self.leaf_map[v])
                self.mark_counter += 1
                
        queue = deque(self.marked)
        
        while queue:                        # contains only d(u)=md(u) nodes
            u = queue.popleft()
            self.marked.remove(u)           # unmark u
            self.unmark_counter += 1
            u.md = 0                        # md(u) <- 0
            if u is not self.T.root:
                w = u.parent                # w <- parent(u)
                if w not in self.marked:
                    self.marked.add(w)      # mark w
                    self.mark_counter += 1
                w.md += 1
                if w.md == len(w.children):
                    queue.append(w)
                    
                if w in self.m_u_children:              # append u to list of
                    self.m_u_children[w].appendleft(u)  # marked and unmarked
                else:                                   # children of w
                    self.m_u_children[w] = deque([u])
                    
        if (self.marked and                             # any vertex remained marked
            len(self.T.root.children) == 1 and 
            self.T.root not in self.marked):
            
            self.marked.add(self.T.root)
            self.mark_counter += 1
    
    
    def _find_lowest(self):
        
        R = self.T.root
        y = "Lambda"
        
        if R not in self.marked:        # R is not marked
            self.error_message = "(iii): R=" + str(R)
            return False                # G+x is not a cograph (iii)
        else:
            if R.md != len(R.children) - 1:
                y = R
            self.marked.remove(R)
            R.md = 0
            u = w = R
        
        while self.marked:              # while there are mark vertices
            u = self.marked.pop()       # choose a arbitrary marked vertex u
            
            if y != "Lambda":
                self.error_message = "(i) or (ii): y=" + str(y)
                return False            # G+x is not a cograph (i) or (ii)
            
            if u.label == "series":
                if u.md != len(u.children) - 1:
                    y = u
                if u.parent in self.marked:
                    self.error_message = "(i) and (vi): u=" + str(u)
                    return False        # G+x is not a cograph (i) and (vi)
                else:
                    t = u.parent.parent
            else:
                y = u
                t = u.parent
            u.md = 0                    # u was already unmarked above
            
            # check if the u-w path is part of the legitimate alternating path
            while t is not w:
                if t is R:
                    self.error_message = "(iv): t=" + str(t)
                    return False        # G+x is not a cograph (iv)
                
                if t not in self.marked:
                    self.error_message = "(iii), (v) or (vi): t=" + str(t)
                    return False        # G+x is not a cograph (iii), (v) or (vi)
                
                if t.md != len(t.children) - 1:
                    self.error_message = "(ii): t=" + str(t)
                    return False        # G+x is not a cograph (ii)
                
                if t.parent in self.marked:
                    self.error_message = "(i): t=" + str(t)
                    return False        # G+x is not a cograph (i)
                
                self.marked.remove(t)   # unmark t
                t.md = 0                # reset md(t)
                t = t.parent.parent
                
            w = u                       # rest w for next choice of marked vertex
        
        return u


if __name__ == "__main__":
    
    from cograph.Cograph import SimpleGraph
    cotree = Cotree.random_cotree(1000)
#    print(cotree.to_newick())
    cograph = cotree.to_cograph()
    
#    # (1,((5,(7,8)<0>)<1>,4)<0>)<1>;
#    cograph = SimpleGraph(initial={1: {8, 4, 5, 7}, 5: {8, 1, 7}, 7: {1, 5}, 8: {1, 5}, 4: {1}})
    #print(cograph.adj_list)
    
    LCD = LCD(cograph)
    new_cotree = LCD.cograph_recognition()
    print("done")
    if new_cotree:
#        print(new_cotree.to_newick())
        new_cograph = new_cotree.to_cograph()
        #print(new_cograph.adj_list)
        print(cograph.graphs_equal(new_cograph))
    else:
        print("Not a cograph!")
    