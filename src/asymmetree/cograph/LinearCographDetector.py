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
import networkx as nx

from asymmetree.cograph.Cograph import Cotree, CotreeNode


__author__ = 'David Schaller'


def linear_cograph_detection(G, return_cotree=True):
    
    lcd = LinearCographDetector(G)
    cotree = lcd.recognition()
    
    if not cotree:
        return False
    else:
        return cotree if return_cotree else True


class LinearCographDetector:
    
    def __init__(self, G):
        
        if not isinstance(G, nx.Graph):
            raise TypeError('not a NetworkX Graph')
        
        self.G = G
        self.V = [v for v in G.nodes()]
        
        self.T = Cotree(None)
        self.already_in_T = set()
        self.leaf_map = {}
        self.node_counter = 0
        
        self.marked = set()
        self.m_u_children = {}              # lists of marked and unmarked children
        self.mark_counter = 0
        self.unmark_counter = 0
        
        self.error_message = ''
    
    
    def recognition(self):
        
        if len(self.V) == 0:
            raise RuntimeError('empty graph in cograph recognition')
            return self.T
        
        elif len(self.V) == 1:
            self.T.root = CotreeNode(self.V[0], label='leaf')
            return self.T
        
        v1, v2 = self.V[0], self.V[1]
        self.already_in_T.update([v1, v2])
        
        R = CotreeNode(None, label='series')
        self.T.root = R
        
        if self.G.has_edge(v1, v2):
            v1_node = CotreeNode(v1, label='leaf')
            v2_node = CotreeNode(v2, label='leaf')
            R.add_child(v1_node)
            R.add_child(v2_node)
            self.node_counter = 3
        else:
            N = CotreeNode(None, label='parallel')
            R.add_child(N)
            v1_node = CotreeNode(v1, label='leaf')
            v2_node = CotreeNode(v2, label='leaf')
            N.add_child(v1_node)
            N.add_child(v2_node)
            self.node_counter = 4
            
        self.leaf_map[v1] = v1_node
        self.leaf_map[v2] = v2_node
        
        if len(self.V) == 2:
            return self.T
        
        for x in self.V[2:]:
            
            # initialization (necessary?)
            self.marked.clear()
            self.m_u_children.clear()
            self.mark_counter = 0
            self.unmark_counter = 0
            self.already_in_T.add(x)        # add x for subsequent iterations
            
            # call procedure _MARK(x)
            self._MARK(x)
            
            # all nodes in T were marked and unmarked
            if self.node_counter == self.unmark_counter:
                R = self.T.root
                x_node = CotreeNode(x, label='leaf')
                R.add_child(x_node)
                self.node_counter += 1
                self.leaf_map[x] = x_node
                continue
            # no nodes in T were marked and unmarked
            elif self.mark_counter == 0:
                # d(R)=1
                if len(self.T.root.children) == 1:
                    N = self.T.root.children[0]
                    x_node = CotreeNode(x, label='leaf')
                    N.add_child(x_node)
                    self.node_counter += 1
                else:
                    R_old = self.T.root
                    R_new = CotreeNode(None, label='series')
                    N = CotreeNode(None, label='parallel')
                    R_new.add_child(N)
                    N.add_child(R_old)
                    self.T.root = R_new
                    
                    x_node = CotreeNode(x, label='leaf')
                    N.add_child(x_node)
                    self.node_counter += 3
                self.leaf_map[x] = x_node
                continue
            
            u = self._find_lowest()
            if not u:
                return False
            
            # label(u)=0 and |A|=1
            if u.label == 'parallel' and len(self.m_u_children[u]) == 1:
                w = self.m_u_children[u][0]
                if w.label == 'leaf':
                    new_node = CotreeNode(None, label='series')
                    u.remove_child(w)
                    u.add_child(new_node)
                    new_node.add_child(w)
                    
                    x_node = CotreeNode(x, label='leaf')
                    new_node.add_child(x_node)
                    self.node_counter += 2
                else:
                    x_node = CotreeNode(x, label='leaf')
                    w.add_child(x_node)
                    self.node_counter += 1 
            
            # label(u)=1 and |B|=1
            elif (u.label == 'series' and 
                  len(u.children) - len(self.m_u_children[u]) == 1):
                set_A = set(self.m_u_children[u])       # auxiliary set bounded by O(deg(x))
                w = None
                for child in u.children:
                    if child not in set_A:
                        w = child
                        break
                if w.label == 'leaf':
                    new_node = CotreeNode(None, label='parallel')
                    u.remove_child(w)
                    u.add_child(new_node)
                    new_node.add_child(w)
                    
                    x_node = CotreeNode(x, label='leaf')
                    new_node.add_child(x_node)
                    self.node_counter += 2
                else:
                    x_node = CotreeNode(x, label='leaf')
                    w.add_child(x_node)
                    self.node_counter += 1
            
            else:
                y = CotreeNode(None, label=u.label)
                for a in self.m_u_children[u]:
                    u.remove_child(a)
                    y.add_child(a)
                    
                if u.label == 'parallel':
                    new_node = CotreeNode(None, label='series')
                    u.add_child(new_node)
                    
                    new_node.add_child(y)
                    x_node = CotreeNode(x, label='leaf')
                    new_node.add_child(x_node)
                else:
                    par = u.parent
                    if par is not None:             # u was the root of T
                        par.remove_child(u)
                        par.add_child(y)
                    else:
                        self.T.root = y             # y becomes the new root
                    
                    new_node = CotreeNode(None, label='parallel')
                    y.add_child(new_node)
                    new_node.add_child(u)
                    x_node = CotreeNode(x, label='leaf')
                    new_node.add_child(x_node)
                self.node_counter += 3
                
            self.leaf_map[x] = x_node
        
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
            u.aux_counter = 0               # md(u) <- 0
            if u is not self.T.root:
                w = u.parent                # w <- parent(u)
                if w not in self.marked:
                    self.marked.add(w)      # mark w
                    self.mark_counter += 1
                w.aux_counter += 1
                if w.aux_counter == len(w.children):
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
        y = 'Lambda'
        
        if R not in self.marked:        # R is not marked
            self.error_message = '(iii): R={}'.format(R)
            return False                # G+x is not a cograph (iii)
        else:
            if R.aux_counter != len(R.children) - 1:
                y = R
            self.marked.remove(R)
            R.aux_counter = 0
            u = w = R
        
        while self.marked:              # while there are mark vertices
            u = self.marked.pop()       # choose a arbitrary marked vertex u
            
            if y != 'Lambda':
                self.error_message = '(i) or (ii): y={}'.format(y)
                return False            # G+x is not a cograph (i) or (ii)
            
            if u.label == 'series':
                if u.aux_counter != len(u.children) - 1:
                    y = u
                if u.parent in self.marked:
                    self.error_message = '(i) and (vi): u={}'.format(u)
                    return False        # G+x is not a cograph (i) and (vi)
                else:
                    t = u.parent.parent
            else:
                y = u
                t = u.parent
            u.aux_counter = 0           # u was already unmarked above
            
            # check if the u-w path is part of the legitimate alternating path
            while t is not w:
                if t is R:
                    self.error_message = '(iv): t={}'.format(t)
                    return False        # G+x is not a cograph (iv)
                
                if t not in self.marked:
                    self.error_message = '(iii), (v) or (vi): t={}'.format(t)
                    return False        # G+x is not a cograph (iii), (v) or (vi)
                
                if t.aux_counter != len(t.children) - 1:
                    self.error_message = '(ii): t={}'.format(t)
                    return False        # G+x is not a cograph (ii)
                
                if t.parent in self.marked:
                    self.error_message = '(i): t={}'.format(t)
                    return False        # G+x is not a cograph (i)
                
                self.marked.remove(t)   # unmark t
                t.aux_counter = 0       # reset md(t)
                t = t.parent.parent
                
            w = u                       # rest w for next choice of marked vertex
        
        return u
    