# -*- coding: utf-8 -*-

import itertools, random

import asymmetree.datastructures.DoublyLinkedList as dll


__author__ = 'David Schaller'


class TreeNode:
    """Class 'TreeNode'.
    
    Components for class 'Tree'. Contains a list of children as well as a
    reference to the parent node.
    """
    
    __slots__ = ('ID', 'label', 'parent', 'parent_dll_element',
                 'children', 'leaves')
    
    
    def __init__(self, ID, label=''):
        
        self.ID = ID
        self.label = label
        
        
        self.parent = None
        # reference to doubly-linked list element in the parents' children
        self.parent_dll_element = None
        
        self.children = dll.DLList()
        
    
    def __str__(self):
        
        return str(self.label)
    
    
    def __repr__(self):
        
        return '<TN: {}>'.format(self.ID)
    
                                            
    def add_child(self, child_node):
        """Add a node as a child of this node."""
        
        # do nothing if child_node is already a child of self
        
        if child_node.parent is None:
            child_node.parent = self
            child_node.parent_dll_element = self.children.append(child_node)
        
        elif child_node.parent is not self:
            child_node.parent.remove_child(child_node)
            child_node.parent = self
            child_node.parent_dll_element = self.children.append(child_node)  
    
    
    def remove_child(self, child_node):
        
        if child_node.parent is self:
            self.children.remove_node(child_node.parent_dll_element)
            child_node.parent = None
            child_node.parent_dll_element = None
        else:
            raise ValueError('{} is not a child of node {}'.format(child_node,
                                                                   self))
            
            
    def detach(self):
        
        if self.parent is not None:
            self.parent.remove_child(self)
        else:
            self.parent = None
            self.parent_dll_element = None
            
    
    def is_leaf(self):
        """Return whether the node is a leaf."""
        
        return not self.children
        

class Tree:
    
    # corresponding node type
    node_type = TreeNode
    
    
    def __init__(self, root):
        
        self.root = root
        
    
    def leaves(self):
        """Generator for leaves of the tree."""
        
        def _leaves(node):
            if not node.children:
                yield node
            else:
                for child in node.children:
                    yield from _leaves(child)
        
        if self.root:
            yield from _leaves(self.root)
        else:
            yield from []
    
    
    def preorder(self):
        """Generator for preorder traversal of the tree."""
        
        def _preorder(node):
            yield node
            for child in node.children:
                yield from _preorder(child)
        
        if self.root:
            yield from _preorder(self.root)
        else:
            yield from []
            
    
    def traverse_subtree(self, u):
        """Generator for preorder traversal of the subtree rooted at u."""
        
        yield u
        for child in u.children:
            yield from self.traverse_subtree(child)
    
    
    def postorder(self):
        """Generator for postorder traversal of the tree."""
        
        def _postorder(node):
            for child in node.children:
                yield from _postorder(child)
            yield node
        
        if self.root:
            yield from _postorder(self.root)
        else:
            yield from []
            
    
    def inner_vertices(self):
        """Generator for inner vertices in preorder."""
        
        def _inner_vertices(node):
            if node.children:
                yield node
                for child in node.children:
                    yield from _inner_vertices(child)
        
        if self.root:
            yield from _inner_vertices(self.root)
        else:
            yield from []
            
    
    def edges(self):
        """Generator for all edges of the tree."""
        
        def _edges(node):
            for child in node.children:
                yield (node, child)
                yield from _edges(child)
        
        if self.root:
            yield from _edges(self.root)
        else:
            yield from []
            
            
    def edges_sibling_order(self):
        """Generator for all edges of the tree with sibling order.
        
        Returns edges uv as tuples (u, v, nr) where nr is the index of v in
        the list of children of node u."""
        
        def _edges_sibling_order(node):
            i = 0
            for child in node.children:
                yield (node, child, i)
                yield from _edges_sibling_order(child)
                i += 1
        
        if self.root:
            yield from _edges_sibling_order(self.root)
        else:
            yield from []
        
    
    def inner_edges(self):
        """Generator for all inner edges of the tree."""
        
        def _inner_edges(node):
            for child in node.children:
                if child.children:
                    yield (node, child)
                    yield from _inner_edges(child)
        
        if self.root:
            yield from _inner_edges(self.root)
        else:
            yield from []
                
    
    def euler_generator(self, id_only=False):
        """Generator for an Euler tour of the tree."""
        
        def _euler_generator(node, id_only):
            if id_only: yield node.ID
            else:       yield node
            
            for child in node.children:
                yield from _euler_generator(child, id_only)
                
                if id_only: yield node.ID
                else:       yield node
        
        if self.root:
            yield from _euler_generator(self.root, id_only)
        else:
            yield from []
            
        
    def euler_and_level(self):
        """Generator for an Euler tour with node levels."""
        
        def _euler_level(node, level):
            yield (node, level)
            
            for child in node.children:
                yield from _euler_level(child, level+1)
                yield (node, level)
        
        if self.root:
            yield from _euler_level(self.root, 0)
        else:
            yield from []
            
    
    def supply_leaves(self):
        """Add the leaves to all nodes that are in the subtree of a specific
        node."""
        
        def _supply_leaves(node):
            node.leaves = []
            
            if not node.children:
                node.leaves.append(node)
            else:
                for child in node.children:
                    node.leaves.extend(_supply_leaves(child))
                    
            return node.leaves
        
        if self.root:
            return _supply_leaves(self.root)
        else:
            return []
        
    
    def contract(self, edges):
        
        contracted = set()
        
        for u, v in edges:
            
            # avoid trying to contract the same edge multiple times
            if v not in contracted:
                self.delete_and_reconnect(v)
            
            contracted.add(v)
        
        
    def get_triples(self, id_only=False):
        """Retrieve a list of all triples of the tree."""
        
        if id_only:
            return [(a.ID, b.ID, c.ID) for a, b, c in self._triple_generator()]
        else:
            return [t for t in self._triple_generator()]
    
    
    def _triple_generator(self):
        
        self.supply_leaves()
        
        for u in self.preorder():
            for v1, v2 in itertools.permutations(u.children, 2):
                if len(v2.leaves) > 1:
                    for c in v1.leaves:
                        for a, b in itertools.combinations(v2.leaves, 2):
                            yield a, b, c
    
    
    def delete_and_reconnect(self, node):
        """Delete a node from the tree and reconnect its parent and children."""
        
        parent = node.parent
        if not parent:
            print("cannot delete and reconnect root '{}'".format(node))
            return False
        else:
            parent.remove_child(node)
            
            # copy list of children to edit edges
            children = [child for child in node.children]
            for child in children:
                parent.add_child(child)
                    
            node.children.clear()
        
        return parent
    
    
    def get_max_ID(self):
        """Returns the maximum of all node IDs."""
        
        max_ID = -1
        
        for v in self.preorder():
            if v.ID is not None:
                max_ID = max(v.ID, max_ID)
                
        return max_ID
        
    
    def to_newick(self, node=None):
        """Recursive Tree --> Newick (str) function."""
        
        def _to_newick(node):
            if not node.children:
                return str(node)
            else:
                s = ''
                for child in node.children:
                    s += _to_newick(child) + ','
                return "({}){}".format(s[:-1], node)
        
        if self.root:
            return _to_newick(self.root) + ';'
        else:
            return ';'
        
    
    def random_leaves(self, proportion):
        """Return a random subset of the leaves."""
        
        if (not isinstance(proportion, (float, int)) or 
            proportion < 0 or proportion > 1):
            raise ValueError('needs a number 0 <= p <= 1')
        
        leaves = [v for v in self.leaves()]
        k = round(proportion * len(leaves))
        
        return random.sample(leaves, k)
      
    
    def get_hierarchy(self):
        """Return the hierarchy set on the leaf labels defined by the tree."""
        
        self.supply_leaves()
        
        hierarchy = set()
        
        for v in self.preorder():
            
            A = [leaf.label for leaf in v.leaves]
            A.sort()
            A = tuple(A)
            hierarchy.add(A)
            
        return hierarchy
    
    
    def compare_topology(self, other):
        """Compare the tree topology based on the hierarchies.
        
        Only works for phylogenetic trees."""
        
        hierarchy1 = sorted(self.get_hierarchy())
        hierarchy2 = sorted(other.get_hierarchy())
        
        if len(hierarchy1) != len(hierarchy2):
            print('Unequal sizes of the hierarchy sets: '\
                  '{} and {}'.format(len(hierarchy1), len(hierarchy2)))
            return False
        
        for i in range(len(hierarchy1)):
            
            if hierarchy1[i] != hierarchy2[i]:
                print('Hierarchies not equal:'\
                      '\n{}\n{}'.format(hierarchy1[i], hierarchy2[i]))
                return False
        
        return True
    
    
    def _assert_integrity(self):
        
        for v in self.preorder():
            for child in v.children:
                if child is v:
                    raise RuntimeError('loop at {}'.format(v))
                if child.parent is not v:
                    raise RuntimeError('Tree invalid for '\
                                       '{} and {}'.format(v, child))
        
        return True
    
    
    @staticmethod
    def random_tree(N, binary=False):
        """Creates a random tree.
        
        The number of leaves is specified by the parameter 'N'. Each non-leaf
        node in the resulting tree will have at least children (property of
        phylogenetic trees).
        
        Keyword arguments:
            binary - forces the tree to be binary; default is False
        """
        
        if not (isinstance(N, int)) or N < 1:
            raise TypeError("N must be an 'int' > 0")
        root = TreeNode(0, label='0')
        tree = Tree(root)
        node_list = [root]
        nr, leaf_count = 1, 1
        
        while leaf_count < N:
            node = random.choice(node_list)
            
            if not node.children: 
                # to be phylogenetic at least two children must be added
                new_child1 = TreeNode(nr, label=str(nr))
                new_child2 = TreeNode(nr+1, label=str(nr+1))
                node.add_child(new_child1)
                node.add_child(new_child2)
                node_list.extend(node.children)
                nr += 2
                leaf_count += 1
            elif node.children and not binary:
                # add only one child if there are already children
                new_child = TreeNode(nr, label=str(nr))
                node.add_child(new_child)
                node_list.append(new_child)
                nr += 1
                leaf_count += 1
                
        return tree
    
# ----------------------------------------------------------------------------
#                       Efficient lca computation
#                               based on
#
# - Bender, M. A., M. Farach-Colton, G. Pemmasani, S. Skiena, and P. Sumazin:
#   Lowest common ancestors in trees and directed acyclic graphs.
#   In: Journal of Algorithms. 57, Nr. 2, November 2005, S. 75–94.
#   ISSN 0196-6774. doi:10.1016/j.jalgor.2005.08.001.
# - https://cp-algorithms.com/data_structures/sparse-table.html
# ----------------------------------------------------------------------------

class LCA:
    """Compute last common ancestors efficiently.
    
    Uses a reduction to the Range minimum query (RMQ) problem and a sparse
    table implementation.
    Preprocessing complexity: O(n * log n)
    Query complexity: O(1)
    where n is the number of vertices in the tree.
    """
    
    def __init__(self, tree):
        
        if not isinstance(tree, Tree):
            raise TypeError("tree must be of type 'Tree'")
            
        self.tree = tree
        
        self.V = [v for v in self.tree.preorder()]
        self.index = {v: i for i, v in enumerate(self.V)}
        
        # store IDs for queries via ID
        self.id_dict = {v.ID: v for v in self.V}
        
        self.euler_tour = []
        # levels of the vertices in the Euler tour
        self.L = []
        
        for v, level in self.tree.euler_and_level():
            self.euler_tour.append(self.index[v])
            self.L.append(level)
        
        # repres. of vertices in the Euler tour (index of first occurence)
        self.R = [None for _ in range(len(self.V))]
        for j, i in enumerate(self.euler_tour):
            if self.R[i] is None:
                self.R[i] = j
                
        # build sparse table for range minimum query (RMQ)
        self._precompute_logs()
        self._RMQ_sparse_table()
        
        
    def __call__(self, a, b):
        """Return the last common ancestor of two nodes."""
        
        return self._get_lca(self._id_to_treenode(a),
                             self._id_to_treenode(b))
        
        
    
    def get(self, a, b):
        """Return the last common ancestor of two nodes."""
        
        return self._get_lca(self._id_to_treenode(a),
                             self._id_to_treenode(b))
    
    
    def displays_triple(self, a, b, c):
        """Return whether the tree displays the rooted triple ab|c (= ba|c)."""
        
        try:
            return self._has_triple(self._id_to_treenode(a),
                                    self._id_to_treenode(b),
                                    self._id_to_treenode(c))
        except KeyError:
            return False
        
    
    def are_comparable(self, u, v):
        """Return whether two nodes/edges are comparable w.r.t. ancestor
        relation."""
        
        return self._are_comparable(self._id_to_treenode(u),
                                    self._id_to_treenode(v))
    
    
    def ancestor_or_equal(self, u, v):
        """Return whether u is equal to or an ancestor of v."""
        
        return self._ancestor_or_equal(self._id_to_treenode(u),
                                       self._id_to_treenode(v))
    
    
    def ancestor_not_equal(self, u, v):
        """Return whether u is a strict ancestor of v."""
        
        u = self._id_to_treenode(u)
        v = self._id_to_treenode(v)
        
        return u != v and self._ancestor_or_equal(u, v)
    
    
    def descendant_or_equal(self, u, v):
        """Return whether u is equal to or a descendant of v."""
        
        return self.ancestor_or_equal(v, u)
    
    
    def descendant_not_equal(self, u, v):
        """Return whether u is a strict descendant of v."""
        
        return self.ancestor_not_equal(v, u)
    
    
    def consistent_triples(self, triples):
        """Return a list with the subset of 'triples' that are displayed."""
        
        return [t for t in triples if self.displays_triple(*t)]
    
    
    def consistent_triple_generator(self, triples):
        """Generator for the items in 'triples' that are displayed."""
        
        for t in triples:
            if self.displays_triple(*t):
                yield t
    
    
    def _precompute_logs(self):
        
        n = len(self.L)
        self.log = [0 for _ in range(n + 1)]
        for i in range(2, n + 1):
            self.log[i] = int(self.log[i//2]) + 1
        
        
    def _RMQ_sparse_table(self):
        
        n = len(self.L)
        K = self.log[n]
        
        # lookup table M
        self.M = [[0 for j in range(K + 1)] for i in range(n)]
        
        # initialize the intervals with length 1
        for i in range(n):
            self.M[i][0] = i
         
        # dynamic programming: compute values from smaller to bigger intervals  
        for j in range(1, K + 1):
            
            # compute minimum value for all intervals with size 2^j
            for i in range(n - (1 << j) + 1):
                
                if (self.L[ self.M[i][j - 1] ] <
                    self.L[ self.M[i + (1 << (j - 1))][j - 1] ]):
                    self.M[i][j] = self.M[i][j - 1]
                else:
                    self.M[i][j] = self.M[i + (1 << (j - 1))][j - 1]
    
    
    def _RMQ_query(self, i, j):  
        
        k = self.log[j - i + 1]
        if self.L[self.M[i][k]] < self.L[self.M[j - (1 << k) + 1][k]]:
            return self.M[i][k]
        else:
            return self.M[j - (1 << k) + 1][k]
        
    
    def _id_to_treenode(self, v):
        
        if isinstance(v, TreeNode):
            return v
        elif isinstance(v, (tuple, list)) and len(v) == 2:
            return (self._id_to_treenode(v[0]),
                    self._id_to_treenode(v[1]))
        else:
            return self.id_dict[v]
        
        
    def _get_lca(self, v1, v2):
        
        if v1 is v2:
            return v1
        
        r1 = self.R[self.index[v1]]
        r2 = self.R[self.index[v2]]
        if r1 > r2:
            r1, r2 = r2, r1
        return self.V[ self.euler_tour[self._RMQ_query(r1, r2)] ]
        
    
    def _has_triple(self, a, b, c):
        
        lca_ab = self._get_lca(a, b)
        return lca_ab is not self._get_lca(lca_ab, c)
    
    
    def _are_comparable(self, u, v):
        
        return self._ancestor_or_equal(u, v) or self._ancestor_or_equal(v, u)
    
    
    def _ancestor_or_equal(self, u, v):
        
        # both are nodes
        if isinstance(u, TreeNode) and isinstance(v, TreeNode):
            return u is self._get_lca(u, v)
        
        # u node, v edge
        elif isinstance(u, TreeNode) and isinstance(v, tuple):
            return u is self._get_lca(u, v[0])
        
        # u edge, v node
        elif isinstance(u, tuple) and isinstance(v, TreeNode):
            return u[1] is self._get_lca(u[1], v)
        
        # both are edges
        elif isinstance(u, tuple) and isinstance(v, tuple):
            return u[1] is v[1] or u[1] is self._get_lca(u[1], v[0])
        

def lcas_naive(tree):
    """Naive computation of the last common ancestor for all leaf pairs."""
    
    lcas = {}
    tree.supply_leaves()
    
    for u in tree.preorder():
        for v1, v2 in itertools.permutations(u.children, 2):
            for x, y in itertools.product(v1.leaves, v2.leaves):
                lcas[x, y] = u
    
    return lcas


if __name__ == '__main__':
    
    # import time
    # tree = Tree.random_tree(4000)
    # leaves = tree.supply_leaves()
    
    # start_time1 = time.time()
    # lca_fast = LCA(tree)
    # end_time1 = time.time()
    
    # start_time2 = time.time()
    # lca_naive = lcas_naive(tree)
    # end_time2 = time.time()
    
    # print(end_time1 - start_time1, end_time2 - start_time2)
    
    # for i in range(10):
    #     v1, v2 = random.sample(leaves, 2)
        
    #     print(v1, v2, lca_fast.get(v1, v2), lca_naive[v1, v2])
    
    tree = Tree.random_tree(10)
    print(tree.to_newick())
    leaves = tree.supply_leaves()
    lca = LCA(tree)
    
    for i in range(10):
        a, b, c = random.sample(leaves, 3)
        
        print(a, b, c, lca.displays_triple(a, b, c),
                       lca.displays_triple(a.ID, b.ID, c.ID))