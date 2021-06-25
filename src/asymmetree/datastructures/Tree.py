# -*- coding: utf-8 -*-

import itertools, random

import asymmetree.datastructures.DoublyLinkedList as dll


__author__ = 'David Schaller'


class TreeNode:
    """Tree nodes for class Tree.
    
    Attributes
    ----------
    ID : int
        Integer identifier of the node, some methods use the convention -1 if
        an identifier is not needed.
    label: str
        Human interpretable label.
    parent: TreeNode
        Parent node of this node.
    children: dll.DLList
        Child nodes of this node in a doubly-linked list.
    leaves: list of TreeNode.
        A list of all leaf nodes in the subtree under this node, only available
        after execution of supply_leaves() function.
    
    See Also
    --------
    Tree
    PhyloTreeNode
    CotreeNode
    """
    
    __slots__ = ('ID', 'label', 'parent', '_par_dll_node',
                 'children', 'leaves')
    
    
    def __init__(self, ID, label=''):
        """Constructor for TreeNode class.
        
        Parameters
        ----------
        ID : int
            Integer identifier of the node, some methods use the convention -1
            if an identifier is not needed.
        label: str, optional
            Human interpretable label (the default is '').
        """
        
        self.ID = ID
        self.label = label
        
        self.parent = None
        # reference to doubly-linked list element in the parents' children
        self._par_dll_node = None
        
        self.children = dll.DLList()
        
    
    def __str__(self):
        
        return str(self.label)
    
    
    def __repr__(self):
        
        return '<TN: {}>'.format(self.ID)
    
                                            
    def add_child(self, child_node):
        """Add a node as a child of this node.
        
        Does nothing if the node is already a child node of this node.
        
        Parameters
        ----------
        child_node : TreeNode
            The node to add as a new child to this node.
        """
        
        # do nothing if child_node is already a child of self
        
        if child_node.parent is None:
            child_node.parent = self
            child_node._par_dll_node = self.children.append(child_node)
        
        elif child_node.parent is not self:
            child_node.parent.remove_child(child_node)
            child_node.parent = self
            child_node._par_dll_node = self.children.append(child_node)  
    
    
    def remove_child(self, child_node):
        """Remove a child node of this node.
        
        Parameters
        ----------
        child_node : TreeNode
            The node to be removed from the list of children.
            
        Raises
        ------
        KeyError
            If the supplied node is not a child of this node.
        """
        
        if child_node.parent is self:
            self.children.remove_node(child_node._par_dll_node)
            child_node.parent = None
            child_node._par_dll_node = None
        else:
            raise KeyError('{} is not a child of node {}'.format(child_node,
                                                                 self))
            
            
    def detach(self):
        """Detach this node from its parent.
        
        The node has no parent afterwards.
        """
        
        if self.parent is not None:
            self.parent.remove_child(self)
        else:
            self.parent = None
            self._par_dll_node = None
            
    
    def is_leaf(self):
        """Return True if the node is a leaf, False otherwise.
        
        Returns
        -------
        bool
            True if the node is a leaf, i.e. it has no children, else False.
        """
        
        return not self.children
        

class Tree:
    """Basic class for trees.
    
    Attributes
    ----------
    root : TreeNode
        The root node of the tree.
    """
    
    # corresponding node type
    node_type = TreeNode
    
    
    def __init__(self, root):
        """Constructor for the class tree.
        
        Parameters
        ----------
        root : TreeNode
            The root node for the newly created tree.
        """
        
        self.root = root
        
    
    def leaves(self):
        """Generator for leaves of the tree.
        
        Yields
        ------
        TreeNode
            The leaf nodes of the tree.
        """
        
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
        """Generator for preorder traversal of the tree.
        
        Yields
        ------
        TreeNode
            All nodes of the tree in pre-order.
        """
        
        def _preorder(node):
            yield node
            for child in node.children:
                yield from _preorder(child)
        
        if self.root:
            yield from _preorder(self.root)
        else:
            yield from []
            
    
    def traverse_subtree(self, u):
        """Generator for pre-order traversal of the subtree rooted at u.
        
        Yields
        ------
        TreeNode
            All nodes in the subtree rooted at node u in pre-order.
        """
        
        yield u
        for child in u.children:
            yield from self.traverse_subtree(child)
    
    
    def postorder(self):
        """Generator for post-order traversal of the tree.
        
        Yields
        ------
        TreeNode
            All nodes in the subtree rooted at node u in post-order.
        """
        
        def _postorder(node):
            for child in node.children:
                yield from _postorder(child)
            yield node
        
        if self.root:
            yield from _postorder(self.root)
        else:
            yield from []
            
    
    def inner_nodes(self):
        """Generator for inner nodes in pre-order.
        
        Yields
        ------
        TreeNode
            All inner nodes of the tree in pre-order.
        """
        
        def _inner_nodes(node):
            if node.children:
                yield node
                for child in node.children:
                    yield from _inner_nodes(child)
        
        if self.root:
            yield from _inner_nodes(self.root)
        else:
            yield from []
            
    
    def edges(self):
        """Generator for all edges of the tree.
        
        Yields
        ------
        tuple of two TreeNode objects
            All edges of the tree.
        """
        
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
        
        Yields
        ------
        tuple of two TreeNode objects and one int
            Edges uv as tuples (u, v, nr) where nr is the index of v in
            the list of children of node u.
        """
        
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
        """Generator for all inner edges of the tree.
        
        Yields
        ------
        tuple of two TreeNode objects
            All inner edges uv of the tree, i.e. edges for which the child v
            of u is not a leaf.
        """
        
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
        """Generator for an Euler tour of the tree.
        
        Parameters
        ----------
        id_only : bool
            If True, the tuples only contain the integer IDs of the nodes
            (the default is False).
        
        Yields
        ------
        TreeNode or int
            Nodes in an Euler tour of the tree.
        """
        
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
        """Generator for an Euler tour with node levels.
        
        Yields
        ------
        tuple of a TreeNode and an int
            Nodes and their level (distance from the root) in an Euler tour of
            the tree.
        """
        
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
        """Leaves in the subtree rooted at each node.
        
        Computes the list of leaves for every node in the tree containing the
        leaf nodes lying in the subtree rooted at the node.
        
        Returns
        -------
        list of TreeNode objects
            The leaves under the root, i.e. the complete list of leaves.
            Returns an empty list if the root is None.
        """
        
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
        """Contract edges in the tree.
        
        Parameters
        ----------
        edges : iterable object of tuples of two TreeNode objects
            The edges to be contracted in the tree.
        """
        
        contracted = set()
        
        for u, v in edges:
            
            # avoid trying to contract the same edge multiple times
            if v not in contracted:
                self.delete_and_reconnect(v)
            
            contracted.add(v)
        
        
    def get_triples(self, id_only=False):
        """Retrieve a list of all triples of the tree.
        
        A tree displays a triple ab|c on the leaf nodes a, b and c if the last
        common ancestor of a and b is a (proper) descendant of the last common
        ancestor of a and c (b and c).
        
        Parameters
        ----------
        id_only : bool
            If True, the triples are represented by the IDs of the nodes.
            
        Returns
        -------
        list of tuples of three TreeNode or int objects
            Each tuple (a, b, c) represents the triple ab|c (=ba|c), i.e. the
            first two items are closer related in the tree.
        """
        
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
        """Delete a node from the tree and reconnect its parent and children.
        
        Parameters
        ----------
        node : TreeNode
            The node to be deleted.
        
        Returns
        -------
        TreeNode or bool
            The parent of the node, if it could be deleted, or False, if the
            node could not be deleted, i.e., it has no parent.
        """
        
        parent = node.parent
        if not parent:
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
        """The maximum of all node IDs.
        
        Returns
        -------
        int
            The maximum of all node IDs in the tree, or -1 if all IDs are None
            or if there are no nodes.
        """
        
        max_ID = -1
        
        for v in self.preorder():
            if v.ID is not None:
                max_ID = max(v.ID, max_ID)
                
        return max_ID
        
    
    def to_newick(self, node=None):
        """Newick representation of the tree.
        
        Parameters
        ----------
        node : TreeNode, optional
            The node whose subtree shall be returned as a Newick string, the
            default is None, in which case the whole tree is returned in Newick
            format.
        
        Returns
        -------
        str
            A newick representation of the (sub)tree.
        """
        
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
        """A random sample of the leaves.
        
        Parameters
        ----------
        proportion : float
            The proportion of the sample w.r.t. the full set of leaves.
        
        Returns
        -------
        list of TreeNode objects
            A random sample of the leaves of the tree.
        
        Raises
        ------
        ValueError
            If `proportion` is not a number between 0 and 1. 
        """
        
        if (not isinstance(proportion, (float, int)) or 
            proportion < 0 or proportion > 1):
            raise ValueError('needs a number 0 <= p <= 1')
        
        leaves = [v for v in self.leaves()]
        k = round(proportion * len(leaves))
        
        return random.sample(leaves, k)
      
    
    def get_hierarchy(self):
        """Hierarchy set on the leaf labels defined by the tree.
        
        Every (phylogenetic) tree can be represented by a hierarchy on the set
        of its leaves.
        The labels of the leaf nodes must be unique.
        
        Returns
        -------
        set of lists of str objects
            Representing the hierarchy where the leaves are represented by
            their labels.
        """
        
        self.supply_leaves()
        
        hierarchy = set()
        
        for v in self.preorder():
            
            A = [leaf.label for leaf in v.leaves]
            A.sort()
            A = tuple(A)
            hierarchy.add(A)
            
        return hierarchy
    
    
    def compare_topology(self, other):
        """Compare the tree topology based on the leaf labels.
        
        Only works for phylogenetic trees with unique leaf labels.
        
        Parameters
        ----------
        other : Tree
            The tree which this tree is compared to.
        
        Returns
        -------
        bool
            True if the topologies are equal, else False.
        """
        
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
    
    
    def copy(self):
        """Return a copy of the tree.
        
        Constructs a deep copy of the tree, i.e. to the level of nodes.
        By default, the node attributes are all immutable data types.
        Hence, the original tree is not affected by operations on the copy.
        
        Returns
        -------
        Tree
            A copy of the tree.
        """
        
        if not self.root:
            return Tree(None)
        
        orig_to_new = {}
        
        for orig in self.preorder():
            new = TreeNode(orig.ID, label=orig.label)
            orig_to_new[orig] = new
            if orig.parent:
                orig_to_new[orig.parent].add_child(new)
        
        return Tree(orig_to_new[self.root])
    
    
    @staticmethod
    def random_tree(N, binary=False):
        """A random tree.
        
        The resulting tree is always phylogenetic, i.e., each inner node has
        at least two children.
        
        Parameters
        ----------
        N : int
            The desired number of leaves.
        binary : bool
            If True, the resulting tree is binary, otherwise it may contain
            multifurcations.
        
        Returns
        -------
        Tree
            A randomly generated tree with `N` leaves.
            
        Raises
        ------
        TypeError
            If `N` is not an integer >= 1.
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
# ----------------------------------------------------------------------------

class LCA:
    """Compute last common ancestors in a tree efficiently.
    
    Uses a reduction to the Range minimum query (RMQ) problem and a sparse
    table implementation.
    Preprocessing complexity: O(n * log n)
    Query complexity: O(1)
    where n is the number of nodes in the tree.
    
    References
    ----------
    .. [1] M. A. Bender, M. Farach-Colton, G. Pemmasani, S. Skiena, P. Sumazin.
       Lowest common ancestors in trees and directed acyclic graphs.
       In: Journal of Algorithms. 57, Nr. 2, November 2005, S. 75–94.
       ISSN 0196-6774. doi:10.1016/j.jalgor.2005.08.001.
    .. [2] https://cp-algorithms.com/data_structures/sparse-table.html
    """
    
    def __init__(self, tree):
        """Constructor for class LCA.
        
        Parameters
        ----------
        tree : Tree
            The Tree instance for which this instance will allow efficient
            last common ancestor queries.
        
        Raises
        ------
        TypeError
            If `tree` is not a Tree instance.
        """
        
        if not isinstance(tree, Tree):
            raise TypeError("tree must be of type 'Tree'")
            
        self._tree = tree
        
        self._V = [v for v in self._tree.preorder()]
        self._index = {v: i for i, v in enumerate(self._V)}
        
        # store IDs for queries via ID
        self._id_dict = {v.ID: v for v in self._V}
        
        self._euler_tour = []
        # levels of the vertices in the Euler tour
        self._L = []
        
        for v, level in self._tree.euler_and_level():
            self._euler_tour.append(self._index[v])
            self._L.append(level)
        
        # repres. of vertices in the Euler tour (index of first occurence)
        self._R = [None for _ in range(len(self._V))]
        for j, i in enumerate(self._euler_tour):
            if self._R[i] is None:
                self._R[i] = j
                
        # build sparse table for range minimum query (RMQ)
        self._precompute_logs()
        self._RMQ_sparse_table()
        
        
    def __call__(self, a, b):
        """Last common ancestor of two nodes.
        
        Parameters
        ----------
        a : TreeNode or int
            A node or its ID in the tree corresponding to this LCA instance.
        b : TreeNode or int
            A node or its ID in the tree corresponding to this LCA instance.
        
        Returns
        -------
        TreeNode
            The last common ancestor of `a` and `b`.
        """
        
        return self._get_lca(self._id_to_treenode(a),
                             self._id_to_treenode(b))
        
        
    
    def get(self, a, b):
        """Last common ancestor of two nodes.
        
        Parameters
        ----------
        a : TreeNode or int
            A node or its ID in the tree corresponding to this LCA instance.
        b : TreeNode or int
            A node or its ID in the tree corresponding to this LCA instance.
        
        Returns
        -------
        TreeNode
            The last common ancestor of `a` and `b`.
        """
        
        return self._get_lca(self._id_to_treenode(a),
                             self._id_to_treenode(b))
    
    
    def displays_triple(self, a, b, c):
        """Return whether the tree displays the rooted triple ab|c (= ba|c).
        
        Parameters
        ----------
        a : TreeNode or int
            A node or its ID in the tree corresponding to this LCA instance.
        b : TreeNode or int
            A node or its ID in the tree corresponding to this LCA instance.
        c : TreeNode or int
            A node or its ID in the tree corresponding to this LCA instance.
        
        Returns
        -------
        bool
            True if the tree displays the triple ab|c (= ba|c).
        """
        
        try:
            return self._has_triple(self._id_to_treenode(a),
                                    self._id_to_treenode(b),
                                    self._id_to_treenode(c))
        except KeyError:
            return False
        
    
    def are_comparable(self, u, v):
        """Returns True if two nodes/edges are comparable in the tree.
        
        Two nodes/edges are comparable if one lies on the unique path from the
        other to the root of the tree.
        
        Parameters
        ----------
        u : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
        v : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
            
        Return
        ------
        bool
            True if `u` and `v` are comparable in the tree, else False.
        """
        
        return self._are_comparable(self._id_to_treenode(u),
                                    self._id_to_treenode(v))
    
    
    def ancestor_or_equal(self, u, v):
        """Return True if u is equal to or an ancestor of v.
        
        Parameters
        ----------
        u : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
        v : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
            
        Return
        ------
        bool
            True if `u` is equal or an ancestor of `v`, else False.
        """
        
        return self._ancestor_or_equal(self._id_to_treenode(u),
                                       self._id_to_treenode(v))
    
    
    def ancestor_not_equal(self, u, v):
        """Return True if u is a strict ancestor of v.
        
        Parameters
        ----------
        u : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
        v : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
            
        Return
        ------
        bool
            True if `u` is a strict ancestor of `v`, else False.
        """
        
        u = self._id_to_treenode(u)
        v = self._id_to_treenode(v)
        
        return u != v and self._ancestor_or_equal(u, v)
    
    
    def descendant_or_equal(self, u, v):
        """Return True if u is equal to or a descendant of v.
        
        Parameters
        ----------
        u : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
        v : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
            
        Return
        ------
        bool
            True if `u` is equal or a descendant of `v`, else False.
        """
        
        return self.ancestor_or_equal(v, u)
    
    
    def descendant_not_equal(self, u, v):
        """Return True if u is a strict descendant of v.
        
        Parameters
        ----------
        u : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
        v : TreeNode or int or tuple of two TreeNode or int objects
            An node or edge in the tree corresponding to this LCA instance.
            
        Return
        ------
        bool
            True if `u` is a strict descendant of `v`, else False.
        """
        
        return self.ancestor_not_equal(v, u)
    
    
    def consistent_triples(self, triples):
        """List with the subset of triples that are displayed by the tree.
        
        Parameters
        ----------
        triples : an iterable object of tuples of three TreeNode or int objects
            Input triples of which each may or may not be displayed by the tree.
        
        Returns
        -------
        list of tuples of three TreeNode of int objects
            Representing the subset of the input list that are displayed by the
            tree.
        """
        
        return [t for t in triples if self.displays_triple(*t)]
    
    
    def consistent_triple_generator(self, triples):
        """Generator for the items in 'triples' that are displayed.
        
        Parameters
        ----------
        triples : an iterable object of tuples of three TreeNode or int objects
            Input triples of which each may or may not be displayed by the tree.
        
        Yields
        -------
        tuple of three TreeNode of int objects
            For each triple in the input list that is displayed by the tree.
        """
        
        for t in triples:
            if self.displays_triple(*t):
                yield t
    
    
    def _precompute_logs(self):
        
        n = len(self._L)
        self.log = [0 for _ in range(n + 1)]
        for i in range(2, n + 1):
            self.log[i] = int(self.log[i//2]) + 1
        
        
    def _RMQ_sparse_table(self):
        
        n = len(self._L)
        K = self.log[n]
        
        # lookup table M
        self._M = [[0 for j in range(K + 1)] for i in range(n)]
        
        # initialize the intervals with length 1
        for i in range(n):
            self._M[i][0] = i
         
        # dynamic programming: compute values from smaller to bigger intervals  
        for j in range(1, K + 1):
            
            # compute minimum value for all intervals with size 2^j
            for i in range(n - (1 << j) + 1):
                
                if (self._L[ self._M[i][j - 1] ] <
                    self._L[ self._M[i + (1 << (j - 1))][j - 1] ]):
                    self._M[i][j] = self._M[i][j - 1]
                else:
                    self._M[i][j] = self._M[i + (1 << (j - 1))][j - 1]
    
    
    def _RMQ_query(self, i, j):  
        
        k = self.log[j - i + 1]
        if self._L[self._M[i][k]] < self._L[self._M[j - (1 << k) + 1][k]]:
            return self._M[i][k]
        else:
            return self._M[j - (1 << k) + 1][k]
        
    
    def _id_to_treenode(self, v):
        
        if isinstance(v, TreeNode):
            return v
        elif isinstance(v, (tuple, list)) and len(v) == 2:
            return (self._id_to_treenode(v[0]),
                    self._id_to_treenode(v[1]))
        else:
            return self._id_dict[v]
        
        
    def _get_lca(self, v1, v2):
        
        if v1 is v2:
            return v1
        
        r1 = self._R[self._index[v1]]
        r2 = self._R[self._index[v2]]
        if r1 > r2:
            r1, r2 = r2, r1
        return self._V[ self._euler_tour[self._RMQ_query(r1, r2)] ]
        
    
    def _has_triple(self, a, b, c):
        
        if a is b:
            return False
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
    """Naive computation of the last common ancestor for all leaf pairs.
    
    .. deprecated:: 1.0.0
        `lcas_naive` may be removed in future version of AsymmeTree.
        
    Parameters
    ----------
    tree : Tree
        The tree for which the last common ancestors are computed.
    
    Returns
    -------
    dict
        Keys are all pairs of TreeNode objects that are leaves, and the values
        correspond to the TreeNode that is their last common ancestor.
    """
    
    lcas = {}
    tree.supply_leaves()
    
    for u in tree.preorder():
        for v1, v2 in itertools.permutations(u.children, 2):
            for x, y in itertools.product(v1.leaves, v2.leaves):
                lcas[x, y] = u
    
    return lcas
