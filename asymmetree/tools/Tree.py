# -*- coding: utf-8 -*-

import itertools

import tools.DoublyLinkedList as dll


class TreeNode:
    """Class 'TreeNode'.
    
    Components for class 'Tree'. Contains a list of children as well as a
    reference to the parent node.
    """
    
    __slots__ = ['ID', 'label', 'parent', 'parent_dll_element',
                 'children', 'leaves']
    
    
    def __init__(self, ID, label=None):
        
        self.ID = ID
        self.label = label
        self.parent = None
        self.children = dll.DLList()
        self.parent_dll_element = None      # reference to doubly-linked list element
                                            # in the parents' children
    
    
    def __repr__(self):
        return "tn" + str(self.ID)
    
                                            
    def add_child(self, child_node):
        
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
            self.children.remove_element(child_node.parent_dll_element)
            child_node.parent = None
            child_node.parent_dll_element = None
        else:
            raise ValueError("Not a child of this node!")
            
            
    def detach(self):
        
        if self.parent is not None:
            self.parent.remove_child(self)
        else:
            self.parent = None
            self.parent_dll_element = None
            
    
#    def __eq__(self, other):
#        return self.ID == other.ID
#    
#    
#    def __lt__(self, other):
#        return self.ID < other.ID
#
#
#    def __hash__(self):
#        return hash(id(self))
        

class Tree:
    
    def __init__(self, root):
        self.root = root
    
    
    def preorder(self, node=None):
        """Generator for preorder traversal of the tree."""
        if not node:
            yield from self.preorder(node=self.root)
        else:
            yield node
            for child in node.children:
                yield from self.preorder(node=child)
    
    
    def postorder(self, node=None):
        """Generator for postorder traversal of the tree."""
        if not node:
            yield from self.postorder(node=self.root)
        else:
            for child in node.children:
                yield from self.postorder(node=child)
            yield node
            
    
    def edges(self, node=None):
        
        """Generator for all edges of the tree."""
        if not node:
            yield from self.edges(node=self.root)
        else:
            for child in node.children:
                yield (node, child)
                yield from self.edges(node=child)
                
    
    def euler_generator(self, node=None, id_only=False):
        """Generator for an Euler tour of the tree."""
        if not node:
            yield from self.euler_generator(node=self.root, id_only=id_only)
        else:
            if id_only: yield node.ID
            else: yield node
            
            for child in node.children:
                yield from self.euler_generator(node=child, id_only=id_only)
                
                if id_only: yield node.ID
                else: yield node
            
    
    def supply_leaves(self, node=None):
        """Add the leaves to all nodes that are in the subtree of a specific node."""
        if not node:
            self.supply_leaves(self.root)
            return self.root.leaves
        else:
            node.leaves = []
            if not node.children:
                node.leaves.append(node)
            else:
                for child in node.children:
                    node.leaves.extend(self.supply_leaves(child))
            return node.leaves
        
        
    def get_triples(self):
        """Retrieve a list of all triples of the tree.
        
        Warning: Algorithm works recursively and can produce extremely 
        large lists.
        """
                
        self.supply_leaves()
        triples = []
        
        for u in self.preorder():
            for v1, v2 in itertools.permutations(u.children, 2):
                if len(v2.leaves) > 1:
                    for t3 in v1.leaves:
                        for t1, t2 in itertools.combinations(v2.leaves, 2):
                            triples.append( (t1, t2, t3) )
        
        return triples
        
    
    def to_newick(self, node=None):
        """Recursive Tree --> Newick (str) function."""
        if node is None:
            return self.to_newick(self.root) + ";"
        elif not node.children:
            return str(node)
        else:
            s = ''
            for child in node.children:
                s += self.to_newick(node=child) + ","
            return f"({s[:-1]}){node}"
        