# -*- coding: utf-8 -*-


import tools.DoublyLinkedList as dll


class SimpleTreeNode:
    
    __slots__ = ['ID', 'label', 'parent', 'parent_dll_element',
                 'children', 'leaves']
    
    
    def __init__(self, ID, label=None):
        
        self.ID = ID
        self.label = label
        self.parent = None
        self.children = dll.DLList()
        self.parent_dll_element = None      # reference to doubly-linked list element
                                            # in the parents' children
                                            
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
        

class SimpleTree:
    
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
        