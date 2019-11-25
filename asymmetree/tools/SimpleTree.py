# -*- coding: utf-8 -*-


class SimpleTreeNode:
    
    __slots__ = ['ID', 'label', 'parent', 'children',
                 'leaves']
    
    
    def __init__(self, ID, label=None, parent=None):
        self.ID = ID
        self.label = label
        self.parent = parent
        self.children = []
        

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
            return "({}){}".format(s[:-1], str(node))
        