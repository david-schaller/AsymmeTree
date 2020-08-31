# -*- coding: utf-8 -*-

"""
AVL-tree implementation.
    
Balanced binary tree.
"""


__author__ = 'David Schaller'


class AVLTreeNode:
    
    __slots__ = ('key', 'value',
                 'parent', 'left', 'right',
                 'height', 'size', 'balance')
    
    def __init__(self, key, value):
        
        self.key = key
        self.value = value
        
        self.parent = None
        self.left = None
        self.right = None
        
        self.height = 1     # height of the subtree
        self.size = 1       # stores number of elements in its subtree
        self.balance = 0
   
    
    def update(self):
        """Update height, size and balance of (the subtree under) the node."""
        
        if self.left and self.right:
            self.height = 1 + max(self.left.height, self.right.height)
            self.size = 1 + self.left.size + self.right.size
            self.balance = self.left.height - self.right.height
        elif self.left:
            self.height = 1 + self.left.height
            self.size = 1 + self.left.size
            self.balance = self.left.height
        elif self.right:
            self.height = 1 + self.right.height
            self.size = 1 + self.right.size
            self.balance = -(self.right.height)
        else:
            self.height = 1
            self.size = 1
            self.balance = 0
        
        return self.balance
    
    
    def __str__(self):
        return '<AVLTreeNode: {}>'.format(self.key)
 

class AVLTree:
    """AVL tree."""
    
    __slots__ = ('root')
    
    def __init__(self, root=None):
        
        self.root = root
    
    
    def __iter__(self):
        
        return AVLTreeIterator(self)


    def __next__(self):
        
        pass
            
            
    def __nonzero__(self):
        
        return True if self.root else False
      
    
    def __len__(self):
        
        return self.root.size if self.root else 0
    
    
    def __contains__(self, item):
        
        return self._find(item) is not None
    
    
    def _find(self, key):
        
        if not self.root:
            return None
        
        current = self.root
        while current:
            if key == current.key:
                return current
            elif key < current.key:
                current = current.left
            else:
                current = current.right
    
    
    def insert(self, key, value=None):
        """Insert a key (and value) into the AVL tree."""
        
        def _insert(node):
            if key == node.key:
                return
            elif key < node.key:
                if node.left:
                    _insert(key, node.left)
                else:
                    node.left = AVLTreeNode(key, value)
                    node.left.parent = node
                    self.root = self.rebalance(node)
            else:
                if node.right:
                    _insert(key, node.right)
                else:
                    node.right = AVLTreeNode(key, value)
                    node.right.parent = node
                    self.root = self.rebalance(node)
        
        if not self.root:
            self.root = AVLTreeNode(key, value)
            return
        
        self._insert(key, self.root)
    
    
    def delete_node(self, node):
        
        if node.left and node.right:
            
            # replace by smallest in right subtree
            subst = self.smallest_in_subtree(node.right)
            to_rebalance = subst.parent
            if to_rebalance is node:
                to_rebalance = subst
            else:
                to_rebalance.left = subst.right
                if subst.right:
                    subst.right.parent = to_rebalance
                subst.right = node.right
                node.right.parent = subst
            subst.left = node.left
            node.left.parent = subst
            par = node.parent
            subst.parent = par
            if par:
                if node is par.left:
                    par.left = subst
                elif node is par.right:
                    par.right = subst
            else:
                self.root = subst
            self.root = self.rebalance(to_rebalance)
            
        elif node.left:
            par = node.parent
            node.left.parent = par
            if par:
                if node is par.left:
                    par.left = node.left
                elif node is par.right:
                    par.right = node.left
                self.root = self.rebalance(par)
            else:
                self.root = node.left
                
        elif node.right:
            par = node.parent
            node.right.parent = par
            if par:
                if node is par.left:
                    par.left = node.right
                elif node is par.right:
                    par.right = node.right
                self.root = self.rebalance(par)
            else:
                self.root = node.right
                
        else:
            par = node.parent
            if par:
                if node is par.left:
                    par.left = None
                elif node is par.right:
                    par.right = None
                self.root = self.rebalance(par)
            else:
                self.root = None
        
        node.parent, node.left, node.right = None, None, None
    
    
    def rebalance(self, node):
        
        while node:
            balance = node.update()
            while abs(balance) > 1:
                if balance > 1:
                    if node.left.balance >= 0:
                        # single rotation needed
                        self.rightrotate(node)
                    else:      
                        # double rotation needed
                        self.leftrotate(node.left)
                        self.rightrotate(node)
                else:
                    if node.right.balance <= 0:
                        # single rotation needed
                        self.leftrotate(node)
                    else:
                        # double rotation needed
                        self.rightrotate(node.right)
                        self.leftrotate(node)
                balance = node.balance
            if not node.parent:
                return node
            node = node.parent
            
            
    def rightrotate(self, y):
        
        x = y.left
        B = x.right
        if y.parent and (y is y.parent.right):
            y.parent.right = x
        elif y.parent and (y is y.parent.left):
            y.parent.left = x
        x.parent = y.parent
        y.parent = x
        x.right = y
        y.left = B
        if B:
            B.parent = y
        y.update()
        x.update()
        
        
    def leftrotate(self, x):
        
        y = x.right
        B = y.left
        if x.parent and (x is x.parent.right):
            x.parent.right = y
        elif x.parent and (x is x.parent.left):
            x.parent.left = y
        y.parent = x.parent
        x.parent = y
        y.left = x
        x.right = B
        if B:
            B.parent = x
        x.update()
        y.update()
    
    
    def smallest_in_subtree(self, node):
        
        current = node
        while current.left:
            current = current.left
        return current
    
    
    def biggest_in_subtree(self, node):
        
        current = node
        while current.right:
            current = current.right
        return current
                
    
    def to_newick(self):
        """Recursive Tree --> Newick (str) function.
        
        Intended for testing purpose.
        """
        
        def _newick(node):
            if not (node.left or node.right):
                return str(node.key)# + "-" + str(node.size)
            else:
                if node.left and node.right:
                    s = '({},{})'.format(_newick(node.left),
                                         _newick(node.right))
                elif node.left:
                    s = '({},-)'.format(_newick(node.left))
                elif node.right:
                    s = '(-,{})'.format(_newick(node.right))
                else:
                    s = ''
                return s + str(node.key)# + "-" + str(node.size)
            
        return _newick(self.root) if self.root else ''
    
    
    def check_integrity(self, node=None):
        """Recursive check of the following properties of AVL trees:
            
        - all children have a correct parent reference
        - the tree is balanced (AVL property)
        - the size is correct in all subtrees
        Intended for testing purpose.
        """
        
        if not node:
            if not self.root:
                print("Tree has no root!")
                return False
            else:
                return self.check_integrity(self.root)
        else:
            height_left, height_right, size = 0, 0, 1
            if node.left:
                if (node is not node.left.parent or 
                    not self.check_integrity(node.left)):
                    print("Check node: '{}' left".format(node))
                    return False
                height_left = node.left.height
                size += node.left.size
            if node.right:
                if (node is not node.right.parent or
                    not self.check_integrity(node.right)):
                    print("Check node:'{}' right".format(node))
                    return False
                height_right = node.right.height
                size += node.right.size
            if abs(node.balance) > 1:
                print("Node '{}' is unbalanced!".format(node))
                return False
            if node.height != 1 + max(height_left, height_right):
                print("Height of node '{}' is incorrect!".format(node))
                return False
            if node.size != size:
                print("Size of node '{}' is incorrect!".format(node))
                return False
            return True
        
        
class AVLTreeIterator:
    """Iterator class for AVL tree."""
    
    def __init__(self, avl_tree):
        
        ######  not yet implemented ######
        self.avl_tree = avl_tree
        self._current = None
        
    
    def __next__(self):
        
        if self._current:
            x = self._current
            self._current = self._current._next
            return x._value
        else:
            raise StopIteration
