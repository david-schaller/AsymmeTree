# -*- coding: utf-8 -*-

"""
AVL-tree implementation.
    
Balanced binary search tree.
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
    
    
    def left_size(self):
        
        return self.left.size if self.left else 0
    
    
    def right_size(self):
        
        return self.right.size if self.right else 0
    
    
    def __str__(self):
        return '<AVLTreeNode: {}>'.format(self.key)
    
    
    def copy(self):
        
        copy = AVLTreeNode(self.key, self.value)
        copy.height = self.height
        copy.size = self.size
        copy.balance = self.balance
        return copy
 

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
    
    
    def __getitem__(self, item):
        
        node = self._find(key)
        
        if not node:
            raise KeyError(str(key))
            
        return node.value
    
    
    def get(self, item, default=None):
        
        node = self._find(key)
        
        if node:
            return node.value
        else:
            return default
        
    
    def keys(self):
        """Return an iterator for the keys."""
        
        return AVLTreeIterator(self, mode=1)
    
    
    def values(self):
        """Return an iterator for the values."""
        
        return AVLTreeIterator(self, mode=2)
    
    
    def items(self):
        """Return an iterator for (key, value) pairs."""
        
        return AVLTreeIterator(self, mode=3)
    
    
    def insert(self, key, value=None):
        """Insert a key (and value) into the AVL tree."""
        
        if not self.root:
            self.root = AVLTreeNode(key, value)
        else:
            node = self._find_insert(key)
            
            if key < node.key:
                node.left = AVLTreeNode(key, value)
                node.left.parent = node
                self.root = self._rebalance(node)
            elif key > node.key:
                node.right = AVLTreeNode(key, value)
                node.right.parent = node
                self.root = self._rebalance(node)
                
                
    def remove(self, key):
        """Remove a key from the tree.
        
         Raises a KeyError if key is not in the tree."""
        
        node = self._find(key)
        
        if not node:
            raise KeyError(str(key))
            
        self._delete_node(node)
    
    
    def discard(self, key):
        """Remove a key from the tree if present."""
        
        node = self._find(key)
        
        if node:
            self._delete_node(node)
        
    
    def pop(self, key):
        """Remove a key from the tree and return its value.
        
         Raises a KeyError if key is not in the tree."""
        
        node = self._find(key)
        
        if not node:
            raise KeyError(str(key))
            
        self._delete_node(node)
        return node.value
        
    
    def clear(self):
        """Removes all items from the tree."""
        
        self.root = None
    
    
    def key_at(self, index):
        """Return the key at the index."""
        
        return self._node_at(index).key
    
    
    def value_at(self, index):
        """Return the value at the index."""
        
        return self._node_at(index).value
    
    
    def key_and_value_at(self, index):
        """Return the (key, value) pair at the index."""
        
        node = self._node_at(index)
        return (node.key, node.value)
    
    
    def remove_at(self, index):
        """Remove node at the index."""
        
        self._delete_node(self._node_at(index))
    
    
    def pop_at(self, index):
        """Remove node at the index and return its (key, value) pair."""
        
        node = self._node_at(index)
        self._delete_node(node)
        return (node.key, node.value)
    
    
    def _node_at(self, index):
        """Return the node at the index."""
        
        if index < 0:
            if index < -self.root.size:
                raise IndexError('index {} is out of range'.format(index))
            else:
                index += self.root.size
        
        if index >= self.root.size:
            raise IndexError('index {} is out of range'.format(index))
        
        current = self.root
        current_sum = 0
        
        while current:
            current_index = current_sum + current.left_size()
            if index == current_index:
                return current
            elif index < current_index:
                current = current.left
            else:
                current = current.right
                current_sum = current_index + 1
                
        raise RuntimeError('could not find node with index {}'.format(index))
    
    
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
                
                
    def _find_insert(self, key):
        
        current = self.root
        while True:
            if key < current.key and current.left:
                current = current.left
            elif key > current.key and current.right:
                current = current.right
            else:
                return current
    
    
    def _delete_node(self, node):
        
        if node.left and node.right:
            
            # replace by smallest in right subtree
            subst = self._smallest_in_subtree(node.right)
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
            self.root = self._rebalance(to_rebalance)
            
        elif node.left:
            par = node.parent
            node.left.parent = par
            if par:
                if node is par.left:
                    par.left = node.left
                elif node is par.right:
                    par.right = node.left
                self.root = self._rebalance(par)
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
                self.root = self._rebalance(par)
            else:
                self.root = node.right
                
        else:
            par = node.parent
            if par:
                if node is par.left:
                    par.left = None
                elif node is par.right:
                    par.right = None
                self.root = self._rebalance(par)
            else:
                self.root = None
        
        node.parent, node.left, node.right = None, None, None
    
    
    def _rebalance(self, node):
        
        while node:
            balance = node.update()
            while abs(balance) > 1:
                if balance > 1:
                    if node.left.balance >= 0:
                        # single rotation needed
                        self._rightrotate(node)
                    else:      
                        # double rotation needed
                        self._leftrotate(node.left)
                        self._rightrotate(node)
                else:
                    if node.right.balance <= 0:
                        # single rotation needed
                        self._leftrotate(node)
                    else:
                        # double rotation needed
                        self._rightrotate(node.right)
                        self._leftrotate(node)
                balance = node.balance
            if not node.parent:
                return node
            node = node.parent
            
            
    def _rightrotate(self, y):
        
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
        
        
    def _leftrotate(self, x):
        
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
    
    
    def _smallest_in_subtree(self, node):
        
        current = node
        while current.left:
            current = current.left
        return current
    
    
    def _biggest_in_subtree(self, node):
        
        current = node
        while current.right:
            current = current.right
        return current
    
    
    def copy(self):
        """Return a copy of the tree."""
        
        def _copy(node, parent=None):
            
            node_copy = node.copy()
            node_copy.parent = parent
            if node.left:
                node_copy.left = _copy(node.left, parent=node_copy)
            if node.right:
                node_copy.right = _copy(node.right, parent=node_copy)
            return node_copy
            
        
        tree_copy = AVLTree()
        if self.root:
            tree_copy.root = _copy(self.root)
        return tree_copy
                
    
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
    
    def __init__(self, avl_tree, mode=1):
        
        self.avl_tree = avl_tree
        self._current = self.avl_tree.root
        
        # Where do I come from?
        # 1 -- up
        # 2 -- left
        # 3 -- right
        self._from = 1
        
        # What is returned?
        # 1 -- key
        # 2 -- value
        # 3 -- (key, value)
        self._mode = mode
        
        
    def __iter__(self):
        
        return self
        
    
    def __next__(self):
        
        while self._current:
            
            # coming from above
            if self._from == 1:
                if self._current.left:
                    self._current = self._current.left
                else:
                    self._from = 2
            # coming from left child --> return this node
            elif self._from == 2:
                x = self._current
                
                if self._current.right:
                    self._current = self._current.right
                    self._from = 1
                else:
                    self._current = self._current.parent
                    if self._current and self._current.left is x:
                        self._from = 2
                    elif self._current:
                        self._from = 3
                
                if self._mode == 1:
                    return x.key
                elif self._mode == 2:
                    return x.value
                else:
                    return (x.key, x.value)
                
            # coming from right child
            else:
                x = self._current
                self._current = self._current.parent
                if self._current and self._current.left is x:
                    self._from = 2
                elif self._current:
                    self._from = 3
        
        # stop when current node becomes None
        raise StopIteration
        

if __name__ == '__main__':
    
    import random
    
    keys = [i for i in range(10000)]
    random.shuffle(keys)
    
    t = AVLTree()
    for key in keys:
        t.insert(key)
    # print(t.to_newick())  
     
    # t = set()
    # for key in keys:
    #     t.add(key)
    
    t = t.copy()
    print(len(t))
    print(t.pop_at(-10000))
    print(len(t))
    print(t.key_and_value_at(980))
    
    l = [key for key in t.items()]
    print(l[-5:])
