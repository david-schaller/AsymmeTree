# -*- coding: utf-8 -*-

import itertools, random

import asymmetree.tools.DoublyLinkedList as dll


class TreeNode:
    """Class 'TreeNode'.
    
    Components for class 'Tree'. Contains a list of children as well as a
    reference to the parent node.
    """
    
    __slots__ = ['ID', 'label', 'parent', 'parent_dll_element',
                 'children', 'leaves']
    
    
    def __init__(self, ID, label=""):
        
        self.ID = ID
        self.label = label
        self.parent = None
        self.children = dll.DLList()
        self.parent_dll_element = None      # reference to doubly-linked list element
                                            # in the parents' children
    
    def __str__(self):
        
        return str(self.label)
    
    
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
        

class Tree:
    
    # corresponding node type
    node_type = TreeNode
    
    
    def __init__(self, root):
        self.root = root
    
    
    def preorder(self, node=None):
        """Generator for preorder traversal of the tree."""
        
        if node is None:
            yield from self.preorder(node=self.root)
        else:
            yield node
            for child in node.children:
                yield from self.preorder(node=child)
    
    
    def postorder(self, node=None):
        """Generator for postorder traversal of the tree."""
        
        if node is None:
            yield from self.postorder(node=self.root)
        else:
            for child in node.children:
                yield from self.postorder(node=child)
            yield node
            
    
    def inner_vertices(self, node=None):
        """Generator for inner vertices in preorder."""
        
        if node is None:
            yield from self.inner_vertices(node=self.root)
        else:
            if node.children:
                yield node
                for child in node.children:
                    yield from self.inner_vertices(node=child)
            
    
    def edges(self, node=None):
        """Generator for all edges of the tree."""
        
        if node is None:
            yield from self.edges(node=self.root)
        else:
            for child in node.children:
                yield (node, child)
                yield from self.edges(node=child)
                
    
    def inner_edges(self, node=None):
        """Generator for all inner edges of the tree."""
        
        if node is None:
            yield from self.inner_edges(node=self.root)
        else:
            for child in node.children:
                if child.children:
                    yield (node, child)
                    yield from self.inner_edges(node=child)
                
    
    def euler_generator(self, node=None, id_only=False):
        """Generator for an Euler tour of the tree."""
        
        if node is None:
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
        
        if node is None:
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
    
    
    def delete_and_reconnect(self, node):
        """Delete a node from the tree and reconnect its parent and children."""
        
        parent = node.parent
        if not parent:
            print("Cannot delete and reconnect root '{node}'!".format(node))
            return False
        else:
            parent.remove_child(node)
            
            # copy list of children to edit edges
            children = [child for child in node.children]
            for child in children:
                parent.add_child(child)
                    
            node.children.clear()
        
        return parent
        
    
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
            return "({}){}".format(s[:-1], node)
      
    
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
        
        Only works for binary trees."""
        
        hierarchy1 = sorted(self.get_hierarchy())
        hierarchy2 = sorted(other.get_hierarchy())
        
        if len(hierarchy1) != len(hierarchy2):
            print("Unequal sizes of the hierarchy sets: {} and {}".format(len(hierarchy1), len(hierarchy2)))
            return False
        
        for i in range(len(hierarchy1)):
            
#            if len(hierarchy1[i]) > 1:
#                print(hierarchy1[i])
#                print(hierarchy2[i])
            
            if hierarchy1[i] != hierarchy2[i]:
                print("Hierarchies not equal:\n{}\n{}".format(hierarchy1[i], hierarchy2[i]))
                return False
        
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
        
        if not (isinstance(N, int)):
            raise TypeError("N must be of type 'int'!")
        root = TreeNode(0, label="0")
        tree = Tree(root)
        node_list = [root]
        break_prob = [0.5, 0.5]
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
                
                # probability to add another child (decraeses exponentially)
                break_prob[0] = 1 / 1.4**len(node.children)
                # probability to choose another node
                break_prob[1] = 1 - break_prob[0]
                if random.choices( (0,1), weights=break_prob )[0] == 1:
                    continue
                
                new_child = TreeNode(nr, label=str(nr))
                node.add_child(new_child)
                node_list.append(new_child)
                nr += 1
                leaf_count += 1
                
        return tree