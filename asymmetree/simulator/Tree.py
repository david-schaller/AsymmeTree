# -*- coding: utf-8 -*-

"""
Tree data structure for phylogentic trees.

Classes in this module:
    - TreeNode
    - Tree
"""

import collections, itertools, random, re
import networkx as nx

__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


class TreeNode():
    """Class 'TreeNode'.
    
    Components for class 'Tree'. Contains a list of children as well as a 
    reference to the parent node.
    """
    
    __slots__ = ('ID', 'label', 'color', 'dist', 'tstamp', 'transferred',
                 'parent', 'children', 'leaves')
    
    def __init__(self, ID, label="", color=None,
                 dist = 1.0, tstamp=None, transferred=0, parent = None):
        self.ID = ID
        self.label = label
        self.color = color
        self.dist = dist
        self.tstamp = tstamp
        self.transferred = transferred
        self.parent = parent
        self.children = []
    
    def __repr__(self):
        return "tn" + str(self.ID)
    
    
    def __str__(self):
        if isinstance(self.color, tuple):
            return "{}<{}-{}>:{}".format(self.label, *self.color, self.dist)
        elif self.color:
            return "{}<{}>:{}".format(self.label, self.color, self.dist)
        else:
            return "{}:{}".format(self.label, self.dist)    
    
    
    def __eq__(self, other):
        return self.ID == other.ID
    
    
    def __lt__(self, other):
        return self.ID < other.ID


    def __hash__(self):
        return hash(id(self))


class Tree:
    
    def __init__(self, root):
        if not isinstance(root, TreeNode):
            raise TypeError("root must be of type 'TreeNode'.")
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
    
    
    def sorted_nodes(self):
        """Return list of sorted nodes by timestamp."""
    
        def quicksort(L, l, r):
            if r <= l:
                return
            i, j = l, r
            piv = L[(l+r)//2]
            while i <= j:
                while L[i].tstamp > piv.tstamp:
                    i += 1
                while L[j].tstamp < piv.tstamp:
                    j -= 1
                if i <= j:
                    L[i], L[j] = L[j], L[i]
                    i += 1
                    j -= 1
            quicksort(L, l, j)
            quicksort(L, i, r)
        
        nodes = [node for node in self.preorder()]
        quicksort(nodes, 0, len(nodes)-1)
        return nodes
    
    
    def sorted_edges(self):
        """Return list of sorted edges (u,v) by timestamp of u."""
    
        def quicksort(L, l, r):
            if r <= l:
                return
            i, j = l, r
            piv = L[(l+r)//2]
            while i <= j:
                while (L[i][0].tstamp > piv[0].tstamp or
                       (L[i][0].tstamp == piv[0].tstamp and
                        L[i][1].tstamp > piv[1].tstamp)):
                    i += 1
                while (L[j][0].tstamp < piv[0].tstamp or
                       (L[j][0].tstamp == piv[0].tstamp and
                        L[j][1].tstamp < piv[1].tstamp)):
                    j -= 1
                if i <= j:
                    L[i], L[j] = L[j], L[i]
                    i += 1
                    j -= 1
            quicksort(L, l, j)
            quicksort(L, i, r)
        
        edges = [edge for edge in self.edges()]
        quicksort(edges, 0, len(edges)-1)
        return edges
    
    
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
    
    
    def delete_and_reconnect(self, node):
        """Delete a node from the tree and reconnect its parent and children."""
        parent = node.parent
        if not parent:
            if node is self.root and len(node.children) == 1:
                self.root = node.children[0]
                self.root.parent = None
                self.root.dist = 0.0
                node.children.clear()
                return self.root
            else:
                print("Cannot delete and reconnect node", node)
                return False
        else:
            parent.children.remove(node)
            parent.children.extend(node.children)
            for child in node.children:
                child.parent = parent
                child.dist += node.dist
                if node.transferred == 1:
                    child.transferred = 1
            node.parent = None
            node.children.clear()
            return parent
    
    
    def supply_leaves(self, node=None, exclude_losses=True):
        """Add the leaves to all nodes that are in the subtree of a specific node."""
        if not node:
            self.supply_leaves(self.root, exclude_losses=exclude_losses)
            return self.root.leaves
        else:
            node.leaves = []
            if not node.children:
                if exclude_losses and node.label != "*":
                    node.leaves.append(node)
                elif not exclude_losses:
                    node.leaves.append(node)
            else:
                for child in node.children:
                    node.leaves.extend(self.supply_leaves(child, exclude_losses=exclude_losses))
            return node.leaves
    
    
    def distances_from_root(self):
        leaf_distances = []
        distance_dict = {}
        for v in self.preorder():
            if not v.parent:
                distance_dict[v] = 0.0
            else:
                depth = distance_dict[v.parent] + v.dist
                distance_dict[v] = depth
            if not v.children:
                leaf_distances.append( (str(v), distance_dict[v]) )
        return distance_dict, leaf_distances
    
    
    def get_triples(self):
        """Retrieve a list of all triples of the tree.
        
        Warning: Algorithm works recursively and can produce extremely 
        large lists.
        """
        def all_triples(node, triples):
            for child1 in node.children:
                for child2 in node.children:
                    if child1 is not child2:
                        for t3 in child1.leaves:
                            if len(child2.leaves) > 1:
                                for t1, t2 in itertools.combinations(child2.leaves, 2):
                                    triples.append( (t1, t2, t3) )
                all_triples(child1, triples=triples)
                
        self.supply_leaves()
        triples = []
        all_triples(self.root, triples)
        return triples

# --------------------------------------------------------------------------
#                          TREE  <--->  NEWICK
# --------------------------------------------------------------------------
        
    def to_newick(self, node=None, distance_only=False):
        """Recursive Tree --> Newick (str) function."""
        if node is None:
            return self.to_newick(self.root, distance_only=distance_only) + ";"
        elif not node.children:
            if not distance_only:
                return str(node)
            else:
                return ":{}".format(node.dist)
        else:
            s = ''
            for child in node.children:
                s += self.to_newick(node=child, distance_only=distance_only) + ","
            if not distance_only:
                return "({}){}".format(s[:-1], node)
            else:
                return "({}):{}".format(s[:-1], node.dist)
    
    
    @staticmethod
    def parse_newick(newick):
        """Parses trees in Newick format into object of type 'Tree'."""
        
        label_col_dist_regex = re.compile(r"'?([a-zA-Z0-9_]*)'?<(.*)>:(-?[0-9]*\.?[0-9]*[Ee]?-?[0-9]+)")  # label<color>:distance
        label_col_regex = re.compile(r"'?([a-zA-Z0-9_]*)'?<(.*)>")                                        # label<color>
        label_dist_regex = re.compile(r"'?([a-zA-Z0-9_]*)'?:(-?[0-9]*\.?[0-9]*[Ee]?-?[0-9]+)")            # label:distance
        
        id_counter = 0
        
        def parse_subtree(subroot, subtree_string):
            """Recursive function to parse the subtrees."""
            nonlocal id_counter
            children = split_children(subtree_string)
            for child in children:
                node = TreeNode(id_counter, parent=subroot)
                id_counter += 1
                subroot.children.append(node)
                end = -1
                if child[0] == '(':                                 # the child has subtrees
                    end = child.rfind(')')
                    if end == -1:
                        raise ValueError("Invalid Newick-String!")
                    parse_subtree(node, child[1:end])               # recursive call 'parse_subtree'
                child = child[end+1:].strip()
                label_col_dist = label_col_dist_regex.match(child)
                if label_col_dist:                                  # CASE 1: label<color>:distance
                    node.label = label_col_dist.group(1)
                    node.color = label_col_dist.group(2)
                    node.dist = float(label_col_dist.group(3))
                else:
                    label_col = label_col_regex.match(child)
                    label_dist = label_dist_regex.match(child)
                    if label_col:                                   # CASE 2: label<color>
                        node.label = label_col.group(1)
                        node.color = label_col.group(2)
                    elif label_dist:                                # CASE 3: label:distance
                        node.label = label_dist.group(1)
                        node.dist = float(label_dist.group(2))
                    else:                                           # CASE 4: label
                        node.label = child
                if node.color and node.color.find('-') != -1:
                    split_color = node.color.split('-')
                    node.color = (split_color[0], split_color[1])
                        
        def split_children(child_string):
            """Splits a given string by all ',' that are not enclosed by parentheses."""
            stack = 0
            children = []
            current = ""
            for c in child_string:
                if (stack == 0) and c == ',':
                    children.append(current)
                    current = ""
                elif c == '(':
                    stack += 1
                    current += c
                elif c == ')':
                    if stack <= 0:
                        raise ValueError("Invalid Newick-String!")
                    stack -= 1
                    current += c
                else:
                    current += c
            children.append(current.strip())
            return children
        
        if not isinstance(newick, str):
            raise ValueError("Newick parser needs a 'str' as input.")
        end = newick.find(";")
        if end != -1:
            newick = newick[:end]
        temp_root = TreeNode(-1)
        parse_subtree(temp_root, newick)
        if temp_root.children:
            root = temp_root.children[0]
            root.dist = 0.0                 # set distance of the root to 0
            root.parent = None              # remove the parent temp_root
                                            # (important for non-recursive to_newick2)
            return Tree(root)
        else:
            raise ValueError("Invalid Newick-String!")
    
    
    def reconstruct_IDs(self):
        """Reconstruct the (leaf) IDs."""
        self.number_of_species = 0
        IDs = set()
        
        for v in self.preorder():
            if not v.children:
                self.number_of_species += 1
            if v.label.isdigit():
                v.ID = int(v.label)
                IDs.add(v.ID)
        
        # assign new IDs to internal nodes
        current_ID = 0
        for v in self.preorder():
            if str(v.ID) != v.label:
                while current_ID in IDs:
                    current_ID += 1
                v.ID = current_ID
                IDs.add(current_ID)
                
                
    def reconstruct_timestamps(self):
        """Reconstruct the timestamps."""
        self.root.tstamp = 1.0
        for v in self.preorder():
            if v.parent:
                v.tstamp = v.parent.tstamp - v.dist
    
# --------------------------------------------------------------------------
#                         TREE  <--->  NETWORKX
# --------------------------------------------------------------------------
            
    def to_nx(self):
        self.check_integrity()
        G = nx.DiGraph()
        
        for u, v in self.edges():
            if u.ID not in G:
                G.add_node(u.ID, label=u.label, color=u.color, tstamp=u.tstamp)
            if v.ID not in G:
                G.add_node(v.ID, label=v.label, color=v.color, tstamp=v.tstamp)
            if u is v or u.ID == v.ID:
                print("Loop at", str(u), str(v), v.transferred)
            G.add_edge(u.ID, v.ID, dist=v.dist, transferred=v.transferred)
            
        return G, self.root.ID
    
    
    @staticmethod
    def parse_nx(G, root):
        
        number_of_leaves = 0
    
        def build_tree(ID, parent=None):
            nonlocal number_of_leaves
            treenode = TreeNode(ID, label=G.nodes[ID]['label'],
                                parent=parent)
            if parent:
                parent.children.append(treenode)
                if 'dist' in G.edges[parent.ID, ID]:
                    treenode.dist = G.edges[parent.ID, ID]['dist']
                if 'transferred' in G.edges[parent.ID, ID]:
                    treenode.transferred = G.edges[parent.ID, ID]['transferred']
            
            if 'label' in G.nodes[ID]:
                treenode.label = G.nodes[ID]['label']
            if 'color' in G.nodes[ID]:
                treenode.color = G.nodes[ID]['color']
            if 'tstamp' in G.nodes[ID]:
                treenode.tstamp = G.nodes[ID]['tstamp']
            
            for child_id in G.neighbors(ID):
                build_tree(child_id, parent=treenode)
            if G.out_degree(ID) == 0:
                number_of_leaves += 1
            
            return treenode
        tree = Tree(build_tree(root))
        tree.number_of_species = number_of_leaves
        
        return tree
# --------------------------------------------------------------------------

    
    @staticmethod
    def copy_tree(tree):
        
        orig_to_new = {}
        
        for orig in tree.preorder():
            new = TreeNode(orig.ID, label=orig.label, 
                           color=orig.color, dist=orig.dist,
                           tstamp=orig.tstamp, transferred=orig.transferred)
            orig_to_new[orig] = new
            if orig.parent:
                new.parent = orig_to_new[orig.parent]
                new.parent.children.append(new)
        
        return Tree(orig_to_new[tree.root])
    
    
    @staticmethod
    def random_colored_tree(N, colors):
        """Creates a random colored tree.
        
        The number of leaves and the color labels are specified in the
        parameters 'N' and 'colors', respectively. Each non-leaf node in the 
        resulting tree will have at least children (property of phylogenetic
        trees).
        """
        if not (isinstance(N, int) and isinstance(colors, collections.Iterable)):
            raise TypeError("N must be of type 'int' and colors must be iterable!")
        root = TreeNode(0, label="0")
        tree = Tree(root)
        node_list = [root]
        break_prob = [0.5, 0.5]
        nr, leaf_count = 1, 1
        while leaf_count < N:
            node = random.choice(node_list)
            if not node.children:                               # to be phylogenetic at least two children must be added
                new_child1 = TreeNode(nr, label=str(nr))
                new_child2 = TreeNode(nr+1, label=str(nr+1))
                new_child1.parent = node
                new_child2.parent = node
                node.children.append(new_child1)
                node.children.append(new_child2)
                node_list.extend(node.children)
                nr += 2
                leaf_count += 1
            elif node.children:                                 # add only one child if there are already children
                break_prob[0] = 1 / 1.4**len(node.children)     # probability to add another child (decraeses exponentially)
                break_prob[1] = 1 - break_prob[0]               # probability to choose another node
                if random.choices( (0,1), weights=break_prob )[0] == 1:
                    continue
                new_child = TreeNode(nr, label=str(nr))
                new_child.parent = node
                node.children.append(new_child)
                node_list.append(new_child)
                nr += 1
                leaf_count += 1
        for node in node_list:                          # assign colors
            if not node.children:                       # to leaves randomly
                node.color = random.choice(colors)
        return tree
    
    
    def check_integrity(self):
        for v in self.preorder():
            for child in v.children:
                if child is v:
                    print("Loop at " + str(v))
                    raise KeyboardInterrupt
                if child.parent is not v:
                    print("Tree invalid for " + str(v) + " and " + str(child))
                    raise KeyboardInterrupt

    
if __name__ == "__main__":
    colors = ("s", "t", "v", "w")
    N = 20
    
    t = Tree.random_colored_tree(N, colors)
    print("------------- Random tree test -------------")
    print( t.to_newick() )
    print("--------------------------------------------")
    
    t2 = Tree.parse_newick(t.to_newick())
    print( t2.to_newick() )
    
    nx_tree, nx_root = t.to_nx()
    t3 = Tree.parse_nx(nx_tree, nx_root)
    print("--------------------------------------------")
    print( t3.to_newick() )
