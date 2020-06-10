# -*- coding: utf-8 -*-

"""
Tree data structure for phylogentic trees.
"""

import collections, itertools, random, re, pickle, json, os
import numpy as np
import networkx as nx

from asymmetree.tools.Tree import Tree, TreeNode


__author__ = 'David Schaller'


class PhyloTreeNode(TreeNode):
    """Class 'PhyloTreeNode'.
    
    Components for class 'PhyloTree'.
    """
    
    __slots__ = ('color', 'dist', 'tstamp', 'transferred')
    
    
    def __init__(self, ID, label='', color=None,
                 dist = 1.0, tstamp=None, transferred=0):
        
        super().__init__(ID, label=label)
        self.color = color                      # color / reconciliation
        self.dist = dist                        # distance from parent
        self.tstamp = tstamp                    # time stamp
        self.transferred = transferred          # 1 if edge parent-self was an 
                                                # HGT edge, and 0 otherwise

    
    def __repr__(self):
        
        return '<PhyloTN: {}>'.format(self.ID)
    
    
    def __str__(self):
        
        if isinstance(self.color, (tuple, list)):
            return '{}<{}-{}>:{}'.format(self.label, *self.color, self.dist)
        elif self.color:
            return '{}<{}>:{}'.format(self.label, self.color, self.dist)
        else:
            return '{}:{}'.format(self.label, self.dist)
        
    
    def is_loss(self):
        
        return self.label == '*'


class PhyloTree(Tree):
    
    # corresponding node type
    node_type = PhyloTreeNode
    
    
    def __init__(self, root):
        
        if root is not None and not isinstance(root, PhyloTreeNode):
            raise TypeError("root must be of type 'PhyloTreeNode'")
        super().__init__(root)
    
    
    def sorted_nodes(self):
        """Return list of sorted nodes by timestamp."""
        
        return sorted(self.preorder(),
                      key=lambda node: node.tstamp,
                      reverse=True)
    
    
    def sorted_edges(self):
        """Return list of sorted edges (u,v) by timestamp of u."""
    
        return sorted(self.edges(),
                      key=lambda e: (e[0].tstamp, e[1].tstamp),
                      reverse=True)
    
    
    def delete_and_reconnect(self, node,
                             add_distances=True,
                             keep_transferred=True):
        """Delete a node from the tree and reconnect its parent and children."""
        
        distance, transferred = node.dist, node.transferred
        parent = super().delete_and_reconnect(node)
        
        if parent:
            for child in parent.children:
                if add_distances:
                    child.dist += distance
                if keep_transferred and transferred:
                    child.transferred = 1
        
        return parent
        
    
    def contract(self, edges):
        
        contracted = set()
        
        for u, v in edges:
            
            # avoid trying to contract the same edge multiple times
            if v not in contracted:
                self.delete_and_reconnect(v)
            
            contracted.add(v)
            
            
    def remove_planted_root(self):
        """Removes the planted root if existent.
        
        The the tree has a planted root, its single child becomes the new root.
        """
        
        if len(self.root.children) == 1:
        
            new_root = self.root.children[0]
            new_root.detach()
            self.root = new_root
            new_root.dist = 0.0
            
    
    def supply_leaves(self, exclude_losses=True):
        """Add the leaves to all nodes under their respective subtrees."""
        
        def _supply_leaves(node, exclude_losses):
        
            node.leaves = []
            
            if not node.children:
                if exclude_losses and not node.is_loss():
                    node.leaves.append(node)
                elif not exclude_losses:
                    node.leaves.append(node)
            else:
                for child in node.children:
                    node.leaves.extend(_supply_leaves(child, exclude_losses))
                    
            return node.leaves
        
        
        if self.root:
            return _supply_leaves(self.root, exclude_losses)
        else:
            return []
    
    
    def color_sorted_leaves(self):
        
        self.supply_leaves()
        color_dict = {}
        
        for leaf in self.root.leaves:
            if leaf.color not in color_dict:
                color_dict[leaf.color] = []
            color_dict[leaf.color].append(leaf)
            
        leaves = []
        for color, leaf_list in color_dict.items():
            for leaf in leaf_list:
                leaves.append(leaf)
        
        return leaves
    
    
    def distance_matrix(self, leaf_order=None):
        """Computes the distance matrix from the phylogenetic tree.
        
        Additionally a list of leaves corresponding to the indices in the matrix
        is returned.
        
        Keyword arguments:
            leaf_order -- list of leaves defining the indices for the matrix;
                default is None, in which case leaves are indexed in sibling order.
        """
        
        distance_dict, _ = self.distances_from_root()
        self.supply_leaves()
        
        if leaf_order:
            if set(leaf_order) != set(self.root.leaves):
                raise ValueError('ordered leaf list does not match with the '\
                                 'leaves in the tree')
            leaves = leaf_order
        else:
            # leaves in sibling order
            leaves = self.root.leaves      
        
        leaf_index = {l: i for i, l in enumerate(leaves)}
        
        D = np.zeros((len(leaves),len(leaves)), dtype=np.float)
        
        for v in self.preorder():
            if v.children:
                for c1, c2 in itertools.combinations(v.children, 2):
                    for x in c1.leaves:
                        x_index = leaf_index[x]
                        x_dist = distance_dict[x] - distance_dict[v]
                        for y in c2.leaves:
                            y_index = leaf_index[y]
                            y_dist = distance_dict[y] - distance_dict[v]
                            D[x_index, y_index] = x_dist + y_dist
                            D[y_index, x_index] = x_dist + y_dist
        
        return leaves, D
    
    
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
    
    
    def topology_only(self, inplace=True):
        """Reset distances, time stamps, transfer status and inner labels."""
        
        if not inplace:
            T = self.copy()
        else:
            T = self
        
        for v in T.preorder():
            if v.children:
                v.label = ""
                v.color = None
            v.dist = 1.0
            v.tstamp = None
            v.transferred = 0
        
        return T
    
    
    def count_node_types(self):
        
        counts = {'S': 0, 'D': 0, 'L': 0, 'H': 0, 'extant': 0}
        
        for v in self.preorder():
            
            if not v.children:
                if v.is_loss():
                    counts['L'] += 1
                else:
                    counts['extant'] += 1
            
            elif v.label == 'S':
                counts['S'] += 1
            elif v.label == 'D':
                counts['D'] += 1
            elif v.label == 'H':
                counts['H'] += 1
        
        return counts        
    

# --------------------------------------------------------------------------
#                          TREE  <--->  NEWICK
# --------------------------------------------------------------------------
        
    def to_newick(self, label=True, color=True, distance=True,
                        label_inner=True, color_inner=False):
        """Recursive PhyloTree --> Newick (str) function."""
        
        def _to_newick(node):
            
            if not node.children:
                token = ''
                if label:
                    token += str(node.label)
                if color and node.color:
                    token += '<{}-{}>'.format(*node.color) if isinstance(node.color, (tuple, list)) else '<{}>'.format(node.color)
                if distance:
                    token += ":{}".format(node.dist)
                return token
            else:
                s = ''
                for child in node.children:
                    s += _to_newick(child) + ','
                token = ''
                if label and label_inner:
                    token += str(node.label)
                if color_inner and node.color:
                    token += '<{}-{}>'.format(*node.color) if isinstance(node.color, (tuple, list)) else '<{}>'.format(node.color)
                if distance:
                    token += ':{}'.format(node.dist)
                return '({}){}'.format(s[:-1], token)
        
        
        if self.root:
            return _to_newick(self.root) + ';'
        else:
            return ';'
    
    
    @staticmethod
    def parse_newick(newick):
        """Parses trees in Newick format into object of type 'PhyloTree'."""
        
        # label<color>:distance
        label_col_dist_regex = re.compile(r"'?([a-zA-Z0-9_]*)'?<(.*)>:(-?[0-9]*\.?[0-9]*[Ee]?-?[0-9]+)")
        # label<color>
        label_col_regex = re.compile(r"'?([a-zA-Z0-9_]*)'?<(.*)>")
        # label:distance
        label_dist_regex = re.compile(r"'?([a-zA-Z0-9_]*)'?:(-?[0-9]*\.?[0-9]*[Ee]?-?[0-9]+)")
        
        id_counter = 0
        
        def parse_subtree(subroot, subtree_string):
            """Recursive function to parse the subtrees."""
            nonlocal id_counter
            children = split_children(subtree_string)
            for child in children:
                node = PhyloTreeNode(id_counter)
                subroot.add_child(node)
                id_counter += 1
                end = -1
                if child[0] == '(':                                 # the child has subtrees
                    end = child.rfind(')')
                    if end == -1:
                        raise ValueError('invalid Newick string')
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
                        raise ValueError('invalid Newick string')
                    stack -= 1
                    current += c
                else:
                    current += c
            children.append(current.strip())
            return children
        
        if not isinstance(newick, str):
            raise ValueError("Newick parser needs a 'str' as input")
        end = newick.find(";")
        if end != -1:
            newick = newick[:end]
        temp_root = PhyloTreeNode(-1)
        parse_subtree(temp_root, newick)
        if temp_root.children:
            root = temp_root.children[0]
            root.dist = 0.0                 # set distance of the root to 0
            root.detach()                   # remove the parent temp_root
                                            # (important for non-recursive to_newick2)
            return PhyloTree(root)
        else:
            raise ValueError('invalid Newick string')
    
    
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
        
        self._assert_integrity()
        G = nx.DiGraph()
        
        if not self.root:
            return G, None
        
        for v in self.preorder():
            G.add_node(v.ID, label=v.label,
                       color=v.color,
                       tstamp=v.tstamp,
                       dist=v.dist,
                       transferred=v.transferred)
        
        for u, v, sibling_nr in self.edges_sibling_order():
            if u is v or u.ID == v.ID:
                raise RuntimeError('loop at {} and {}'.format(u, v))
            G.add_edge(u.ID, v.ID)
            G.nodes[v.ID]['sibling_nr'] = sibling_nr
            
        return G, self.root.ID
    
    
    @staticmethod
    def parse_nx(G, root):
        
        number_of_leaves = 0
        
        if root is None:
            return PhyloTree(None)
    
        def build_tree(ID, parent=None):
            
            nonlocal number_of_leaves
            
            treenode = PhyloTreeNode(ID, label=G.nodes[ID]['label'])
            
            if parent:
                parent.add_child(treenode)
                
                # edge attributes are deprecated,
                # for compatibilty with older versions of AsymmeTree:
                if 'dist' in G.edges[parent.ID, ID]:
                    treenode.dist = G.edges[parent.ID, ID]['dist']
                if 'transferred' in G.edges[parent.ID, ID]:
                    treenode.transferred = G.edges[parent.ID, ID]['transferred']
            
            if 'label' in G.nodes[ID]:
                treenode.label = G.nodes[ID]['label']
            if 'color' in G.nodes[ID]:
                if isinstance(G.nodes[ID]['color'], list):
                    treenode.color = tuple(G.nodes[ID]['color'])
                else:
                    treenode.color = G.nodes[ID]['color']
            if 'tstamp' in G.nodes[ID]:
                treenode.tstamp = G.nodes[ID]['tstamp']
            if 'dist' in G.nodes[ID]:
                treenode.dist = G.nodes[ID]['dist']
            if 'transferred' in G.nodes[ID]:
                treenode.transferred = G.nodes[ID]['transferred']
            
            try:
                children_ids = sorted(G.neighbors(ID),
                                  key=lambda item: G.nodes[item]['sibling_nr'])
            except KeyError:
                # for compatibilty with older versions of AsymmeTree:
                children_ids = [item for item in G.neighbors(ID)]
            
            for child_id in children_ids:
                build_tree(child_id, parent=treenode)
            if G.out_degree(ID) == 0:
                number_of_leaves += 1
            
            return treenode
        
        tree = PhyloTree(build_tree(root))
        tree.number_of_species = number_of_leaves
        
        return tree

# --------------------------------------------------------------------------
#                           SERIALIZATION
# --------------------------------------------------------------------------
    
    @staticmethod
    def _infer_serialization_mode(filename):
        
        _, file_ext = os.path.splitext(filename)
        
        if file_ext in ('.json', '.JSON'):
            return 'json'
        elif file_ext in ('.pickle', '.PICKLE'):
            return 'pickle'
        else:
            raise ValueError('serialization format is not supplied and could '\
                             'not be inferred from file extension')
            
            
    def serialize(self, filename, mode=None):
        """Serialize the tree using pickle."""
        
        if not mode:
            mode = PhyloTree._infer_serialization_mode(filename)
        
        tree_nx, root_id = self.to_nx()
        
        if mode == 'json':
            data = nx.readwrite.json_graph.tree_data(tree_nx, root=root_id)
            
            with open(filename, 'w') as f:
                f.write( json.dumps(data) )
                
        elif mode == 'pickle':
            pickle.dump( (tree_nx, root_id), open(filename, "wb") )
            
        else:
            raise ValueError("serialization mode '{}' not supported".format(mode))
    
    
    @staticmethod
    def load(filename, mode=None):
        """Load a phylogenetic tree from a file.
        
        Using either the Python module json or pickle."""
        
        if not mode:
            mode = PhyloTree._infer_serialization_mode(filename)
        
        if mode == 'json':
            with open(filename, 'r') as f:
                data = json.loads( f.read() )
                
            tree_nx = nx.readwrite.json_graph.tree_graph(data)
            
            root_id = None
            for v in tree_nx:
                if tree_nx.in_degree(v) == 0:
                    root_id = v
                    break
                
            if root_id is None:
                raise RuntimeError('could not identify root')
                
        elif mode == 'pickle':
            tree_nx, root_id = pickle.load( open(filename, "rb") )
            
        else:
            raise ValueError("serialization mode '{}' not supported".format(mode))
        
        tree = PhyloTree.parse_nx(tree_nx, root_id)
        
        return tree
    
# --------------------------------------------------------------------------
        
    
    def copy(self):
        
        if not self.root:
            return PhyloTree(None)
        
        orig_to_new = {}
        
        for orig in self.preorder():
            new = PhyloTreeNode(orig.ID, label=orig.label, 
                                color=orig.color, dist=orig.dist,
                                tstamp=orig.tstamp,
                                transferred=orig.transferred)
            orig_to_new[orig] = new
            if orig.parent:
                orig_to_new[orig.parent].add_child(new)
        
        return PhyloTree(orig_to_new[self.root])
    
    
    @staticmethod
    def random_colored_tree(N, colors, binary=False):
        """Creates a random colored tree.
        
        The number of leaves and the color labels are specified in the
        parameters 'N' and 'colors', respectively. Each non-leaf node in the 
        resulting tree will have at least children (property of phylogenetic
        trees).
        
        Keyword arguments:
            binary - forces the tree to be binary; default is False
        """
        
        if not isinstance(N, int):
            raise TypeError("N must be of type 'int'")
            
        if isinstance(colors, int):
            colors = [i+1 for i in range(colors)]
        elif not isinstance(colors, collections.abc.Iterable):
            raise TypeError("'colors' must be of type 'int' or iterable")
        
        root = PhyloTreeNode(0, label='0')
        tree = PhyloTree(root)
        node_list = [root]
        nr, leaf_count = 1, 1
        
        while leaf_count < N:
            node = random.choice(node_list)
            
            if not node.children: 
                # to be phylogenetic at least two children must be added
                new_child1 = PhyloTreeNode(nr, label=str(nr))
                new_child2 = PhyloTreeNode(nr+1, label=str(nr+1))
                node.add_child(new_child1)
                node.add_child(new_child2)
                node_list.extend(node.children)
                nr += 2
                leaf_count += 1
            elif node.children and not binary:
                # add only one child if there are already children
                new_child = PhyloTreeNode(nr, label=str(nr))
                node.add_child(new_child)
                node_list.append(new_child)
                nr += 1
                leaf_count += 1
                
        for node in node_list:                          # assign colors
            if not node.children:                       # to leaves randomly
                node.color = random.choice(colors)
                
        return tree
   
     
# --------------------------------------------------------------------------
#                         TREE MANIPULATION
# -------------------------------------------------------------------------- 

def delete_losses_and_contract(tree, inplace=False):
    
    if not inplace:
        tree = tree.copy()
    
    loss_nodes = []
    for node in tree.postorder():
        if not node.children and node.is_loss():
            loss_nodes.append(node)
    
    # traverse from loss node to root delete if degree <= 1
    for loss_node in loss_nodes:
        current = tree.delete_and_reconnect(loss_node,
                                            add_distances=True,
                                            keep_transferred=True)
        
        while len(current.children) < 2 and current.parent:
            current = tree.delete_and_reconnect(current,
                                                add_distances=True,
                                                keep_transferred=True)
    
    return tree


def remove_planted_root(tree, inplace=True):
    
    if not inplace:
        tree = tree.copy()
        
    # delete the root if the tree is planted
    tree.remove_planted_root()
        
    if not tree.root.children and not tree.root.label:
        # no surviving genes --> return empty tree
        tree.root = None
    
    return tree