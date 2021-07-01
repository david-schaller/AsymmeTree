# -*- coding: utf-8 -*-

"""
Tree data structure for phylogentic trees.
"""

import collections, itertools, random, re, pickle, json, os
import numpy as np
import networkx as nx

from tralda.datastructures.Tree import Tree, TreeNode


__author__ = 'David Schaller'


class PhyloTreeNode(TreeNode):
    """Tree nodes for class PhyloTree.
    
    Inherits from the basic class TreeNode.
    New attributes here are:
    
    Attributes
    ----------
    color : int, tuple of two int objects, or None
        If not None, represent the ID of the node or the IDs in an edge in a
        corresponding species tree to which this node is reconciled.
    dist : float
        Evolutionary distance from the parent node.
    tstamp : float or None
        Time stamp.
    transferred : int
        Equals 1 if the edge from the parent was a transfer edge (parent must
        be an HGT event).
    
    See Also
    --------
    PhyloTree
    TreeNode (tralda package)
    CotreeNode (tralda package)
    """
    
    __slots__ = ('color', 'dist', 'tstamp', 'transferred')
    
    
    def __init__(self, ID, label='', color=None,
                 dist = 1.0, tstamp=None, transferred=0):
        """Constructor for class PhyloTreeNode.
        
        Parameters
        ----------
        ID : int
            Integer identifier of the node, some methods use the convention -1
            if an identifier is not needed.
        label: str, optional
            Human interpretable label (the default is '').
        color : int, tuple of two int objects, or None, optional
            If not None, represent the ID of the node or the IDs in an edge in
            a corresponding species tree to which this node is reconciled (the
            default is None.
        dist : float, optional
            Evolutionary distance from the parent node (the default is 1.0).
        tstamp : float or None, optional
            Time stamp (the default is None).
        transferred : int, optional
            Equals 1 if the edge from the parent was a transfer edge (the
            default is 0).
        """
        
        super().__init__(ID, label=label)
        self.color = color
        self.dist = dist
        self.tstamp = tstamp
        self.transferred = transferred

    
    def __repr__(self):
        
        return '<PhyloTN: {}>'.format(self.ID)
    
    
    def __str__(self):
        """String representation of the node.
        
        Returns
        -------
        str
            String representation of the node including color (if defined)
            and dist.
        """
        
        if isinstance(self.color, (tuple, list)):
            return '{}<{}-{}>:{}'.format(self.label, *self.color, self.dist)
        elif self.color:
            return '{}<{}>:{}'.format(self.label, self.color, self.dist)
        else:
            return '{}:{}'.format(self.label, self.dist)
        
    
    def is_loss(self):
        """Return True if the node is a loss event.
        
        AsymmeTree uses the node label '*' to mark losses.
        
        Returns
        -------
        bool
            True if the node is a loss event, else False.
        """
        
        return self.label == '*'


class PhyloTree(Tree):
    """Class for phylogentic trees.
    
    Phylogentic trees are usually defined as trees where every inner node has
    at least two children. However, the class does not force this property.
    Simulated trees are often "planted" i.e. the root has only a single child,
    and this edge represents the ancestral lineage.
    
    Attributes
    ----------
    root : PhyloTreeNode
        Inherited from class Tree, but must be a PhyloTreeNode.
    
    See Also
    --------
    PhyloTreeNode
    Tree (tralda package)
    Cotree (tralda package)
    """
    
    # corresponding node type
    node_type = PhyloTreeNode
    
    
    def __init__(self, root):
        """Constructor for class PhyloTree.
        
        Parameters
        ----------
        root : PhyloTreeNode
            The root of the PhyloTree.
        
        Raises
        ------
        TypeError
            If `root` is not an instance of PhyloTreeNode (or None).
        """
        
        if root is not None and not isinstance(root, PhyloTreeNode):
            raise TypeError("root must be of type 'PhyloTreeNode'")
        super().__init__(root)
    
    
    def sorted_nodes(self, oldest_to_youngest=True):
        """List of nodes sorted by timestamp.
        
        Return a list of all nodes in the tree sorted by time stamp beginning
        with the oldest (highest time stamp). Optionally the order can be
        reversed.
        
        Parameters
        ----------
        oldest_to_youngest : bool, optional
            If True, the nodes are sorted from oldest to youngest, otherwise
            from youngest to oldest (the default is True).
        """
        
        return sorted(self.preorder(),
                      key=lambda node: node.tstamp,
                      reverse=oldest_to_youngest)
    
    
    def sorted_edges(self):
        """List of edges (u,v) sorted by timestamp of u.
        
        Return a list of all edges (u,v) sorted by timestamp of u beginning
        with the oldest (highest time stamp)."""
    
        return sorted(self.edges(),
                      key=lambda e: (e[0].tstamp, e[1].tstamp),
                      reverse=True)
    
    
    def delete_and_reconnect(self, node,
                             add_distances=True,
                             keep_transferred=True):
        """Delete a node from the tree and reconnect its parent and children.
        
        Parameters
        ----------
        node : PhyloTreeNode
            The node in the tree to be deleted.
        add_distances : bool, optional
            When the node v is deleted and its children are reconnected to its
            parent, add to the `dist` parameter of the children `v.dist`, i.e.
            the distance of v from its parent (the default is True).
        keep_transferred : bool, optional
            When the edge of the deleted node v from its parent u was a
            transfer edge, make all edges from u to the children of v transfer
            edges (the default is True).
            
        Returns
        -------
        PhyloTreeNode or bool
            The parent of the node, if it could be deleted, or False, if the
            node could not be deleted, i.e., it has no parent.
        """
        
        if not node.parent:
            return False
        
        for child in node.children:
            if add_distances:
                child.dist += node.dist
            if keep_transferred and node.transferred:
                child.transferred = 1
        
        return super().delete_and_reconnect(node)
    
    
    def add_planted_root(self):
        """Add a new root that has the original root as its single child.
        
        Returns
        -------
        PhyloTreeNode
            The newly added root.
        """
        
        old_root = self.root
        self.root = PhyloTreeNode(-1)
        
        if old_root:
            self.root.add_child(old_root)
        
        return self.root
            
    
    def remove_planted_root(self):
        """Removes the planted root if existent.
        
        If the the tree has a planted root, its single child becomes the new
        root. Otherwise, nothing happens.
        """
        
        if len(self.root.children) == 1:
        
            new_root = self.root.children[0]
            new_root.detach()
            self.root = new_root
            new_root.dist = 0.0
            
    
    def supply_leaves(self, exclude_losses=False):
        """Leaves in the subtree rooted at each node.
        
        Computes the list of leaves for every node in the tree containing the
        leaf nodes lying in the subtree rooted at the node.
        
        Optionally, loss nodes can be excluded.
        
        Parameters
        ----------
        exclude_losses : bool, optional
            If True, all loss leaves are excluded (the default is False).
        
        Returns
        -------
        list of TreeNode objects
            The leaves under the root, i.e. the complete list of leaves.
            Returns an empty list if the root is None.
        """
        
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
    
    
    def color_sorted_leaves(self, return_list=False):
        """Sort leaves by their color.
        
        Parameters
        ----------
        return_list : bool, optional
            If True, also return a list of leaves such that leaves of the same
            color appear as a coherent sequence (the default is False).
        
        Returns
        -------
        dict
            A dictionary with colors as keys and a list of corresponding nodes
            as values.
        list, optional
            List of leaves such that leaves of the same color appear as a
            coherent sequence.
        """
        
        color_dict = {}
        
        for leaf in self.leaves():
            if leaf.color not in color_dict:
                color_dict[leaf.color] = []
            color_dict[leaf.color].append(leaf)
        
        if not return_list:
            return color_dict
        else:
            leaves = []
            for color, leaf_list in color_dict.items():
                for leaf in leaf_list:
                    leaves.append(leaf)
            
            return color_dict, leaves
    
    
    def distance_matrix(self, leaf_order=None):
        """Distance matrix on the leaf set of the phylogenetic tree.
        
        Computes a distance matrix on the set of leaves of the tree where each
        distances is the sum of the distances (`dist`) on the path connecting
        the pair of leaves.
        
        Additionally a list of leaves corresponding to the indices in the
        matrix is returned.
        
        Parameters
        ----------
        leaf_order : list, optional
            A list of all leaves in the tree defining the indices for the
            matrix (the default is None, in which case leaves are indexed in
            sibling order).
        
        Returns
        -------
        list of PhyloTreeNode objects
            Represents the order for the lines/columns in the distance matrix.
        numpy.ndarray (dtype=numpy.float)
            The distance matrix.
        """
        
        distance_dict = self.distances_from_root()
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
        """The distances of each node to the root of the tree.
        
        Returns
        -------
        dict
            The keys are PhyloTreeNode objects and the values their distances
            (sum of `dist`) to the root.
        """
        
        distance_dict = {}
        
        for v in self.preorder():
            if not v.parent:
                distance_dict[v] = 0.0
            else:
                depth = distance_dict[v.parent] + v.dist
                distance_dict[v] = depth
                
        return distance_dict
    
    
    def topology_only(self, inplace=True):
        """Reset colors, distances, time stamps, transfer status and inner
        labels.
        
        Parameters
        ----------
        inplace : bool
            If True, reset the attributes of this tree instance, otherwise
            make a copy first and modify the copy.
        
        Returns
        -------
        PhyloTree
            The original or a copy of the tree instance with reset attributes.
        """
        
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
        """Count speciations, duplication, losses, HGTs and surviving genes.
        
        Returns
        -------
        dict
            With the event counts as values. Key are 'S' (speciations), 'D'
            (duplication), 'L' (losses), 'H' (HGTs) and 'extant' (surviving
            genes).
        """
        
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
        """Return a Newick representation of the tree.
        
        This function overrides the function of the parent class.
        
        Parameters
        ----------
        label : bool, optional
            If True, the Newick str contains the labels of the nodes (the
            default is True).
        color : bool, optional
            If True, the Newick str contains the colors of the nodes in
            <[...]> brackets (the default is True).
        distance : bool, optional
            If True, the Newick str contains the distances of the nodes in
            standard :[...] notation (the default is True).
        label_inner : bool, optional
            If True, the Newick str also contains the labels of the inner 
            nodes (the default is True).
        color_inner : bool, optional
            If True, the Newick str contains the colors of the inner nodes (the
            default is False).
        
        Returns
        -------
        str
            A Newick representation of the tree.
        """
        
        def _to_newick(node):
            
            if not node.children:
                token = ''
                if label:
                    token += str(node.label)
                if color and node.color:
                    token += '<{}-{}>'.format(*node.color) \
                        if isinstance(node.color, (tuple, list)) \
                        else '<{}>'.format(node.color)
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
                    token += '<{}-{}>'.format(*node.color) \
                        if isinstance(node.color, (tuple, list)) \
                        else '<{}>'.format(node.color)
                if distance:
                    token += ':{}'.format(node.dist)
                return '({}){}'.format(s[:-1], token)
        
        
        if self.root:
            return _to_newick(self.root) + ';'
        else:
            return ';'
    
    
    @staticmethod
    def parse_newick(newick):
        """Parses trees in Newick format into object of type 'PhyloTree'.
        
        Parameters
        ----------
        newick : str
            A tree in Newick format.
        
        Returns
        -------
        PhyloTree
            The parsed tree.
        
        Raises
        ------
        TypeError
            If the input is not a string.
        ValueError
            If the input is not a valid Newick string.
        
        Notes
        -----
        Do not use this function for serialization and reloading PhyloTree
        objects. Use the `serialize()` function instead.
        IDs are set arbitrarily and do not necessarily match the labels.
        The colors (if present in <...> in the string) are parsed as strings
        and need to be converted to integers afterwards if necessary.
        """
        
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
            raise TypeError("Newick parser needs a 'str' as input")
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
    
# --------------------------------------------------------------------------
#                         TREE  <--->  NETWORKX
# --------------------------------------------------------------------------
            
    def to_nx(self):
        """Convert a PhyloTree into a NetworkX Graph.
        
        The attributes correspond to the node attributes in the resulting graph.
        
        Returns
        -------
        networkx.Graph
            A graph representation of the tree.
        int
            The ID of the root in order to be able to completely reconstruct
            the tree.
        """
        
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
        """Convert a NetworkX Graph version back into a PhyloTree.
        
        Parameters
        ----------
        G : networkx.Graph
            A tree represented as a Networkx Graph.
        root : int
            The ID of the root.
        
        Returns
        -------
        PhyloTree
            The reconstructed tree.
        """
        
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
        
        if file_ext.lower() == '.json':
            return 'json'
        elif file_ext.lower() == '.pickle':
            return 'pickle'
        else:
            raise ValueError('serialization format is not supplied and could '\
                             'not be inferred from file extension')
            
            
    def serialize(self, filename, mode=None):
        """Serialize the tree using pickle or json.
        
        Parameters
        ----------
        filename : str
            The filename (including the path) of the file to be created.
        mode : str or None, optional
            The serialization mode. Supported are pickle and json. The default
            is None in which case the mode is inferred from the file extension.
        
        Raises
        ------
        ValueError
            If the serialization mode is unknown or could not be inferred.
        """
        
        if not mode:
            mode = PhyloTree._infer_serialization_mode(filename)
        
        tree_nx, root_id = self.to_nx()
        
        if mode == 'json':
            data = nx.readwrite.json_graph.tree_data(tree_nx, root=root_id)
            
            with open(filename, 'w') as f:
                f.write( json.dumps(data) )
                
        elif mode == 'pickle':
            pickle.dump( (tree_nx, root_id), open(filename, 'wb') )
            
        else:
            raise ValueError("serialization mode '{}' not supported".format(mode))
    
    
    @staticmethod
    def load(filename, mode=None):
        """Reload a PhyloTree from a file (pickle or json).
        
        Parameters
        ----------
        filename : str
            The filename (including the path) of the file to be loaded.
        mode : str or None, optional
            The serialization mode. Supported are pickle and json. The default
            is None in which case the mode is inferred from the file extension.
        
        Returns
        -------
        PhyloTree
            The tree reloaded from file.
        
        Raises
        ------
        ValueError
            If the serialization mode is unknown or could not be inferred.
        """
        
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
            tree_nx, root_id = pickle.load( open(filename, 'rb') )
            
        else:
            raise ValueError("serialization mode '{}' not supported".format(mode))
        
        tree = PhyloTree.parse_nx(tree_nx, root_id)
        
        return tree
    
# --------------------------------------------------------------------------
#                    RECONSTRUCTION OF INFORMATION
# --------------------------------------------------------------------------

    def reconstruct_IDs(self):
        """Reconstruct the IDs from the (str) labels.
        
        If the labels of the nodes can be converted into integers, then the
        nodes receive as ID their converted labels.
        After the first traversal, all remaining nodes that still have an
        undefined ID (-1) obtain are unique integer ID that has not yet been
        assigned in the first traversal.
        """
        
        self.number_of_species = 0
        IDs = set()
        
        for v in self.preorder():
            if not v.children:
                self.number_of_species += 1
            
            # convention for undefined ID is -1
            if not isinstance(v.ID, int) or v.ID < 0:
                if v.label.isdigit():
                    v.ID = int(v.label)
                    IDs.add(v.ID)
            elif isinstance(v.ID, int) and v.ID >= 0:
                IDs.add(v.ID)
        
        # assign new IDs to remaining nodes
        current_ID = 0
        for v in self.preorder():
            if not isinstance(v.ID, int) or v.ID < 0:
                while current_ID in IDs:
                    current_ID += 1
                v.ID = current_ID
                IDs.add(current_ID)
    
    
    def reconstruct_info_from_graph(self, G):
        """Reconstruct the labels and colors from a NetworkX Graph.
        
        Parameters
        ----------
        G : networkx.Graph
            The graph from which labels and colors shall be reconstructed.
        """
        
        for v in self.preorder():
            if v.ID in G:
                if 'label' in G.nodes[v.ID]:
                    v.label = G.nodes[v.ID]['label']
                if 'color' in G.nodes[v.ID]:
                    if isinstance(G.nodes[v.ID]['color'], list):
                        v.color = tuple(G.nodes[v.ID]['color'])
                    else:
                        v.color = G.nodes[v.ID]['color']
                
                
    def reconstruct_timestamps(self):
        """Reconstruct the timestamps.
        
        Make the time stamps matching with the distance attribute. The root
        obtains time stamp 1.0, and all other node smaller time stamps such
        that the difference to the parent's time stamp is exactly `dist`.
        """
        
        self.root.tstamp = 1.0
        for v in self.preorder():
            if v.parent:
                v.tstamp = v.parent.tstamp - v.dist
    
# --------------------------------------------------------------------------
        
    
    def copy(self, mapping=False):
        """Return a copy of the tree.
        
        Constructs a deep copy of the tree, i.e. to the level of nodes.
        By default, the node attributes are all immutable data types (except
        `leaves` which is not copied, and may be recomputed later for the copy).
        Hence, the original tree is not affected by operations on the copy.
        
        Parameters
        ----------
        mapping : bool
            If True, additionally return the mapping from original to copied
            nodes as dictionary.
        
        Returns
        -------
        PhyloTree or tuple of PhyloTree and dict
            A copy of the tree and optionally the mapping from original to 
            copied nodes.
        """
        
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
        
        if mapping:
            return PhyloTree(orig_to_new[self.root]), orig_to_new
        else:
            return PhyloTree(orig_to_new[self.root])
    
    
    @staticmethod
    def to_phylotree(tree, mapping=False):
        """Convert a Tree instance into a PhyloTree.
        
        By default, the node attributes are all immutable data types (except
        `leaves` which is not copy, and may be recomputed later for the copy).
        Hence, the original tree is not affected by operations on the copy.
        
        Parameters
        ----------
        tree : Tree
        mapping : bool
            If True, additionally return the mapping from original to the
            newly created nodes (of type 'PhyloTreeNode') as dictionary.
        
        Returns
        -------
        PhyloTree or tuple of PhyloTree and dict
            A PhyloTree version of the tree and optionally the mapping from
            original to copied nodes.
        """
        
        if not isinstance(tree, Tree) or isinstance(tree, PhyloTree):
            raise TypeError("requires an instance of type 'Tree' (but not "\
                            "'PhyloTree', use copy() for this purpose)")
        
        if not tree.root:
            return PhyloTree(None)
        
        orig_to_new = {}
        
        for orig in tree.preorder():
            new = PhyloTreeNode(orig.ID, label=orig.label)
            orig_to_new[orig] = new
            if orig.parent:
                orig_to_new[orig.parent].add_child(new)
        
        if mapping:
            return PhyloTree(orig_to_new[tree.root]), orig_to_new
        else:
            return PhyloTree(orig_to_new[tree.root])
    
    
    @staticmethod
    def random_colored_tree(N, colors, binary=False, force_all_colors=False):
        """Create a random colored tree.
        
        The number of leaves and the color labels are specified in the
        parameters `N` and `colors`, respectively. Each non-leaf node in the 
        resulting tree will have at least children (property of phylogenetic
        trees).
        
        Parameters
        ----------
        N : int
            The desired number of leaves.
        colors : int or list
            The list of colors, or the desired number of colors in which case
            the colors {1, ..., colors} are used.
        binary : bool, optional
            If True, forces the tree to be binary (the default is False).
        force_all_colors : bool
            If True, the resulting tree is guaranteed to have at least one leaf
            of each color (the default is False).
        
        Returns
        -------
        PhyloTree
            A random tree with `N` leaves and random leaf coloring.
        
        Raises
        ------
        TypeError
            If `N` is not an integer > 0.
        ValueError
            If the number of colors is greater than `N` and `force_all_colors`
            is true.
        """
        
        if not isinstance(N, int) or N < 1:
            raise TypeError("N must be of type 'int' and > 0")
            
        if isinstance(colors, int):
            colors = [i+1 for i in range(colors)]
        elif not isinstance(colors, collections.abc.Iterable):
            raise TypeError("'colors' must be of type 'int' or iterable")
            
        if len(colors) > N and force_all_colors:
            raise ValueError('cannot force all colors since #colors > N')
        
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
        
        leaves = [node for node in node_list if not node.children]
        
        if force_all_colors:
            # use every color at least once
            permutation = np.random.permutation(len(leaves))
            for i in range(len(leaves)):
                if i < len(colors):
                    leaves[permutation[i]].color = colors[i]
                else:
                    # color the remaining leaves randomly
                    leaves[permutation[i]].color = random.choice(colors)
        else:
            # assign colors completely randomly
            for leaf in leaves:
                leaf.color = random.choice(colors)
                
        return tree
   
     
# --------------------------------------------------------------------------
#                         TREE MANIPULATION
# -------------------------------------------------------------------------- 

def delete_losses_and_contract(tree, inplace=False):
    """Delete all branches leading to loss leaves only.
    
    Nodes that would have only a single child afterwards are suppressed (except
    possibly the planted root), i.e. their children are recursively reconnected
    to the parents. Distances are cumulated in this process and the transferred
    status is kept in the sense that an edge is a transfer edge if at least
    one edge on the contracted path to which it corresponds was a transfer edge.
    
    Parameters
    ----------
    tree : PhyloTree
        The tree in which loss branches shall be removed.
    inplace : bool, optional
        If True, the tree is directly manipulated. The default is False in
        which case a copy of the tree is created which gets manipulated while
        the original tree remains untouched.
    
    Returns
    -------
    PhyloTree
        The tree with all loss branches removed (original instance or a new one
        depending on the `inplace` parameter).
    """
    
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
    """Remove the planted root of the tree (if existent).
    
    Parameters
    ----------
    tree : PhyloTree
        The tree in which the planted root shall be removed.
    inplace : bool, optional
        If True, the tree is directly manipulated, otherwise a copy is created
        and the original tree remains untouched (the default is True).
    
    Returns
    -------
    PhyloTree
        The tree with the planted root removed (original instance or a new one
        depending on the `inplace` parameter).
    """
    
    if not inplace:
        tree = tree.copy()
        
    # delete the root if the tree is planted
    tree.remove_planted_root()
        
    if not tree.root.children and not tree.root.label:
        # no surviving genes --> return empty tree
        tree.root = None
    
    return tree