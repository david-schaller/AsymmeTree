# -*- coding: utf-8 -*-

"""HGT, xenology and Fitch graphs.

This module contains several functions for the analysis of horizontal gene
transfer events in the simulated scenarios. In particular, the directed and
undirected Fitch graph can be extracted, as well as the pairs of genes that
diverged later than the respective species in which they reside, i.e. the
later-divergence-time (LDT) graph, see [1]. The latter situation is indicative 
or the presence of HGT events in the scenario.

References
----------
.. [1] D. Schaller, M. Lafond, P.F. Stadler, N. Wieseke, M. Hellmuth.
   Indirect identification of horizontal gene transfer.
   In: Journal of Mathematical Biology, 2021, 83(1):10.
   doi: 10.1007/s00285-021-01631-0.
"""

import itertools

import networkx as nx

from tralda.datastructures import Tree, TreeNode, LCA
from tralda.tools.GraphTools import independent_sets
from tralda.cograph.Cograph import (to_cotree, paths_of_length_2)
from tralda.supertree import Build

from asymmetree.tools.PhyloTreeTools import (add_planted_root,
                                             phylo_tree_attributes,
                                             assign_missing_labels)
from asymmetree.treeevolve.SpeciesTree import distance_from_timing


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                     Transfer edges and Fitch graphs
# --------------------------------------------------------------------------

def true_transfer_edges(T):
    """Returns a set containing v if (u, v) is labeled as a transfer edge.
    
    Parameters
    ----------
    T : Tree
        A tree whose nodes have the 'transferred' attribute.
    
    Returns
    -------
    set
        The subset of TreeNode instances in the tree for which transferred is 1.
    """
    
    return {v for _, v in T.edges() if v.transferred}


def rs_transfer_edges(T, S, lca_S=None):
    """Transfer edges in T according to the relaxed scenario definition.
    
    An edge (u,v) in T is an (rs-)transfer edge if u and v are mapped to
    incomparable nodes/edges in the species tree S.
    
    Parameters
    ----------
    T : Tree
        The gene tree.
    S : Tree
        The species tree corresponding to the gene tree.
    lca_S : tralda.datastructures.Tree.LCA, optional
        Instance of LCA corresponding to the species tree, default is None
        in which case a new instance is created and used.
    
    Returns
    -------
    set
        The subset of TreeNode instances in the tree for which the 'reconc'
        attributes of its incident edges refer to incomparable branches in the
        species tree.
    """
    
    if not isinstance(lca_S, LCA):
        lca_S = LCA(S)
        
    transfer_edges = set()
    
    for u, v in T.edges():
        if not lca_S.are_comparable(u.reconc, v.reconc):
            transfer_edges.add(v)
    
    return transfer_edges


def fitch(tree, transfer_edges, supply_undirected=False, lca_T=None):
    """Returns the (directed) Fitch graph.
    
    Parameters
    ----------
    tree : Tree
        The gene tree.
    transfer_edges : iterable
        The subset of TreeNode instances in the tree that are transfer edges.
    supply_undirected : bool, optional
        If True, additionally return the undirected Fitch graph, default is
        False.
    lca_T : tralda.datastructures.Tree.LCA, optional
        Instance of LCA corresponding to the tree, default is False in which
        case a new instance is created and used.
        
    Returns
    -------
    networkx.DiGraph or tuple
        The directed Fitch graph or a tuple containing the directed Fitch graph
        (networkx.DiGraph) and the undirected Fitch graph (networkx.Graph).
    """
    
    if not isinstance(lca_T, LCA):
        lca_T = LCA(tree)
        
    if not isinstance(transfer_edges, (set, dict)):
        transfer_edges = set(transfer_edges)
    
    leaves = [l for l in tree.leaves()]
    fitch = nx.DiGraph()
    
    # store for each leaf the first transfer edge on the way to the root
    first_transfer = {}
    
    for x in leaves:
        fitch.add_node(x.label, color=x.reconc)
        
        current = x
        while current:
            if current in transfer_edges:
                first_transfer[x] = current
                break
            current = current.parent
    
    for x, y in itertools.permutations(leaves, 2):
        
        if (y in first_transfer and
            lca_T.ancestor_not_equal(lca_T(x, y), first_transfer[y])):
            fitch.add_edge(x.label, y.label)
    
    if not supply_undirected:
        return fitch
    else:
        return fitch, fitch.to_undirected()
    

def undirected_fitch(tree, transfer_edges, lca_T=None):
    """Returns the undirected Fitch graph.
    
    Parameters
    ----------
    tree : Tree
        The gene tree.
    transfer_edges : iterable
        The subset of TreeNode instances in the tree that are transfer edges.
    lca_T : tralda.datastructures.Tree.LCA, optional
        Instance of LCA corresponding to the tree, default is False in which
        case a new instance is created and used.
        
    Returns
    -------
    networkx.Graph
        The undirected Fitch graph (networkx.Graph).
    """
    
    return fitch(tree, transfer_edges, supply_undirected=True, lca_T=lca_T)[1]


def _rs_fitch_aux_graph(G, color_set, ind_sets, ignore_set):
    
    aux_graph = nx.Graph()
    aux_graph.add_nodes_from(color_set)
    
    for i, ind_set in enumerate(ind_sets):
        
        # skip the independent set that shall not be included
        if ignore_set == i:
            continue
        
        for x, y in itertools.combinations(ind_set, 2):
            X = G.nodes[x]['color']
            Y = G.nodes[y]['color']
            
            if X != Y:
                aux_graph.add_edge(X, Y)
    
    return aux_graph


def is_rs_fitch(G, color_set=None):
    """Checks whether a given graph is an rs-Fitch graph.
    
    I.e. whether the graph is the Fitch graph of some relaxed scenario (rs),
    see [1].
    
    Parameters
    ----------
    G : networkx.Graph
        A graph whose nodes have the 'color' attribute.
    color_set : set, optional
        The set of colors appearing in the graph. The default is None, in which
        case such a set is created internally.
    
    Returns
    -------
    bool
        True if the graph is an rs-Fitch graph as defined in [1].
        
    References
    ----------
    .. [1] D. Schaller, M. Lafond, P.F. Stadler, N. Wieseke, M. Hellmuth.
       Indirect identification of horizontal gene transfer.
       In: Journal of Mathematical Biology, 2021, 83(1):10.
       doi: 10.1007/s00285-021-01631-0.
    """
    
    if color_set is None:
        color_set = set()
        for v in G.nodes():
            color_set.add(G.nodes[v]['color'])
    
    ind_sets = independent_sets(G)
    
    # not a complete multipartite graph
    if ind_sets is False:
        return False
    
    # trivial number of independent sets
    k = len(ind_sets)
    if k <= 1:
        return True
    
    # more than one independent set
    for ignore_set in range(-1, k):
        
        # ignore_set = -1 means that no set is ignored
        aux_graph = _rs_fitch_aux_graph(G, color_set, ind_sets, ignore_set)
        if not nx.is_connected(aux_graph):
            return True
    
    # if nothing was True
    return False


# --------------------------------------------------------------------------
#                   Later-divergence-time (LDT) graphs
# --------------------------------------------------------------------------

def below_equal_above(T, S, lca_T=None, lca_S=None):
    """Detemine the pairs of genes that diverged later, at the same time, or 
    earlier as their species.
    
    Returns three graphs with the leaves of the gene tree T as vertex set, and 
    edges ab if and only if a and b diverged later, at the same time, or later,
    resp., than the corresponding species A and B in the species tree S.
    
    Parameters
    ----------
    T : Tree
        The gene tree.
    S : Tree
        The corresponding species tree.
    lca_T : tralda.datastructures.Tree.LCA, optional
        Instance of LCA corresponding to the gene tree, default is False in
        which case a new instance is created and used.
    lca_S : tralda.datastructures.Tree.LCA, optional
        Instance of LCA corresponding to the species tree, default is False in
        which case a new instance is created and used.
    
    Returns
    -------
    tuple of three networkx.Graph instances
        The below, equal, and above relation as described above.
    """
    
    L_T = [l for l in T.leaves()]
    L_S = {l.label: l for l in S.leaves()}
    
    below = nx.Graph()
    equal = nx.Graph()
    above = nx.Graph()
    for u in L_T:
        below.add_node(u.label, color=u.reconc)
        equal.add_node(u.label, color=u.reconc)
        above.add_node(u.label, color=u.reconc)
    
    if not isinstance(lca_T, LCA):
        lca_T = LCA(T)
    if not isinstance(lca_S, LCA):
        lca_S = LCA(S)
    
    for a, b in itertools.combinations(L_T, 2):
        
        t_ab = lca_T.get(a, b).tstamp
        t_AB = lca_S.get(L_S[a.reconc], L_S[b.reconc]).tstamp
        
        if t_ab < t_AB:
            below.add_edge(a.label, b.label)
        elif t_ab == t_AB:
            equal.add_edge(a.label, b.label)
        else:
            above.add_edge(a.label, b.label)
    
    
    return below, above, equal


def ldt_graph(T, S, lca_T=None, lca_S=None):
    """Later-divergence-time graph, see [1].
    
    Returns a graph with the leaves of the gene tree T as vertex set, and 
    edges ab if and only if a and b diverged later than the corresponding
    species A and B in the species tree S.
    
    Parameters
    ----------
    T : Tree
        The gene tree.
    S : Tree
        The corresponding species tree.
    lca_T : tralda.datastructures.Tree.LCA, optional
        Instance of LCA corresponding to the gene tree, default is False in
        which case a new instance is created and used.
    lca_S : tralda.datastructures.Tree.LCA, optional
        Instance of LCA corresponding to the species tree, default is False in
        which case a new instance is created and used.
    
    References
    ----------
    .. [1] D. Schaller, M. Lafond, P.F. Stadler, N. Wieseke, M. Hellmuth.
       Indirect identification of horizontal gene transfer.
       In: Journal of Mathematical Biology, 2021, 83(1):10.
       doi: 10.1007/s00285-021-01631-0.
    """
    
    return below_equal_above(T, S, lca_T=lca_T, lca_S=lca_S)[0]


def _simulate_timing(tree):
    """Simulates a timing for the tree.
    
    It is t(root) = 1 and t(x) = 0 for x in L(S)."""
    
    max_depth = {}
    
    for v in tree.postorder():
        if not v.children:
            max_depth[v] = 0
        else:
            max_depth[v] = max(max_depth[c] for c in v.children) + 1
    
    for v in tree.preorder():
        if not v.children:
            v.tstamp = 0.0
        elif not v.parent:
            v.tstamp = 1.0
        else:
            v.tstamp = v.parent.tstamp * (1 - 1 / (max_depth[v]+1))


class RsScenarioConstructor:
    """Construct an rs-scenario for a given LDT graph.
    
    Implementation of the algorithm presented in [1].
    
    References
    ----------
    .. [1] D. Schaller, M. Lafond, P.F. Stadler, N. Wieseke, M. Hellmuth.
       Indirect identification of horizontal gene transfer.
       In: Journal of Mathematical Biology, 2021, 83(1):10.
       doi: 10.1007/s00285-021-01631-0.
    """
    
    def __init__(self, colored_cograph, color_set=None):
        """
        Parameters
        ----------
        colored_cograph : networkx.Graph
            A cograph whose nodes have the 'color' attribute.
        color_set : set, optional
            The set of colors appearing in the graph. The default is None,
            in which case such a set is created internally.
        """
        
        self.G = colored_cograph
        
        self.L = [x for x in self.G.nodes()]
        
        if color_set is None:
            color_set = set()
            for v in self.L:
                color_set.add(self.G.nodes[v]['color'])
        elif not isinstance(color_set, (set, dict)):
            color_set = set(color_set)
        self.color_set = color_set
        
        
    def run(self):
        """Construct an rs-scenario consisting of a species and gene tree.
        
        Returns
        -------
        tuple of two Tree instances or bool
            The species and gene tree, or False if the supplied graph was 
            not an LDT graph.
        """
        
        self.S = self._species_tree()
        if not self.S:
            return False
        
        self.epsilon = float('inf')
        for u, v in self.S.edges():
            self.epsilon = min(self.epsilon,
                               abs(u.tstamp - v.tstamp))
        self.epsilon /= 3.0
        
        self.S_leaves = self.S.leaf_dict()
        
        # top-level call on full leaf set and rho_S
        self.T = Tree(self._build_gene_tree(self.L, self.S.root.children[0]))
        
        # --- post-processing of T ---
        
        # planted root
        add_planted_root(self.T)
        self.T.root.tstamp = self.S.root.tstamp
        self.T.root.reconc = self.S.root.label
        
        # suppress all vertices with a single child except the planted root
        to_suppress = []
        for u in self.T.preorder():
            if u.parent and len(u.children) == 1:
                to_suppress.append((u.parent, u))
        self.T.contract(to_suppress)
        
        distance_from_timing(self.T)
        
        return self.S, self.T
        
    
    def _species_tree(self):
        
        cotree = to_cotree(self.G)
        
        if not cotree:
            return False
        
        triples = set()
        
        for a, b, c in paths_of_length_2(cotree):
            
            A = self.G.nodes[a.label]['color']
            B = self.G.nodes[b.label]['color']
            C = self.G.nodes[c.label]['color']
            
            if A != B and A != C and B != C:
                
                # sort to avoid redundancy in triple set
                if A > C:
                    A, C = C, A
                
                triples.add( (A, C, B) )
        
        build = Build(triples, self.color_set, mincut=False)
        S = build.build_tree()
        
        if not S:
            return False
        else:
            phylo_tree_attributes(S, inplace=True)
            add_planted_root(S)
            assign_missing_labels(S)
            _simulate_timing(S)
            distance_from_timing(S)
            return S
    
    
    def _build_gene_tree(self, L, u_S):
        
        u_T = TreeNode(label='', event='D',
                       reconc=(u_S.parent.label, u_S.label),
                       tstamp=u_S.tstamp + self.epsilon)
        
        # u_S is a leaf
        if not u_S.children:
            for x in L:
                leaf = TreeNode(label=x,
                                reconc=self.G.nodes[x]['color'],
                                tstamp=0.0)
                u_T.add_child(leaf)
        
        # u_S is an inner vertex
        else:
            
            # maps color to the respective child of u_S
            color_to_v_S = {}
            for v_S in u_S.children:
                for x in self.S_leaves[v_S]:
                    color_to_v_S[x.label] = v_S
                    
            # connected components C of G[L']
            for C in nx.connected_components(self.G.subgraph(L)):
                
                # choose v_S_star s.t. L(S(v_S_star)) \cap sigma(C) is non-empty
                v_S_star = color_to_v_S[self.G.nodes[next(iter(C))]['color']]
                v_T = TreeNode(label='', event='H',
                               reconc=(u_S.label, v_S_star.label),
                               tstamp=u_S.tstamp - self.epsilon)
                u_T.add_child(v_T)
                
                # aux. graph for equivalence classes
                aux_graph = nx.Graph()
                aux_graph.add_nodes_from(C)
                
                for x, y in itertools.combinations(C, 2):
                    if (color_to_v_S[self.G.nodes[x]['color']] is
                        color_to_v_S[self.G.nodes[y]['color']]):
                        aux_graph.add_edge(x, y)
                
                # for each equivalence class K that is a subset of C
                for K in nx.connected_components(aux_graph):
                    
                    # choose v_S s.t. \sigma(K) \subseteq L(S(v_S))
                    v_S = color_to_v_S[self.G.nodes[next(iter(K))]['color']]
                    
                    w_K = self._build_gene_tree(K, v_S)
                    v_T.add_child(w_K)
                    
                    if v_S is not v_S_star:
                        w_K.transferred = 1
        
        return u_T
    
    
# --------------------------------------------------------------------------
#                  Compatibility of trees and partitions
# --------------------------------------------------------------------------

def is_compatible(T, partition, lca=None):
    """Checks whether a tree and a partition are compatible.
    
    A tree with leaf set L and a partition of L are compatible if there is a
    subset of tree edges such that the forest obtained by their removal induces
    the given partition (in terms of membership to the same connected
    component).
    In this case, the vertices can be uniquely colored by set A of the
    partition if they lie on a path connecting two elements from A.
    
    Parameters
    ----------
    T : Tree
        A tree with unique leaf label.
    partition : list or tuple
        A partition of the labels of the trees' leaves.
    lca : LCA, optional
        Precomputed LCA datastructure for the tree. The default is None in
        which case it is computed internally.
    
    Returns
    -------
    bool
        Returns False if the tree and the partition are incompatible, and
        otherwise a tuple of a dict and a list.
    dict
        A dictionary containing vertices of the tree that lie on a path
        connecting two elements from the same set of the partition with the
        index of this set as value.
    list
        A list containing the last common ancestor in the tree for each set
        in the partition.
    """
    
    if not isinstance(T, Tree):
        raise TypeError("T must be of type 'Tree'")
    
    if not isinstance(partition, (list, tuple)):
        raise TypeError("partition must be of type 'list' or 'tuple'")
    
    if not lca:
        lca = LCA(T)
    
    # map labels to leaves of the tree
    label_to_leaf = {v.label: v for v in T.leaves()}
    
    # color each vertex if lies on the path connecting two vertices from the
    # same set
    vertex_coloring = {}
    
    # last common ancestors of the sets in the partition
    lcas = []
    
    for i, A in enumerate(partition):
        
        A_iterator = iter(A)
        x_1 = label_to_leaf[next(A_iterator)]
        current_lca = x_1
        visited = {x_1}
        vertex_coloring[x_1] = i
        
        for x_label in A_iterator:
            x = label_to_leaf[x_label]
            new_lca = lca(x, current_lca)
            
            starts = [x]
            if new_lca is not current_lca:
                starts.append(current_lca.parent)
                current_lca = new_lca
            
            for current in starts:
            
                while current:
                    if current in visited:
                        break
                    
                    if current in vertex_coloring and \
                        vertex_coloring[current] != i:
                        return False
                    else:
                        vertex_coloring[current] = i
                        visited.add(current)
                    
                    if current is new_lca:
                        break
                    
                    current = current.parent
        
        lcas.append(current_lca)
    
    return vertex_coloring, lcas


def is_refinement_compatible(T, partition, lca=None):
    """Checks whether a tree and a partition are refinement-compatible.
    
    A tree T with leaf set L and a partition of L are refinement-compatible if
    there is a refinement T' of T for which a subset of edges of T' exists
    such that the forest obtained by their removal induces the given partition
    (in terms of membership to the same connected component).
    In this case, the edges can be uniquely colored by set A of the partition
    if they lie on a path connecting two elements from A.
    
    Parameters
    ----------
    T : Tree
        A tree with unique leaf label.
    partition : list or tuple
        A partition of the labels of the trees' leaves.
    lca : LCA, optional
        Precomputed LCA datastructure for the tree. The default is None in
        which case it is computed internally.
    
    Returns
    -------
    bool
        Returns False if the tree and the partition are not refinement-
        compatible, and otherwise a tuple of a dict and a list.
    dict
        A dictionary containing nodes v of the tree such that the edge
        (parent(v), v) lies on a path connecting two elements from the same
        set of the partition with the index of this set as value.
    list
        A list containing the last common ancestor in the tree for each set
        in the partition.
    """
    
    if not isinstance(T, Tree):
        raise TypeError("T must be of type 'Tree'")
    
    if not isinstance(partition, (list, tuple)):
        raise TypeError("partition must be of type 'list' or 'tuple'")
    
    if not lca:
        lca = LCA(T)
    
    # map labels to leaves of the tree
    label_to_leaf = {v.label: v for v in T.leaves()}
    
    # color each vertex v if the edge (v.parent, v) lies on the path
    # connecting two vertices from the same set
    edge_coloring = {}
    
    # last common ancestors of the sets in the partition
    lcas = []
    
    for i, A in enumerate(partition):
        
        A_iterator = iter(A)
        current_lca = label_to_leaf[next(A_iterator)]
        visited = set()
        
        for x_label in A_iterator:
            x = label_to_leaf[x_label]
            new_lca = lca(x, current_lca)
            
            starts = [x]
            if new_lca is not current_lca:
                starts.append(current_lca)
                current_lca = new_lca
            
            for current in starts:
            
                while current:
                    
                    if current in edge_coloring and \
                        edge_coloring[current] != i:
                        return False
                    else:
                        edge_coloring[current] = i
                        visited.add(current)
                        
                    # parent is always defined since at least the first edge
                    # must be added and if the root is reached, it must equal
                    # new_lca
                    current = current.parent
                    if current in visited or current is new_lca:
                        break
        
        lcas.append(current_lca)
    
    return edge_coloring, lcas


def fitch_orientation(T, partition, lca=None):
    """(Un)ambiguously present or absent edges in the directed Fitch graph.
    
    The ordered pair (x, y) of two leaves x and y is in the Fitch relation
    if the path from lca(x, y) to y contains a transfer edge. An undirected
    Fitch relation ((x, y) if (x, y) or (y, x) in the directed relation) can
    be represented by a partition (the maximal independent sets).
    In combination with a tree, this partition determines to some extent the
    directed Fitch relation.
    
    Parameters
    ----------
    T : Tree
        A tree with unique leaf label.
    partition : list or tuple
        A partition of the labels of the trees' leaves.
    lca : LCA, optional
        Precomputed LCA datastructure for the tree. The default is None in
        which case it is computed internally.
    
    Returns
    -------
    bool
        Returns False if the tree and the partition are incompatible.
    2-dimensional array (list of lists)
        Stores at indices i and j corresponding to sets A and B, resp., of the
        partition whether (A, B) is 'essential', 'forbidden', or 'ambiguous'.
    """
    
    if not isinstance(T, Tree):
        raise TypeError("T must be of type 'Tree'")
    
    if not isinstance(partition, (list, tuple)):
        raise TypeError("partition must be of type 'list' or 'tuple'")
    
    if not lca:
        lca = LCA(T)
    
    compatibility = is_compatible(T, partition, lca=lca)
    if compatibility:
        vertex_coloring, lcas = compatibility
    else:
        return False
    
    # compute for each vertex its lowest colored (strict) ancestor
    lowest_colored_ancestor = {}
    for v in T.preorder():
        if not v.parent:
            lowest_colored_ancestor[v] = None
        elif v.parent in vertex_coloring:
            lowest_colored_ancestor[v] = v.parent
        else:
            lowest_colored_ancestor[v] = lowest_colored_ancestor[v.parent]
    
    n = len(partition)
    matrix = [[None for j in range(n)] for i in range(n)]
    
    for i in range(n):
        for j in range(n):
            
            if i == j:
                matrix[i][j] = 'forbidden'
                continue
            
            u = lca(lcas[i], lcas[j])       # lca(lca(A), lca(B))
            w = lowest_colored_ancestor[lcas[j]]
            
            # lca(B) < u and there is a colored vertex w with lca(B) < w <= u
            if (lca.ancestor_not_equal(u, lcas[j]) and
                w and 
                lca(w, u) is u):
                matrix[i][j] = 'essential'
            # lca(A) < lca(B)
            elif lca.ancestor_not_equal(lcas[j], lcas[i]):
                matrix[i][j] = 'forbidden'
            else:
                matrix[i][j] = 'ambiguous'
    
    return matrix


def _find_child(u, v, lca):
    """For v < u, find the child w of u such that v <= w."""
    
    for w in u.children:
        if lca(v, w) is w:
            return w
    
    raise RuntimeError('could not find corresponding child of u '\
                       '(u not an ancestor of v?)')


def fitch_orientation_for_refinements(T, partition, lca=None):
    """(Un)ambiguously present or absent edges in the directed Fitch graph of
    all compatible refinements.
    
    The ordered pair (x, y) of two leaves x and y is in the Fitch relation
    if the path from lca(x, y) to y contains a transfer edge. An undirected
    Fitch relation ((x, y) if (x, y) or (y, x) in the directed relation) can
    be represented by a partition (the maximal independent sets).
    In combination with a tree, this partition determines to some extent the
    directed Fitch relation w.r.t. the possible refinements of the tree.
    
    Parameters
    ----------
    T : Tree
        A tree with unique leaf label.
    partition : list or tuple
        A partition of the labels of the trees' leaves.
    lca : LCA, optional
        Precomputed LCA datastructure for the tree. The default is None in
        which case it is computed internally.
    
    Returns
    -------
    bool
        Returns False if the tree and the partition are not refinement-
        compatible.
    2-dimensional array (list of lists)
        Stores at indices i and j corresponding to sets A and B, resp., of the
        partition whether (A, B) is 'essential', 'forbidden', or 'ambiguous'
        for all refinements of the tree that ae compatible with the partition.
    """
    
    if not isinstance(T, Tree):
        raise TypeError("T must be of type 'Tree'")
    
    if not isinstance(partition, (list, tuple)):
        raise TypeError("partition must be of type 'list' or 'tuple'")
    
    if not lca:
        lca = LCA(T)
    
    r_compatibility = is_refinement_compatible(T, partition, lca=lca)
    if r_compatibility:
        edge_coloring, lcas = r_compatibility
    else:
        return False
    
    # compute for each vertex its lowest colored (strict) ancestor
    lowest_colored_edge = {}
    for v in T.preorder():
        if not v.parent:
            lowest_colored_edge[v] = None
        elif v in edge_coloring:
            lowest_colored_edge[v] = v
        else:
            lowest_colored_edge[v] = lowest_colored_edge[v.parent]
    
    n = len(partition)
    matrix = [[None for j in range(n)] for i in range(n)]
    
    for i in range(n):
        for j in range(n):
            
            if i == j:
                matrix[i][j] = 'forbidden'
                continue
            
            u = lca(lcas[i], lcas[j])       # lca(lca(A), lca(B))
            v_A = lcas[i]                   # lca(A)
            v_B = lcas[j]                   # lca(B)
            lce_B = lowest_colored_edge[v_B]
            
            if (lca.ancestor_not_equal(u, v_B) and
                lce_B and lca.ancestor_not_equal(u, lce_B)):
                matrix[i][j] = 'essential'
            elif (u.parent and
                  lca.ancestor_not_equal(u, v_A) and
                  u in edge_coloring and
                  edge_coloring[u] == \
                      edge_coloring.get(_find_child(u, v_A, lca))):
                matrix[i][j] = 'essential'
            elif (lca.ancestor_not_equal(v_B, v_A) and
                  edge_coloring.get(_find_child(u, v_A, lca)) == j):
                matrix[i][j] = 'forbidden'
            else:
                matrix[i][j] = 'ambiguous'
    
    return matrix
    
