"""Orthology graph, BMG, RBMG.

This module provides classes concerning colored graphs and colored (phylogenetic) trees, including
a BMG(Tree)-function, and computation of anLRT from a (not necessarily valid) BMG or from a
(pruned) gene tree.

Phylogenetic best matches of a gene x of species X are defined as those genes y of another species
Y that share the lowest common ancestor with x in the gene tree among all genes in that species.
In contrast, two genes are orthologs if their last common ancestor was a speciation event.
Orthology and reciprocal best matches are closely related. The module asymmetree.analysis.
BestMatches contains methods for extracting BMGs from trees, least resolved trees of BMG, and
orthology graphs from trees.
"""

from __future__ import annotations

import itertools

import networkx as nx
import numpy as np

from tralda.datastructures import Tree
from tralda.datastructures import TreeNode
from tralda.utils.graph_tools import graphs_equal
from tralda.utils.graph_tools import is_properly_colored
from tralda.utils.graph_tools import sort_by_colors
from tralda.supertree import Build

from asymmetree.utils.phylogenetic_trees import reconstruct_reconc_from_graph
from asymmetree.utils.phylogenetic_trees import topology_only


# --------------------------------------------------------------------------------------------------
#                            True Orthology Graph (TOG), BMG & RBMG
#
#                                    (given a known tree)
# --------------------------------------------------------------------------------------------------


def orthology_from_tree(tree: Tree) -> nx.Graph:
    """Constructs the true orthology graph from an event-labeled tree.

    Args:
        tree: A tree with leaves that have the `label`, and inner vertices that have the `event`
            attribute set.

    Returns:
        The orthology graph explained by the tree, i.e., a graph on the set of leaf labels and edges
        between two of them if and only if their last common ancestor was a speciation event.
    """
    if not isinstance(tree, Tree):
        raise TypeError("not of type 'Tree'")

    leaves = tree.leaf_dict()
    TOG = nx.Graph()

    for v in leaves[tree.root]:
        TOG.add_node(v.label, color=v.reconc)

    for node in tree.preorder():
        if node.event != "S":
            continue

        for c1, c2 in itertools.combinations(node.children, 2):
            for u in leaves[c1]:
                for v in leaves[c2]:
                    TOG.add_edge(u.label, v.label)

    return TOG


def bmg_from_tree(
    tree: Tree,
    supply_rbmg: bool = False,
) -> nx.DiGraph | tuple[nx.DiGraph, nx.Graph]:
    """Construct a BMG (and optionally RBMG) from a given tree.

    Args:
        tree: A tree with leaves that have the `label` and `reconc` attribute set.
        supply_rbmg: If True, also return the symmetric part of the constructed BMG.

    Returns:
        The constructed BMG and optionally its symmetric part.

    References:
        1. M. Geiß, E. Chávez, M. González Laffitte, A. López Sánchez, B. M. R. Stadler, D. I.
           Valdivia, M. Hellmuth, M. Hernández Rosales, and P. F.Stadler. Best match graphs.
           In: Journal of Mathematical Biology, 78:2015-2057, 2019. doi:10.1007/s00285-019-01332-9.
    """
    if not isinstance(tree, Tree):
        raise TypeError("not of type 'Tree'")

    leaves = tree.leaf_dict()
    bmg = nx.DiGraph()
    colors = set()

    for v in leaves[tree.root]:
        colors.add(v.reconc)
        bmg.add_node(v.label, color=v.reconc)

    for u in leaves[tree.root]:
        remaining = colors - set([u.reconc])  # colors to which no best match has yet been found
        parent = u.parent  # start with direct parent of each node

        while remaining and parent:
            colors_here = set()
            for v in leaves[parent]:
                if v.reconc in remaining:  # best match found
                    colors_here.add(v.reconc)
                    bmg.add_edge(u.label, v.label)
            remaining -= colors_here
            parent = parent.parent

    if not supply_rbmg:
        return bmg
    else:
        return bmg, bmg.to_undirected(reciprocal=True)


def extended_best_hits(
    leaves: list[TreeNode],
    D: np.ndarray,
    epsilon: float = 1e-8,
    supply_rbmg: bool = False,
) -> nx.DiGraph | tuple[nx.DiGraph, nx.Graph]:
    """Compute BMG and RBMG from a distances matrix D.

    Apporximation of best matches as the genes within a tolerance range from the gene with the
    smallest distance, see [1].

    Args:
        leaves: List of TreeNode instances with `label` attribute.
        D: The distances of the leaves where the order of the rows/columns corresponds to that in
            the list.
        epsilon: Relative threshold for the inclusion of suboptimal hits, i.e., (x,y) in BMG if
            D(x,y) <= (1 + epsilon) * min d(x,y').
        supply_rbmg: If True, also return the symmetric part of the constructed BMG.

    Returns:
        The constructed BMG and optionally its symmetric part.

    References:
        1. P. F. Stadler, M. Geiß, D. Schaller, A. López Sánchez, M. Gonzalez Laffitte, D. I.
           Valdivia, M. Hellmuth, and M. Hernández Rosales. From pairs of most similar sequences to
           phylogenetic best matches. In: Algorithms for Molecular Biology, 15:5, 2020.
           doi: 10.1186/s13015-020-00165-2.
    """
    bmg = nx.DiGraph()
    colors = set()
    relative_threshold = 1 + epsilon

    for v in leaves:
        bmg.add_node(v.label, color=v.reconc)
        colors.add(v.reconc)

    # build BMG
    for i, u in enumerate(leaves):
        minima = {c: float("inf") for c in colors}
        for j, v in enumerate(leaves):
            if D[i, j] < minima[v.reconc]:
                minima[v.reconc] = D[i, j]

        for j, v in enumerate(leaves):
            if u.reconc != v.reconc and D[i, j] <= relative_threshold * minima[v.reconc]:
                bmg.add_edge(u.label, v.label, distance=D[i, j])

    if not supply_rbmg:
        return bmg
    else:
        return bmg, bmg.to_undirected(reciprocal=True)


# --------------------------------------------------------------------------------------------------
#                                   Least resolved trees
# --------------------------------------------------------------------------------------------------


def informative_triples(
    graph: nx.DiGraph,
    color_dict: dict | None = None,
) -> list[tuple]:
    """Compute the informative triples of a colored digraph.

    Args:
        graph: A digraph whose nodes have the 'color' attribute.
        color_dict: A dict containing the colors that appear in the graph as keys and the list of
            nodes with the respective color as values. If not provided, this dict is computed
            internally.

    Returns:
        A list where each item is a tuple (a, b1, b2) of three nodes such that a b1 | b2 is an
        informative triple in the graph, see [1].

    References:
        1. M. Geiß, E. Chávez, M. González Laffitte, A. López Sánchez, B. M. R. Stadler, D. I.
           Valdivia, M. Hellmuth, M. Hernández Rosales, and P. F. Stadler. Best match graphs.
           In: Journal of Mathematical Biology, 78:2015-2057, 2019. doi:10.1007/s00285-019-01332-9.
    """
    if not isinstance(graph, nx.DiGraph):
        raise TypeError("must be a NetworkX 'Digraph'")

    if not color_dict:
        color_dict = sort_by_colors(graph)

    R = []

    for c1, c2 in itertools.permutations(color_dict.keys(), 2):
        for a in color_dict[c1]:
            for b1, b2 in itertools.permutations(color_dict[c2], 2):
                if graph.has_edge(a, b1) and (not graph.has_edge(a, b2)):
                    R.append((a, b1, b2))

    return R


def forbidden_triples(graph: nx.DiGraph, color_dict: dict | None = None) -> list[tuple]:
    """Compute the forbidden triples of a colored digraph.

    Args:
        graph: A digraph whose nodes have the 'color' attribute.
        color_dict: A dict containing the colors that appear in the graph as keys and the list of
            nodes with the respective color as values. If not provided, this dict is computed
            internally.

    Returns:
        A list where each item is a tuple (a, b1, b2) of three nodes such that a b1 | b2 is a
        forbidden triple in the graph, see [1].

    References:
        1. D. Schaller, P. F. Stadler, and M. Hellmuth. Complexity of modification problems for
           best match graphs. In: Theoretical Computer Science, 865:63-84, 2021.
           doi:10.1016/j.tcs.2021.02.037.
    """
    if not isinstance(graph, nx.DiGraph):
        raise TypeError("must be a NetworkX 'Digraph'")

    if not color_dict:
        color_dict = sort_by_colors(graph)

    F = []

    for c1, c2 in itertools.permutations(color_dict.keys(), 2):
        for a in color_dict[c1]:
            for b1, b2 in itertools.combinations(color_dict[c2], 2):
                if graph.has_edge(a, b1) and graph.has_edge(a, b2):
                    F.append((a, b1, b2))
                    F.append((a, b2, b1))

    return F


def informative_forbidden_triples(
    graph: nx.DiGraph, color_dict: dict | None = None
) -> tuple[list[tuple], list[tuple]]:
    """Compute the informative and forbidden triples of a colored digraph.

    Args:
        graph: A digraph whose nodes have the 'color' attribute.
        color_dict: A dict containing the colors that appear in the graph as keys and the list of
            nodes with the respective color as values. If not provided, this dict is computed
            internally.

    Returns:
        Tuple of two lists. Each item in the first (resp. second) list is a tuple (a, b1, b2) of
        three nodes such that a b1 | b2 is a informative (resp. forbidden) triple in the graph, see
        [1].

    References:
        1. D. Schaller, P. F. Stadler, and M. Hellmuth. Complexity of modification problems for
           best match graphs. In: Theoretical Computer Science, 865:63-84, 2021.
           doi:10.1016/j.tcs.2021.02.037.
    """
    if not isinstance(graph, nx.DiGraph):
        raise TypeError("must be a NetworkX 'Digraph'")

    if not color_dict:
        color_dict = sort_by_colors(graph)

    R, F = [], []

    for c1, c2 in itertools.permutations(color_dict.keys(), 2):
        for a in color_dict[c1]:
            for b1, b2 in itertools.combinations(color_dict[c2], 2):
                if graph.has_edge(a, b1) and graph.has_edge(a, b2):
                    F.append((a, b1, b2))
                    F.append((a, b2, b1))
                elif graph.has_edge(a, b1):
                    R.append((a, b1, b2))
                elif graph.has_edge(a, b2):
                    R.append((a, b2, b1))

    return R, F


def binary_explainable_triples(
    graph: nx.DiGraph,
    color_dict: dict | None = None,
) -> list[tuple]:
    """Extended informative triple set for binary-explainable graphs.

    Args:
        graph: A digraph whose nodes have the 'color' attribute.
        color_dict: A dict containing the colors that appear in the graph as keys and the list of
            nodes with the respective color as values. If not provided, this dict is computed
            internally.

    Returns:
        List of tuples where each item is a tuple (a, b, c) of three nodes such that ab|c is a
        triple displayed by every binary tree that explains the graph (if existent), see [1].

    References:
        1. D. Schaller, M. Geiß, M. Hellmuth, and P. F. Stadler. Best match graphs with binary
           trees. In: C. Martín-Vide, M. A. Vega-Rodríguez, and T. Wheeler, editors, Algorithms for
           Computational Biology, 8th AlCoB, volume 12715 of Lecture Notes in Computer Science,
           pages 82-93, 2021. doi: 10.1007/978-3-030-74432-8_6.
    """
    if not isinstance(graph, nx.DiGraph):
        raise TypeError("must be a NetworkX 'Digraph'")

    if not color_dict:
        color_dict = sort_by_colors(graph)

    R_binary = []

    for c1, c2 in itertools.permutations(color_dict.keys(), 2):
        for a in color_dict[c1]:
            for b1, b2 in itertools.combinations(color_dict[c2], 2):
                if graph.has_edge(a, b1) and graph.has_edge(a, b2):
                    try:
                        if b2 < b1:
                            b1, b2 = b2, b1
                    except NotImplemented:
                        pass
                    R_binary.append((b1, b2, a))
                elif graph.has_edge(a, b1):
                    R_binary.append((a, b1, b2))
                elif graph.has_edge(a, b2):
                    R_binary.append((a, b2, b1))

    return R_binary


def _finalize(tree: Tree | TreeNode | None, G: nx.Graph | nx.DiGraph) -> Tree | None:
    """Finalize the construction of a tree from a graph.

    This method is used in the construction of least resolved trees from BMGs and RBMGs. It performs
    the following steps: 1) If the input is a TreeNode, it is converted to a Tree.
    2) Reconciliations are assigned to the leaves based on the input graph.

    Args:
        tree: A Tree, TreeNode, or None. If None, the method returns None.
        G: A graph used to assign reconciliations to the leaves of the tree.

    Returns:
        A Tree with reconciliations assigned to the leaves, or None if the input tree was None.
    """
    if not tree:
        return None

    if isinstance(tree, TreeNode):
        tree = Tree(tree)

    # assign reconciliationss to the leaves
    reconstruct_reconc_from_graph(tree, G)

    return tree


# --------------------------------------------------------------------------------------------------
#                      LRT FROM OBSERVABLE GENE TREE
# --------------------------------------------------------------------------------------------------


def lrt_from_tree(T: Tree) -> Tree:
    """Computes the Least Resolved Tree from a tree.

    The unique Least Resolved Tree from a (pruned) gene tree is computed by contraction of all
    redundant edges.

    Args:
        T: A tree whose leaf nodes have the 'reconc' attribute.

    Returns:
        The unique least resolved tree (w.r.t. best matches), see [1].

    References:
        1. D. Schaller, M. Geiß, P. F. Stadler, and M. Hellmuth. Complete characterization of
           incorrect orthology assignments in best match graphs. In: Journal of Mathematical
           Biology, 82:20, 2021. doi: 10.1007/s00285-021-01564-8.
    """
    lrt = T.copy()
    if not lrt.root:
        return lrt

    # remove planted root if existent
    if len(lrt.root.children) == 1:
        new_root = lrt.root.children[0]
        new_root.detach()
        lrt.root = new_root
        new_root.dist = 0.0

    # assign list of leaves to each node
    leaves = lrt.leaf_dict()

    subtree_reconcs = {}
    for v in lrt.preorder():
        subtree_reconcs[v] = {leaf.reconc for leaf in leaves[v]}

    arc_colors = _arc_colors(lrt, leaves, subtree_reconcs)
    redundant_edge_list = redundant_edges(lrt, subtree_reconcs, arc_colors)
    lrt.contract(redundant_edge_list)
    lrt = topology_only(lrt)

    return lrt


def _arc_colors(
    T: Tree,
    leaves: dict[TreeNode, list[TreeNode]],
    subtree_reconcs: dict[TreeNode, set],
) -> dict[TreeNode, set]:
    """Color sets relevant for redundant edge computation.

    Computes for all inner vertices v the color set of y such that y with (x,y) is an arc in the
    BMG and lca(x,y) = v.

    Args:
        T: A tree whose leaf nodes have the 'reconc' attribute.
        leaves: A dict mapping each node to the list of its leaf descendants.
        subtree_reconcs: A dict mapping each node to the set of reconciliations of its leaf
            descendants.

    Returns:
        A dict mapping each inner vertex to the set of colors s such that there is an arc (x,y) in
        the BMG with lca(x,y) = v and y has color s.
    """
    all_colors = subtree_reconcs[T.root]

    # color sets for all v
    arc_colors = {v: set() for v in T.preorder()}

    for u in leaves[T.root]:
        # colors to which no best match has yet been found
        remaining = all_colors - {u.reconc}

        # start with direct parent of each node
        current = u.parent

        while remaining and current:
            colors_here = set()
            for v in leaves[current]:
                # best match found
                if v.reconc in remaining:
                    colors_here.add(v.reconc)

            arc_colors[current].update(colors_here)
            remaining -= colors_here
            current = current.parent

    return arc_colors


def redundant_edges(
    T: Tree,
    subtree_reconcs: dict[TreeNode, set],
    arc_colors: dict[TreeNode, set],
) -> list[tuple[TreeNode, TreeNode]]:
    """Compute the redundant edges of a tree.

    An inner edge (u,v) is redundant if there is no arc (x,y) in the BMG with lca(x,y) = u and y
    has a color that is not in the subtree of v, see [1].

    Args:
        T: A tree whose leaf nodes have the 'reconc' attribute.
        subtree_reconcs: A dict mapping each node to the set of reconciliations of its leaf
            descendants.
        arc_colors: A dict mapping each inner vertex to the set of colors s such that there is an
            arc (x,y) in the BMG with lca(x,y) = v and y has color s.

    Returns:
        A list of tuples representing the redundant edges of the tree.
    """
    redundant_edge_list = []

    for u, v in T.inner_edges():
        # colors s in sigma( L(T(u) \ T(v)) )
        aux_set = set()

        for v2 in u.children:
            if v2 is not v:
                aux_set.update(subtree_reconcs[v2])

        if not arc_colors[v].intersection(aux_set):
            redundant_edge_list.append((u, v))

    return redundant_edge_list


# --------------------------------------------------------------------------------------------------
#                                       LRT FROM BMG
# --------------------------------------------------------------------------------------------------


def is_bmg(G: nx.DiGraph) -> Tree | None:
    """Determine whether a colored digraph is a BMG.

    If the graph is a BMG, then its LRT is returned, see [1].

    Args:
        G: A digraph whose nodes have the 'color' attribute.

    Returns:
        The least-resolved tree if the graph was a BMG, and None otherwise.

    References:
        1. D. Schaller, M. Geiß, E. Chávez, M. González Laffitte, A. López Sánchez, B. M. R.
           Stadler, D. I. Valdivia, M. Hellmuth, M. Hernández Rosales, and P. F. Stadler.
           Corrigendum to 'Best match graphs'. In: Journal of Mathematical Biology, 82:47, 2021.
           doi: 10.1007/s00285-021-01601-6.
    """
    if not isinstance(G, nx.DiGraph):
        raise TypeError("not a digraph")
    if not is_properly_colored(G):
        raise RuntimeError("not a properly colored digraph")

    subtrees = []
    colors = None

    for wcc in nx.weakly_connected_components(G):
        sg = G.subgraph(wcc).copy()
        color_dict = sort_by_colors(sg)

        # the color sets of all components must be equal
        if colors is None:
            colors = set(color_dict)
        elif colors != set(color_dict):
            return

        L = {v for v in sg.nodes()}
        R = informative_triples(sg, color_dict=color_dict)
        build = Build(R, L, mincut=False)
        subtree = build.build_tree()

        if not subtree:
            return
        else:
            # a digraph is a BMG iff its equal to the BMG of BUILD(R)
            reconstruct_reconc_from_graph(subtree, sg)
            if not graphs_equal(sg, bmg_from_tree(subtree)):
                return
            subtrees.append(subtree)

    if len(subtrees) == 1:
        root = subtrees[0].root
    else:
        root = TreeNode()
        for subtree in subtrees:
            root.add_child(subtree.root)

    return _finalize(root, G)


def lrt_from_colored_graph(
    G: nx.DiGraph,
    mincut: bool = False,
    weighted_mincut: bool = False,
) -> Tree | None:
    """Compute the least resolved tree from a colored digraph.

    If the graph is a BMG, then the unique least resolved tree is returned. If not, then a heuristic
    least resolved tree is returned if 'mincut' is True, and None otherwise.

    Args:
        G: A digraph whose nodes have the 'color' attribute.
        mincut: Handle inconsistencies that occur if the supplied graph is not a BMG.
        weighted_mincut: If True, apply a weighted mincut heuristic.

    Returns:
        The least resolved tree if the graph was a BMG, a heuristic least resolved tree if 'mincut'
        is True, and None otherwise.
    """
    L = {v for v in G.nodes()}
    R = informative_triples(G)

    build = Build(R, L, mincut=mincut, weighted_mincut=weighted_mincut)
    tree = build.build_tree()

    return _finalize(tree, G)


def correct_bmg(G: nx.DiGraph) -> nx.DiGraph:
    """Build the LRT (using min cut in BUILD algorithm) and return its BMG.

    Args:
        G: A digraph whose nodes have the 'color' attribute.

    Returns:
        The graph edited to be a BMG.

    Raises:
        RuntimeError: If the mincut heuristic fails to construct a tree.
    """
    if G.order() == 0:
        return nx.DiGraph()

    subtrees = []
    for sg in (G.subgraph(c) for c in nx.weakly_connected_components(G)):
        tree = lrt_from_colored_graph(sg, mincut=True)
        if not tree:
            raise RuntimeError("mincut heuristic failed to construct a tree")

        subtrees.append(tree)

    if len(subtrees) == 1:
        tree = subtrees[0]
    else:
        tree = Tree(TreeNode())
        for subtree in subtrees:
            tree.root.add_child(subtree.root)

    return bmg_from_tree(tree)


# --------------------------------------------------------------------------------------------------
#                                   LRT FROM 2-col. BMG
# --------------------------------------------------------------------------------------------------


class TwoColoredLRT:
    """Construct the least-resolved tree of a 2-colored best match graph.

    Implementation of the efficient algorithm in [1].

    References:
        1. D. Schaller, M. Geiß, M. Hellmuth, and P. F. Stadler. Least resolved trees for
           two-colored best match graphs. In: Journal of Graph Algorithms and Applications,
           25(1):397-416, 2021. doi: 10.7155/jgaa.00564.
    """

    def __init__(self, digraph: nx.DiGraph) -> None:
        """Construct the least-resolved tree of a 2-colored best match graph.

        Args:
            digraph: A digraph whose nodes have the 'color' attribute.
        """
        if not isinstance(digraph, nx.DiGraph):
            raise TypeError("input must be a digraph")

        self.digraph = digraph
        self.color_dict = sort_by_colors(digraph)

    def build_tree(self) -> Tree | None:
        """Construct the least-resolved tree.

        Returns:
            The least-resolved tree if the graph was a 2-colored BMG, and None otherwise.
        """
        if not is_properly_colored(self.digraph):
            raise RuntimeError("not a properly colored digraph")
        if len(self.color_dict) > 2:
            raise RuntimeError("more than 2 colors")

        # star tree if there is only one color
        if len(self.color_dict) == 1:
            root = TreeNode()
            for v in self.digraph.nodes():
                root.add_child(TreeNode(label=v))
        # 2 colors
        else:
            subtrees = []
            for wcc in nx.weakly_connected_components(self.digraph):
                if len(wcc) == 1:
                    return

                subroot = self._build_tree(self.digraph.subgraph(wcc).copy())

                if not subroot:
                    return
                else:
                    subtrees.append(subroot)

            if len(subtrees) == 1:
                root = subtrees[0]
            else:
                root = TreeNode()
                for subroot in subtrees:
                    root.add_child(subroot)

        return _finalize(root, self.digraph)

    def _build_tree(self, G: nx.DiGraph) -> TreeNode | None:
        """Recursive helper method for build_tree.

        Args:
            G: A digraph whose nodes have the 'color' attribute.

        Returns:
            A TreeNode representing the least-resolved tree of G if G is a 2-colored BMG, and None
            otherwise.
        """
        color_count = {c: 0 for c in self.color_dict.keys()}
        for v in G.nodes():
            color_count[G.nodes[v]["color"]] += 1
        other_color = {color: G.order() - count for color, count in color_count.items()}

        umbrella = {v for v in G.nodes() if other_color[G.nodes[v]["color"]] == G.out_degree(v)}
        S_1 = {v for v in umbrella if umbrella.issuperset(G.predecessors(v))}
        S_2 = {v for v in S_1 if S_1.issuperset(G.predecessors(v))}

        if not S_2 or len(S_1) != len(S_2):
            return

        node = TreeNode()
        for v in S_2:
            node.add_child(TreeNode(label=v))
            G.remove_node(v)

        for wcc in nx.weakly_connected_components(G):
            if len(wcc) == 1:
                return

            child = self._build_tree(G.subgraph(wcc).copy())

            if not child:
                return
            else:
                node.add_child(child)

        return node


def lrt_from_2bmg(G: nx.DiGraph) -> Tree | None:
    """Effieciently constructs the LRT for a 2-BMG via the support vertices.

    Implementation of the efficient algorithm in [1].

    Args:
        G: A digraph whose nodes have the 'color' attribute.

    Returns:
        The least-resolved tree if the graph was a 2-colored BMG, and None otherwise.

    References:
        1. D. Schaller, M. Geiß, M. Hellmuth, and P. F. Stadler. Least resolved trees for
           two-colored best match graphs. In: Journal of Graph Algorithms and Applications,
           25(1):397-416, 2021. doi: 10.7155/jgaa.00564.
    """
    builder = TwoColoredLRT(G)

    return builder.build_tree()


# --------------------------------------------------------------------------------------------------
#                               Binary-refinable tree (BRT)
# --------------------------------------------------------------------------------------------------


def binary_refinable_tree(
    G: nx.DiGraph,
    mincut: bool = False,
    weighted_mincut: bool = False,
) -> Tree | None:
    """Compute the binary-refinable tree (BRT) for a graph.

    The BRT of a BMG can, if it exists, be refined arbitrarily and such that it still explains the
    BMG. In particular, all binary explanations are refinements of the BRT.

    Args:
        G: A digraph whose nodes have the 'color' attribute.
        mincut: If True, handle inconsistencies that occur if the supplied graph is not a
            binary-explainable BMG.
        weighted_mincut: If True, apply a weighted mincut heuristic.

    Returns:
        The unique binary-resolvable tree if the graph is a binary-explainable BMG, see [1], or a
        heuristic tree if it is not but 'mincut' is True, or None if neither is the case.

    References:
        1. D. Schaller, M. Geiß, M. Hellmuth, and P. F. Stadler. Best match graphs with binary
           trees. In: C. Martín-Vide, M. A. Vega-Rodríguez, and T. Wheeler, editors, Algorithms for
           Computational Biology, 8th AlCoB, volume 12715 of Lecture Notes in Computer Science,
           pages 82-93, 2021. doi: 10.1007/978-3-030-74432-8_6.
    """
    L = {v for v in G.nodes()}
    R = binary_explainable_triples(G)

    build = Build(R, L, mincut=mincut, weighted_mincut=weighted_mincut)
    brt = build.build_tree()

    return _finalize(brt, G)


# --------------------------------------------------------------------------------------------------
#                                     Augmented tree
# --------------------------------------------------------------------------------------------------


def augment_and_label(tree: Tree, inplace: bool = False) -> Tree:
    """Augment tree and add event labeling based on color intersections.

    Args:
        tree: A tree whose nodes have the 'reconc' attribute.
        inplace: If True, the supplied tree instance is directly manipulated, otherwise the tree
            is copied first.

    Returns:
        The unique augmented and duplication/loss-labeled tree as defined in [1].

    References:
        1. D. Schaller, M. Geiß, P. F. Stadler, and M. Hellmuth. Complete characterization of
           incorrect orthology assignments in best match graphs. In: Journal of Mathematical
           Biology, 82:20, 2021. doi: 10.1007/s00285-021-01564-8.
    """
    if not inplace:
        tree = tree.copy()

    # precompute since tree will be modified
    inner_nodes = [u for u in tree.inner_nodes()]
    leaves = tree.leaf_dict()

    for u in inner_nodes:
        # should not appear except if tree is planted, otherwise ignore
        if len(u.children) == 1:
            continue

        C = _color_intersection_components(u, leaves)

        if len(C) == 1:
            u.event = "D"

        else:
            u.event = "S"

            for cc in C:
                if len(cc) > 1:
                    w = TreeNode(event="D")
                    u.add_child(w)
                    for v in cc:
                        u.remove_child(v)
                        w.add_child(v)

    return tree


def _color_intersection_components(
    u: TreeNode,
    leaves: dict[TreeNode, list[TreeNode]],
) -> list[list[TreeNode]]:
    """Compute the color intersection components for an inner vertex.

    For an inner vertex u, the color intersection components are the connected components of the
    auxiliary graph with vertex set the children of u and an edge between v and w if there is a
    color that appears in both the subtree of v and the subtree of w.

    Args:
        u: An inner vertex of a tree whose leaf nodes have the 'reconc' attribute.
        leaves: A dict mapping each node to the list of its leaf descendants.

    Returns:
        A list of lists of TreeNodes, where each inner list corresponds to a color intersection
        component.
    """
    aux_graph = nx.Graph()

    for v in u.children:
        aux_graph.add_node(v)

        for x in leaves[v]:
            aux_graph.add_edge(v, x.reconc)

    result = []

    for cc in nx.connected_components(aux_graph):
        # discard the colors in the components
        result.append([v for v in cc if isinstance(v, TreeNode)])

    return result


# --------------------------------------------------------------------------------------------------
#                         Analysis of the structure of (R)BMGs
# --------------------------------------------------------------------------------------------------


def classify_good_ugly(bmg: nx.DiGraph, rbmg: nx.Graph, fp: nx.Graph) -> nx.Graph:
    """Classify edges in the false-positives graph as middle_in_good, first_in_ugly, first_in_bad.

    The false positives in an RBMG (w.r.t. to orthology) can be classified into several categories,
    including middle edges of good quartets, first edges of ugly quartets, and first edges of bad
    quartets, see [1].

    Args:
        bmg: The best match graph.
        rbmg: The reciprocal best match graph.
        fp: The false-positives graph.

    Returns:
        The updated false-positives graph with edges classified as middle_in_good, first_in_ugly,
        and first_in_bad.

    References:
        1. D. Schaller, M. Geiß, P. F. Stadler, and M. Hellmuth. Complete characterization of
           incorrect orthology assignments in best match graphs. In: Journal of Mathematical
           Biology, 82:20, 2021. doi: 10.1007/s00285-021-01564-8.
    """
    for x, y in fp.edges():
        fp.edges[x, y]["middle_in_good"] = False
        fp.edges[x, y]["first_in_ugly"] = False
        fp.edges[x, y]["first_in_bad"] = False

    color_dict = sort_by_colors(bmg)

    for x, y in fp.edges():
        # middle edges of good quartets
        for c in color_dict:
            if c == bmg.nodes[x]["color"] or c == bmg.nodes[y]["color"]:
                continue

            for z1, z2 in itertools.permutations(color_dict[c], 2):
                if (
                    rbmg.has_edge(z1, x)
                    and rbmg.has_edge(z2, y)
                    and bmg.has_edge(z1, y)
                    and not bmg.has_edge(y, z1)
                    and bmg.has_edge(z2, x)
                    and not bmg.has_edge(x, z2)
                ):
                    fp.edges[x, y]["middle_in_good"] = True

        # first edges of ugly/bad quartets
        for _ in range(2):
            for x2 in color_dict[bmg.nodes[x]["color"]]:
                if x == x2:
                    continue

                for z in rbmg.neighbors(x2):
                    if (
                        bmg.nodes[z]["color"] == bmg.nodes[x]["color"]
                        or bmg.nodes[z]["color"] == bmg.nodes[y]["color"]
                    ):
                        continue
                    # first edges of ugly quartets
                    if rbmg.has_edge(y, x2) and not rbmg.has_edge(z, y) and not rbmg.has_edge(z, x):
                        fp.edges[x, y]["first_in_ugly"] = True

                    # first edges of bad quartets
                    elif (
                        not rbmg.has_edge(y, x2) and rbmg.has_edge(z, y) and not bmg.has_edge(x, z)
                    ):
                        fp.edges[x, y]["first_in_bad"] = True

            # swap for second iteration
            x, y = y, x

    return fp


def count_good_ugly(bmg: nx.DiGraph, rbmg: nx.Graph, fp: nx.Graph) -> tuple[int, int, int]:
    """Count number of edge types in the false-positives graph.

    Args:
        bmg: The best match graph.
        rbmg: The reciprocal best match graph.
        fp: The false-positives graph.

    Returns:
        A tuple containing the number of edges in the false-positives graph that are middle edges of
        good quartets, first edges of ugly quartets, and both.
    """
    classify_good_ugly(bmg, rbmg, fp)

    good, ugly, both = 0, 0, 0

    for x, y in fp.edges():
        if fp.edges[x, y]["middle_in_good"] and fp.edges[x, y]["first_in_ugly"]:
            good += 1
            ugly += 1
            both += 1
        elif fp.edges[x, y]["middle_in_good"]:
            good += 1
        elif fp.edges[x, y]["first_in_ugly"]:
            ugly += 1

    return good, ugly, both


def count_good_ugly_bad(
    bmg: nx.DiGraph,
    rbmg: nx.Graph,
    fp: nx.Graph,
) -> tuple[int, int, int, int, int, int, int]:
    """Count number of edge types in the false-positives graph.

    Args:
        bmg: The best match graph.
        rbmg: The reciprocal best match graph.
        fp: The false-positives graph.

    Returns:
        A tuple containing the number of edges in the false-positives graph that are middle edges of
        good quartets, first edges of ugly quartets, first edges of bad quartets, and all
        combinations thereof in in order: (gub, gu, gb, ub, g, u, b).
    """
    classify_good_ugly(bmg, rbmg, fp)

    gub, gu, gb, ub, g, u, b = 0, 0, 0, 0, 0, 0, 0

    for x, y in fp.edges():
        if (
            fp.edges[x, y]["middle_in_good"]
            and fp.edges[x, y]["first_in_ugly"]
            and fp.edges[x, y]["first_in_bad"]
        ):
            gub += 1
            gu += 1
            gb += 1
            ub += 1
            g += 1
            u += 1
            b += 1

        elif fp.edges[x, y]["middle_in_good"] and fp.edges[x, y]["first_in_ugly"]:
            gu += 1
            g += 1
            u += 1
        elif fp.edges[x, y]["middle_in_good"] and fp.edges[x, y]["first_in_bad"]:
            gb += 1
            g += 1
            b += 1
        elif fp.edges[x, y]["first_in_ugly"] and fp.edges[x, y]["first_in_bad"]:
            ub += 1
            u += 1
            b += 1
        elif fp.edges[x, y]["middle_in_good"]:
            g += 1
        elif fp.edges[x, y]["first_in_ugly"]:
            u += 1
        elif fp.edges[x, y]["first_in_bad"]:
            b += 1

    return gub, gu, gb, ub, g, u, b


# --------------------------------------------------------------------------------------------------
#                                          P4 EDITING
# --------------------------------------------------------------------------------------------------


def list_path(t: int, k: int, path: list, G: nx.Graph, P4_list: list) -> None:
    """List all paths of length k in G.

    The paths are stored in P4_list as tuples of the form (v1, v2, v3, v4).

    Args:
        t: The current length of the path.
        k: The desired length of the path.
        path: The current path being explored.
        G: The graph in which to find paths.
        P4_list: A list to store the found paths.
    """
    for u in G[path[-1]]:
        G.nodes[u]["vis"] += 1

    for u in G[path[-1]]:
        if G.nodes[u]["vis"] == 1:
            if t < k:
                list_path(t + 1, k, path + [u], G, P4_list)
            elif path[0] < u:
                P4_list.append(tuple(path) + (u,))

    for u in G[path[-1]]:
        G.nodes[u]["vis"] -= 1


def find_all_P4(G: nx.Graph) -> list[tuple]:
    """Find all induced P4 in G.

    Args:
        G: A graph in which to find induced P4.

    Returns:
        A list of tuples of the form (v1, v2, v3, v4) representing all induced P4 in G.
    """
    P4_list = []

    for v in G:
        G.nodes[v]["vis"] = 0

    for v in G:
        G.nodes[v]["vis"] = 1
        list_path(2, 4, [v], G, P4_list)
        G.nodes[v]["vis"] = 0

    return P4_list


def is_good_quartet(path: tuple, bmg: nx.Graph) -> bool:
    """Determine whether an induced P4 is a good quartet in the bmg.

    Args:
        path: A tuple representing the induced P4.
        bmg: The best match graph.

    Returns:
        True if the induced P4 is a good quartet, False otherwise.
    """
    if (
        bmg.has_edge(path[0], path[2])
        and bmg.has_edge(path[3], path[1])
        and bmg.nodes[path[0]]["color"] == bmg.nodes[path[3]]["color"]
    ):
        return True
    else:
        return False


def remove_all_P4(rbmg: nx.Graph, bmg: nx.Graph, P4_list: list[tuple] = None) -> nx.Graph:
    """Find all good quartets and removes the middle edges.

    Args:
        rbmg: The reciprocal best match graph.
        bmg: The best match graph.
        P4_list: A list of tuples representing all induced P4 in the rbmg. If not provided, it is
            computed internally.

    Returns:
        The graph obtained from the RBMG by removing the middle edges of all good quartets.
    """
    graph = rbmg.copy()

    if P4_list is None:
        P4_list = find_all_P4(graph)

    for path in P4_list:
        if is_good_quartet(path, bmg) and graph.has_edge(path[1], path[2]):
            graph.remove_edge(path[1], path[2])

    return graph


# --------------------------------------------------------------------------------------------------
#                                         C6 DETECTION
# --------------------------------------------------------------------------------------------------


def find_all_C6(G: nx.Graph, P4_list: list[tuple] | None = None) -> list[tuple]:
    """Find induces C6 of the form <x1 y1 z1 x2 y2 z2>.

    Args:
        G: A graph in which to find induced C6.
        P4_list: A list of tuples representing all induced P4 in G. If not provided, it is computed
            internally.

    Returns:
        A list of tuples of the form (x1, y1, z1, x2, y2, z2) representing all induced C6 in G.

    Raises:
        RuntimeError: If a detected C6 does not satisfy the properties of an induced C6 in G.
    """
    if P4_list is None:
        P4_list = find_all_P4(G)

    endpoints: dict[tuple, dict[tuple, list[tuple]]] = {}
    for path in P4_list:
        a, b, c, d = path if path[0] < path[3] else path[::-1]
        a_col, b_col, c_col, d_col = (
            G.nodes[a]["color"],
            G.nodes[b]["color"],
            G.nodes[c]["color"],
            G.nodes[d]["color"],
        )

        # avoid multiple recording for the 3 different colors
        if a_col == d_col and a_col < b_col and a_col < c_col:
            if (a, d) not in endpoints:
                endpoints[(a, d)] = {}
            if (b_col, c_col) not in endpoints[(a, d)]:
                endpoints[(a, d)][(b_col, c_col)] = []
            endpoints[(a, d)][(b_col, c_col)].append(path)

    C6_list = []
    for colors in endpoints.values():
        for color_pair, paths in colors.items():
            y_color, z_color = color_pair

            # avoid multiple recording for the 2 directions:
            if y_color >= z_color or (z_color, y_color) not in colors:
                continue

            for path1 in paths:
                for path2 in colors[(z_color, y_color)]:
                    if (not G.has_edge(path1[1], path2[1])) and (
                        not G.has_edge(path1[2], path2[2])
                    ):
                        C6_list.append(path1 + (path2[2], path2[1]))
                        if not _validate_C6(G, C6_list[-1]):
                            raise RuntimeError("problem in C6 detection")

    return C6_list


def graph_type(G: nx.Graph) -> tuple[str, list[tuple], list[tuple]]:
    """Classify the graph as type A, B, or C based on the presence of induced P4 and C6.

    A graph is of type A if it contains no induced P4, of type B if it contains an induced P4 but no
    induced C6, and of type C if it contains an induced C6.

    Args:
        G: A graph in which to classify the structure.

    Returns:
        A tuple containing the type of the graph ("A", "B", or "C"), a list of all induced P4 in the
        graph, and a list of all induced C6 in the graph.
    """
    P4_list = find_all_P4(G)
    C6_list = find_all_C6(G, P4_list=P4_list)

    if C6_list:
        return "C", P4_list, C6_list
    elif P4_list:
        return "B", P4_list, C6_list
    else:
        return "A", P4_list, C6_list


def _validate_C6(G: nx.Graph, C6_path: list) -> bool:
    """Validate whether a given path forms a C6 in the graph.

    Args:
        G: The graph in which to validate the C6.
        C6_path: A list of nodes representing the path to validate.

    Returns:
        True if the path forms a C6, False otherwise.
    """
    C6 = G.subgraph(C6_path)
    if C6.size() != 6:
        return False

    for n in C6.nodes():
        if not C6.degree(n) == 2:
            return False

        for neighbor in G.neighbors(n):
            if G.nodes[n]["color"] == G.nodes[neighbor]["color"]:
                return False

    return True
