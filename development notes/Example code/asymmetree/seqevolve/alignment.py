"""Module for constructing the true alignment of simulated sequences.

In the simulated sequences, the true history of substitutions and indels is known. This module uses
this information to construct the true alignment of the sequences, which can for instance be used
for evaluating the performance of alignment methods.
"""

from __future__ import annotations

import numpy as np
import networkx as nx
from tralda.datastructures import Tree
from tralda.datastructures import TreeNode

from asymmetree.seqevolve.evolving_sequence import EvoSeq


class AlignmentBuilder:
    """Construction of the true alignment of simulated sequences."""

    def __init__(
        self,
        tree: Tree,
        sequence_dict: dict[TreeNode, EvoSeq],
        alphabet: str,
        include_inner: bool = True,
    ) -> None:
        """Constructor for the AlignmentBuilder class.

        Args:
            tree: The tree along which the sequences where simulated.
            sequence_dict: A dict containing the TreeNode instances in the tree as keys and the
                simulated sequences as values.
            alphabet: The alphabet of the substitution model that was used.
            include_inner: If True, the alignment also contains the sequences of the inner nodes of
                the tree.
        """
        self.tree = tree
        self.sequence_dict = sequence_dict
        self.include_inner = include_inner

        # add gap character to the alphabet if not already present
        self.alphabet = alphabet if alphabet[-1] == "-" else alphabet + "-"

    def build(self) -> dict[TreeNode, str]:
        """Build the true alignment.

        Returns:
            A `dict` containing the `TreeNode` instances in the tree as keys and the str sequences
                as values that include the necessary gaps.
        """
        self._get_preorder()
        self._sort_sites()
        self._alignment_matrix()

        return self._sequences()

    def _get_preorder(self) -> None:
        """Get the preorder of the tree."""
        self.nodes = [v for v in self.tree.preorder() if self.include_inner or not v.children]
        self.node_index = {item: index for index, item in enumerate(self.nodes)}

    def _sort_sites(self) -> None:
        """Sort the sites of the sequences in a way that they can be correctly aligned."""
        G = nx.DiGraph()

        for node in self.nodes:
            seq = self.sequence_dict[node]

            # handle case |seq| = 1 correctly
            if len(seq) == 1:
                G.add_node(seq[0].site_id)

            # add nodes and edges to auxiliary graph
            for first, second in seq.element_pairs():
                G.add_edge(first.site_id, second.site_id)

        self.sites = list(nx.topological_sort(G))

    def _alignment_matrix(self) -> None:
        """Construct the alignment matrix."""
        self.alignment = np.zeros((len(self.nodes), len(self.sites)), dtype=np.int8)
        positions = [self.sequence_dict[v]._first for v in self.nodes]

        for j in range(len(self.sites)):
            for i in range(len(positions)):
                if positions[i] and self.sites[j] == positions[i].site_id:
                    self.alignment[i, j] = positions[i]._value
                    positions[i] = positions[i]._next
                else:
                    self.alignment[i, j] = -1

    def _sequences(self) -> dict[TreeNode, str]:
        """Convert the alignment matrix to a dictionary of sequences.

        Returns:
            A `dict` containing the `TreeNode` instances in the tree as keys and the str sequences
                as values that include the necessary gaps.
        """
        result = {}

        for i in range(len(self.nodes)):
            sequence = "".join(self.alphabet[x] for x in self.alignment[i, :])
            result[self.nodes[i]] = sequence

        return result
