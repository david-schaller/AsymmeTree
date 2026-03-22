"""Evolve sequences along a tree using models for substitutions, indels, and rate heterogeneity.

References:
    .. [1] Z. Yang. Computational molecular evolution. Oxford series in ecology and evolution.
       Oxford University Press, 2006. ISBN 978-0-19-856699-1 978-0-19-856702-8.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from tralda.datastructures import Tree
from tralda.datastructures import TreeNode

from asymmetree.io.sequences import write_alignment
from asymmetree.io.sequences import write_fasta
from asymmetree.seqevolve.alignment import AlignmentBuilder
from asymmetree.seqevolve.evolving_sequence import EvoSeq
from asymmetree.seqevolve.evolving_sequence import State
from asymmetree.seqevolve.indels import IndelModel
from asymmetree.seqevolve.heterogeneity import HetModel
from asymmetree.seqevolve.substitution import SubstModel


class Evolver:
    """Class for evolving sequences along a tree."""

    def __init__(
        self,
        subst_model: SubstModel,
        indel_model: IndelModel | None = None,
        het_model: HetModel | None = None,
        gillespie: str | bool = "auto",
        **kwargs,
    ) -> None:
        """Constructor for the Evolver class.

        Args:
            subst_model: Substitution model.
            indel_model: Model for insertions and deletions. If not provided, no indels are
                simulated.
            het_model: Model for substitution rate heterogeneity. If not provided, no rate
                heterogeneity is simulated.
            gillespie: If True, the Gillespie algorithm is used instead of constructing the exchange
                probability matrix P via matrix diagonalization. The Gillespie algorithm may be
                faster if the rate heterogeneity is large. The default is 'auto', in which case the
                exchange probability matrix is used unless the heterogeneity model has more than one
                class or is sitewise, in which case the Gillespie algorithm is used.
        """
        self.subst_model = subst_model

        self.indel_model = indel_model
        self.het_model = het_model

        if gillespie == "auto":
            if het_model and (het_model.classes > 1 or het_model.sitewise):
                self._gillespie = True
            else:
                self._gillespie = False
        else:
            self._gillespie = bool(gillespie)

    def evolve_along_tree(
        self,
        tree: Tree,
        start_length: int = 200,
        start_seq: str | None = None,
    ) -> dict[TreeNode, EvoSeq]:
        """Evolve sequence along a tree.

        If 'start_seq' is not supplied, then the sequence at the root is constructed randomly using
        the stationary probabilities of the specified substitution model.

        Args:
            tree: The tree along which sequences are simulated.
            start_length: The length of the sequence at the root. The default is 200. This parameter
                is ignored when 'start_seq' is provided.
            start_seq: The start sequence for the root of the tree. The default is None, in which
                case a sequence is constructed at random.

        Returns:
            A dictionary containing the TreeNode instances in the tree as keys and the simulated
            sequences as values (instances of type EvoSeq).
        """
        self.tree = tree
        self.site_counter = 0
        self.sequences = {}

        if start_seq is None:
            root_seq = self._random_sequence(start_length)
        else:
            root_seq = self._initialize_with_sequence(start_seq)

        if self.het_model:
            self.het_model.assign(root_seq)

        for v in tree.preorder():
            if v.parent is None:
                self.sequences[v] = root_seq
            else:
                self.sequences[v] = self._evolve(self.sequences[v.parent], v.dist)

        return self.sequences

    def get_sequences(self, include_inner: bool = True) -> dict[TreeNode, str]:
        """Return the simulated sequence as strings.

        Args:
            include_inner: If True, the sequences of the inner nodes of the tree are also included.

        Returns:
            A dictionary with the tree nodes as keys and the str sequence (without gaps) as values.
        """
        if include_inner:
            return {v: self.subst_model.to_sequence(seq) for v, seq in self.sequences.items()}

        return {
            v: self.subst_model.to_sequence(seq) for v, seq in self.sequences.items() if v.is_leaf()
        }

    def write_sequences(self, file_path: str | Path, include_inner: bool = True) -> None:
        """Write the simulated sequences into a file in fasta format.

        Args:
            file_path: The path and filename to the file into which the sequences shall be written.
            include_inner: If True, the sequences of the inner nodes of the tree are also written
                to the file.
        """
        write_fasta(file_path, self.get_sequences(include_inner=include_inner))

    def true_alignment(
        self,
        include_inner: bool = True,
        write_to: str | Path | None = None,
        alignment_format: str = "phylip",
    ) -> dict[TreeNode, str]:
        """Construct the 'true' alignment of the simulated sequences.

        The true alignment is constructed using the information about the history of substitutions
        and indels that is available for the simulated sequences.

        Args:
            include_inner: If True, the alignment also contains the sequences of the inner nodes of
                the tree.
            write_to: The path and filename to the file into which the alignment shall be written.
                The default is None, in which case nothing is written to a file.
            alignment_format: Format of the alignment that is written to a file. The default is
                'phylip'.

        Returns:
            A dictionary containing the TreeNode instances in the tree as keys and the str sequences
            as values that include the necessary gaps.
        """
        alg_builder = AlignmentBuilder(
            self.tree, self.sequences, self.subst_model.alphabet, include_inner=include_inner
        )

        alignment = alg_builder.build()

        if write_to is not None:
            write_alignment(write_to, alignment, alignment_format=alignment_format)

        return alignment

    def _random_positions(self, n: int) -> np.ndarray:
        """Draw n random positions using the stationary probabilities of the substitution model.

        Args:
            n: The number of random positions to draw.

        Returns:
            An array of length n containing the drawn positions as indices corresponding to the
            characters in the alphabet of the substitution model.
        """
        return np.random.choice(len(self.subst_model.alphabet), n, p=self.subst_model.freqs)

    def _random_sequence(self, n: int) -> EvoSeq:
        """Generate a random sequence using the stationary probabilities of the substitution model.

        Args:
            n: The length of the random sequence to generate.

        Returns:
            An instance of type EvoSeq containing the generated sequence.
        """
        seq = EvoSeq()

        for x in self._random_positions(n):
            seq.append(x, State.ROOT, self.site_counter)
            self.site_counter += 1

        return seq

    def _initialize_with_sequence(self, sequence: str) -> EvoSeq:
        """Initialize the sequence at the root of the tree with a given sequence.

        Args:
            sequence: The sequence to initialize the root with. This should be a string containing
                characters from the alphabet of the substitution model.

        Returns:
            An instance of type EvoSeq containing the initialized sequence.
        """
        seq = EvoSeq()

        for x in self.subst_model.to_indices(sequence):
            seq.append(x, State.ROOT, self.site_counter)
            self.site_counter += 1

        return seq

    def _evolve(self, parent_seq: EvoSeq, distance: float) -> EvoSeq:
        """Evolve a sequence along a branch of the tree.

        Args:
            parent_seq: The sequence at the parent node of the branch along which the sequence is
                evolved.
            distance: The evolutionary distance along the branch.

        Returns:
            An instance of type EvoSeq containing the evolved sequence.
        """
        child_seq = parent_seq.clone()

        if self.indel_model:
            self._generate_indels(child_seq, distance)

            if self.het_model:
                self.het_model.assign(child_seq)

        if not self._gillespie:
            # use transition prob. matrix
            self._substitute(child_seq, distance)
        else:
            # use Gillespie algorithm
            self._substitute_gillespie(child_seq, distance)

        return child_seq

    # ----------------------------------------------------------------------------------------------
    #                       Substitution by transition probability matrix
    # ----------------------------------------------------------------------------------------------

    def _substitute(self, sequence: EvoSeq, t: float) -> None:
        """Substitute characters in a sequence using the transition probability matrix.

        Args:
            sequence: The sequence to be evolved.
            t: The evolutionary time.
        """
        # (cumulative) transition probability matrices
        P = {}

        r = np.random.random(len(sequence))

        for i, site in enumerate(sequence):
            if site.status == State.INHERITED:
                # mutate according to matrix of the corresponding class

                d = site.rate_factor * t

                if d not in P:
                    if d > 0.0:
                        P[d] = np.cumsum(self.subst_model.transition_prob_matrix(d), axis=1)
                    else:
                        continue

                site._value = np.argmax(P[d][site._value, :] > r[i])

            else:
                # choose random character
                site._value = np.argmax(self.subst_model.freqs_cumulative > r[i])

    # ----------------------------------------------------------------------------------------------
    #                            Substitution by Gillespie algorithm
    # ----------------------------------------------------------------------------------------------

    def _substitute_gillespie(self, sequence: EvoSeq, t: float) -> None:
        """Substitute characters in a sequence using the Gillespie algorithm.

        Args:
            sequence: The sequence to be evolved.
            t: The evolutionary time.
        """
        for site in sequence:
            # choose random characters for insertions
            if site.status == State.INSERTION:
                site._value = np.argmax(self.subst_model.freqs_cumulative > np.random.random())

            # Gillespie algorithm for inherited positions
            else:
                current_time = 0.0
                while current_time < t:
                    # use negative value on the diagonal of the rate matrix Q
                    rate = -site.rate_factor * self.subst_model.Q[site._value, site._value]

                    current_time += np.random.exponential(1 / rate) if rate > 0.0 else float("inf")

                    if current_time >= t:
                        break

                    r_mutation = -np.random.random() * self.subst_model.Q[site._value, site._value]
                    current_sum = 0.0

                    for i in range(len(self.subst_model.alphabet)):
                        # skip the entry on diagonal
                        if i == site._value:
                            continue

                        current_sum += self.subst_model.Q[site._value, i]

                        if current_sum > r_mutation:
                            mutation = i
                            break

                    site._value = mutation

    # ----------------------------------------------------------------------------------------------
    #                                          Indels
    # ----------------------------------------------------------------------------------------------

    def _generate_indels(self, sequence: EvoSeq, t: float) -> None:
        """Generates indels using the Gillespie algorithm.

        Args:
            sequence: The sequence to be evolved.
            t: The evolutionary time.
        """
        current_time = 0.0

        while current_time < t:
            ins_rate, del_rate = self.indel_model.get_rates(len(sequence))
            total_rate = ins_rate + del_rate

            current_time += (
                np.random.exponential(1 / total_rate) if total_rate > 0.0 else float("inf")
            )

            if current_time >= t:
                break

            if np.random.random() < ins_rate / total_rate:
                self._insertion(sequence)
            else:
                self._deletion(sequence)

    def _insertion(self, sequence: EvoSeq) -> None:
        """Insert characters into a sequence.

        Args:
            sequence: The sequence to be evolved.
        """
        d = self.indel_model.draw_length()
        if d == 0:
            return

        pos = np.random.randint(-1, high=len(sequence))

        # initialize insertion before the first item
        if pos == -1:
            current_site = sequence.append_left(None, State.INSERTION, self.site_counter)
            self.site_counter += 1
            d -= 1

        # go to the element at position pos
        else:
            current_site = sequence.node_at(pos)

        for _ in range(d):
            current_site = sequence.insert_right_of(
                current_site, None, State.INSERTION, self.site_counter
            )

            self.site_counter += 1

    def _deletion(self, sequence: EvoSeq) -> None:
        """Delete characters from a sequence.

        Args:
            sequence: The sequence to be evolved.
        """
        d = self.indel_model.draw_length()
        if d == 0:
            return

        pos = np.random.randint(-d + 1, high=len(sequence))

        # deletion begins before or at the start of the sequence
        if pos <= 0:
            d += pos
            pos = 0

        sequence.remove_range(pos, d)
