"""Indel models for sequence evolution.

Indel events (insertions and deletions) are an important aspect of sequence evolution. This module
provides a class for modeling indel events, which can be used to simulate sequence evolution along a
tree. The class supports different distributions for the length of the indels.

References:
    1. Z. Yang. Computational molecular evolution. Oxford series in ecology and evolution.
       Oxford University Press, 2006. ISBN 978-0-19-856699-1 978-0-19-856702-8.
    2. R. A. Cartwright. DNA assembly with gaps (Dawg): Simulating sequence evolution.
       In: Bioinformatics, 21(Suppl 3):iii31-iii38, November 2005.
       doi:10.1093/bioinformatics/bti1200.
"""

from __future__ import annotations

from asymmetree.utils.sampling import Sampler


class IndelModel:
    """Indel model"""

    def __init__(
        self,
        insertion_rate: float,
        deletion_rate: float,
        length_distr: tuple = ("zipf", 1.821),  # (Chang and Benner 2004)
        min_length: int = 1,
        max_length: int | None = None,
        **kwargs,
    ) -> None:
        """Constructor for the IndelModel class.

        Args:
            insertion_rate: Insertion rate w.r.t. units of evolutionary distance.
            deletion_rate: Deletion rate w.r.t. units of evolutionary distance.
            length_distr: Distribution of the length of the indels, see documentation for available
                options. The default is a Zipf distribution as observed for empirical data in [1].
            min_length: The minimal length of an indel, default is 1.
            max_length: The maximal length of an indel. The default is None in which case the
                distribution is not truncated (on the upper side).

        Raises:
            ValueError: If insertion_rate or deletion_rate is negative.

        References:
            1. M. S. Chang and S. A. Benner. Empirical Analysis of Protein Insertions and
               Deletions Determining Parameters for the Correct Placement of Gaps in Protein
               Sequence Alignments. In: Journal of Molecular Biology, 341(2):617-631, August 2004.
               doi:10.1016/j.jmb.2004.05.045.
        """
        if insertion_rate < 0.0 or deletion_rate < 0.0:
            raise ValueError("insertion and deletion rates must be non-negative")

        self._ins_rate = insertion_rate  # insertion rate per site
        self._del_rate = deletion_rate  # deletion rate per site

        if isinstance(length_distr, Sampler):
            self.sampler = length_distr
        else:
            self.sampler = Sampler(
                length_distr, minimum=min_length, maximum=max_length, discrete=True
            )

    def get_rates(self, seq_length: int) -> tuple[float, float]:
        """Return the current insertion and deletion rate.

        Computes the current insertion and deletion rate according to the model and the length of
        the sequence.

        Args:
            seq_length: The length of the sequence.

        Returns:
            The total insertion rate for the given sequence length, see 1.
            The total deletion rate for the given sequence length, see 1.

        References:
            1. R. A. Cartwright. DNA assembly with gaps (Dawg): Simulating sequence evolution.
               In: Bioinformatics, 21(Suppl 3):iii31-iii38, November 2005.
               doi:10.1093/bioinformatics/bti1200.
        """
        # expected value may be infinite for zipf distribution, in which case the minimal length is
        # used, see Sampler class for details
        return (
            (seq_length + 1) * self._ins_rate,
            (seq_length + self.sampler._exp_val - 1) * self._del_rate,
        )

    def draw_length(self) -> int:
        """Draw the length for an indel.

        Returns:
            The drawn length for the indel.
        """
        return self.sampler()
