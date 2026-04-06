"""Heterogeneity models for sequence evolution.

Not all sites in a sequence evolve at the same rate. This module provides a class for modeling rate
heterogeneity across sites, which is a common feature of sequence evolution. The class supports
rate heterogeneity based on a Gamma distribution ('+Gamma'-model) as well as invariant sites
('+I'-model), see [1].

References:
    1. Z. Yang. Computational molecular evolution. Oxford series in ecology and evolution.
       Oxford University Press, 2006. ISBN 978-0-19-856699-1 978-0-19-856702-8.
"""

from __future__ import annotations

import numpy as np

from asymmetree.seqevolve.evolving_sequence import EvoSeq
from asymmetree.seqevolve.evolving_sequence import State


class HetModel:
    """Heterogeneity model.

    Supports rate heterogeneity based on a Gamma distribution ('+Gamma'-model) as well as invariant
    sites ('+I'-model), see [1].

    References:
    1. Z. Yang. Computational molecular evolution. Oxford series in ecology and evolution.
       Oxford University Press, 2006. ISBN 978-0-19-856699-1 978-0-19-856702-8.
    """

    def __init__(
        self,
        alpha: float,
        classes: int = 5,
        sitewise: bool = False,
        invariant: float = 0.0,
    ) -> None:
        """Constructor for the HetModel class.

        Args:
            alpha: Parameter of the Gamma distribution (with mean 1) from which the rate factors are
                drawn.
            classes: Number of classes such that sites in the same class share the same rate factor.
            sitewise: If True, ignore the 'classes' attribute and treat each site as its own class.
            invariant: Expected proportion of invariant sites, i.e., sites where no mutations
                happen.

        Raises:
            ValueError: If alpha is not a float > 0.0, if classes is not an int > 0, or if
                invariant is not a float in [0.0, 1.0].
        """
        if not isinstance(alpha, float) or alpha <= 0.0:
            raise ValueError("heterogeneity parameter alpha must be a float > 0.0")

        self._alpha = alpha

        if not isinstance(classes, int) or classes <= 0:
            raise ValueError("number of classes must be an int >0")

        self.classes = classes
        self.sitewise = sitewise

        if not isinstance(invariant, float) or invariant < 0.0 or invariant > 1.0:
            raise ValueError("proportion of invariant sites must be in [0.0, 1.0]")

        self._invariant = invariant

        if not self.sitewise:
            self._initialize_classes()

    def assign(self, sequence: EvoSeq, exclude_inherited: bool = True) -> None:
        """Assign rate classes and rate factors to the sites of a sequence.

        Args:
            sequence: The sequence to which the rate classes and rate factors should be assigned.
            exclude_inherited: If True, do not assign new classes and rates to inherited sites,
                i.e., only if the sequence corresponds to the root of a tree along which is
                simulated or if new sites are added as a result of an insertion event.
        """
        n = len(sequence)
        if exclude_inherited:
            n -= sequence.count_status(State.INHERITED)

        rate_classes, rate_factors = self._draw(n)

        pos = 0
        for site in sequence:
            if not exclude_inherited or site.status != State.INHERITED:
                site.rate_class = (
                    rate_classes[pos] if rate_classes[pos] is not None else site.site_id
                )
                site.rate_factor = rate_factors[pos]
                pos += 1

    def _initialize_classes(self) -> None:
        """Initialize the rate factors for the classes."""
        if self.classes == 1:
            self._class_rates = [1.0]
        else:
            self._class_rates = np.random.gamma(
                self._alpha, scale=1 / self._alpha, size=self.classes
            )

    def _draw(self, n) -> tuple[list[int | str | None], np.ndarray]:
        """Draw rate classes and rate factors for n sites.

        Args:
            n: The number of sites for which to draw the rate classes and rate factors.

        Returns:
            A list of rate classes.
            An array of rate factors for the `n` sites.
        """
        # mode 1: sitewise heterogeneity
        if self.sitewise:
            rate_classes = [None for _ in range(n)]
            rate_factors = np.random.gamma(self._alpha, scale=1 / self._alpha, size=n)

        # mode 2: one or multiple  classes
        else:
            rate_classes = np.random.randint(self.classes, size=n).tolist()
            rate_factors = np.asarray([self._class_rates[c] for c in rate_classes])

        if self._invariant:
            drawn_variable = np.random.random(n) > self._invariant

            for i in range(n):
                if not drawn_variable[i]:
                    rate_classes[i] = "invariant"

            rate_factors *= drawn_variable

        return rate_classes, rate_factors
