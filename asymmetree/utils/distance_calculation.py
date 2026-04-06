"""Calculation of sequence distances for various models.

The distance between two sequences is defined as the expected number of substitutions per site that
have occurred since the two sequences diverged from their common ancestor. This module provides
functions to calculate distances under various substitution models, including the Jukes-Cantor 1969
model, the Kimura 1980 model, and maximum likelihood distances under general substitution models.

References:
    1. Z. Yang. Computational molecular evolution. Oxford series in ecology and evolution.
    Oxford University Press, 2006. ISBN 978-0-19-856699-1 978-0-19-856702-8.
"""

from __future__ import annotations

import numpy as np
import scipy.optimize

from asymmetree.seqevolve import SubstModel


MIN_LENGTH = 1e-6
MAX_LENGTH = 20.0


def maximum_likelihood_distance(
    seq1: str,
    seq2: str,
    subst_model: SubstModel = None,
    model_type: str | None = None,
    model_name: str | None = None,
) -> float:
    """Calculate the maximum likelihood distance between two aligned sequences.

    Args:
        seq1: First aligned sequence.
        seq2: Second aligned sequence.
        subst_model: Substitution model to use for distance calculation. If None, model_type and
            model_name must be provided to construct a SubstModel instance.
        model_type: Type of substitution model (e.g., "n" for nucleotide, "a" for amino acid).
            Required if subst_model is None.
        model_name: Name of substitution model (e.g., "JC69", "K80"). Required if subst_model is
            None.

    Returns:
        Estimated distance between the two sequences under the specified substitution model.

    Raises:
        ValueError: If subst_model is None and either model_type or model_name is not provided.
    """
    if subst_model is None:
        if model_type is None or model_name is None:
            raise ValueError("no substitution model specified")

        subst_model = SubstModel(model_type, model_name)

    seqs = _to_indices(seq1, seq2, subst_model)

    x0 = _initial_guess(seqs, subst_model.model_type)

    # distance 0.0 if sequences are equal
    if x0 == 0.0:
        return x0

    opt_result = scipy.optimize.minimize(
        _likelihood,
        x0,
        args=(seqs, subst_model),
        bounds=((MIN_LENGTH, None),),
    )

    if opt_result.success:
        return opt_result.x[0]
    else:
        return float("nan")


def _initial_guess(seqs: np.ndarray, model_type: str) -> float:
    """Calculate an initial guess for the distance optimization.

    The initial guess is based on the p-distance and the Jukes-Cantor 1969 transformation, which
    provides a simple correction for multiple substitutions. If the p-distance is too high (i.e.,
    above the threshold for the Jukes-Cantor model), the initial guess is set to half of the maximum
    allowed distance to avoid starting the optimization in a region where the likelihood is
    undefined.

    Args:
        seqs: A 2xN array of integer indices representing the aligned sequences.
        model_type: The type of substitution model (e.g., "n" for nucleotide, "a" for amino acid) to
            determine the appropriate Jukes-Cantor transformation.

    Returns:
        An initial guess for the distance between the two sequences.
    """
    p = np.sum(seqs[0] != seqs[1]) / seqs.shape[1]

    x0 = _JC69_transform(p, amino_acid=(model_type == "a"))

    x0 = x0 if x0 != float("nan") else MAX_LENGTH / 2

    return x0


def _likelihood(x, seqs: np.ndarray, subst_model: SubstModel) -> float:
    """Calculate the negative log-likelihood of the observed sequence differences given a distance.

    Args:
        x: The distance (expected number of substitutions per site) to evaluate the likelihood at.
        seqs: A 2xN array of integer indices representing the aligned sequences.
        subst_model: The substitution model to use for calculating the transition probabilities.

    Returns:
        The negative log-likelihood of the observed sequence differences given the distance x.
    """
    P = subst_model.transition_prob_matrix(x)

    # - log-likelihood
    return -np.sum(np.log(P[seqs[0], seqs[1]]))


def _to_indices(seq1: str, seq2: str, subst_model: SubstModel) -> np.ndarray:
    """Convert sequences to integer indices based on the substitution model's alphabet.

    Args:
        seq1: First aligned sequence.
        seq2: Second aligned sequence.
        subst_model: The substitution model whose alphabet is used for indexing.

    Returns:
        A 2xN array of integer indices corresponding to the aligned sequences, excluding positions
            with gaps in either sequence.

    Raises:
        ValueError: If the input sequences have different lengths.
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"unequal sequence lengths: {len(seq1)} and {len(seq2)}")

    seq1_indeces, seq2_indeces = [], []

    for i in range(len(seq1)):
        if seq1[i] != "-" and seq2[i] != "-":
            seq1_indeces.append(subst_model.alphabet_dict[seq1[i]])
            seq2_indeces.append(subst_model.alphabet_dict[seq2[i]])

    return np.asarray([seq1_indeces, seq2_indeces], dtype=int)


def p_distance(seq1: str, seq2: str, exclude_gaps: bool = True) -> tuple[float, int]:
    """Calculate the p-distance of two aligned sequences (normalized Hamming distance).

    The p-distance is the proportion of sites at which the two sequences differ. It equals the
    Hamming distance divided by the number of valid columns (i.e., columns without gaps if
    exclude_gaps is True).

    Args:
        seq1: First aligned sequence.
        seq2: Second aligned sequence.
        exclude_gaps: If True, columns with a gap in one sequence are ignored. Columns with gaps in
            both sequences are always ignored.

    Returns:
        The p-distance.
        The number of valid columns used in the calculation.

    Raises:
        ValueError: If the input sequences have different lengths.
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"unequal sequence lengths: {len(seq1)} and {len(seq2)}")

    diffs, valid_columns = 0, 0

    for i in range(len(seq1)):
        if seq1[i] != "-" and seq2[i] != "-":
            valid_columns += 1
            if seq1[i] != seq2[i]:
                diffs += 1

        elif not exclude_gaps and (seq1[i] != "-" or seq2[i] != "-"):
            valid_columns += 1
            diffs += 1

    p = diffs / valid_columns if valid_columns > 0 else float("inf")

    return p, valid_columns


def JC69_distance(
    seq1: str,
    seq2: str,
    exclude_gaps: bool = True,
    amino_acid: bool = False,
    variance: bool = False,
) -> float | tuple[float, float]:
    """Jukes Cantor 1969 distance.

    The Jukes-Cantor 1969 model assumes equal base frequencies and equal substitution rates among
    all nucleotides.

    Args:
        seq1: First aligned sequence.
        seq2: Second aligned sequence.
        exclude_gaps: If True, columns with a gap in one sequence are ignored. Columns with gaps in
            both sequences are always ignored.
        amino_acid: If True, the distance is calculated under the Jukes-Cantor model for amino
            acids, which assumes equal frequencies and substitution rates among all amino acids.
        variance: If True, also calculate the variance of the distance estimate.

    Returns:
        If `variance` is False, returns the Jukes-Cantor distance. If `variance` is True, returns a
            tuple containing the distance and its variance.
    """
    p, valid_columns = p_distance(seq1, seq2, exclude_gaps=exclude_gaps)

    d = _JC69_transform(p, amino_acid=amino_acid)

    if not variance:
        return d
    else:
        var_d = _JC69_distance_var(p, valid_columns, amino_acid=amino_acid)
        return d, var_d


def _JC69_transform(p: float, amino_acid: bool = False) -> float:
    """Jukes Cantor 1969 transformation of p-distance.

    Args:
        p: The p-distance.
        amino_acid: If True, use the amino acid model.

    Returns:
        The Jukes-Cantor distance.
    """
    a = 3 if not amino_acid else 19  # numerator
    b = 4 if not amino_acid else 20  # denominator

    if p >= a / b:
        d = float("nan")

    else:
        d = -(a / b) * np.log(1 - (b / a) * p)

    return d


def _JC69_distance_var(p: float, n: int, amino_acid: bool = False) -> float:
    """Jukes Cantor 1969 distance variance.

    Args:
        p: The p-distance.
        n: The number of valid columns.
        amino_acid: If True, use the amino acid model.

    Returns:
        The variance of the Jukes-Cantor distance.
    """
    a = 3 if not amino_acid else 19  # numerator
    b = 4 if not amino_acid else 20  # denominator

    if p >= a / b:
        var = float("nan")

    else:
        var = p * (1 - p) / (n * (1 - (b / a) * p) ** 2)

    return var


def K80_distance(
    seq1: str, seq2: str, variance: bool = False
) -> tuple[float, float] | tuple[float, float, float]:
    """Kimura 1980 distance and transition/transversion ratio.

    The Kimura 1980 model distinguishes between transitions (substitutions between purines or
    pyrimidines) and transversions (substitutions between a purine and a pyrimidine). It assumes
    equal base frequencies and equal substitution rates among all nucleotides, but allows
    for different rates between transitions and transversions.

    Args:
        seq1: First aligned sequence.
        seq2: Second aligned sequence.
        variance: If True, also calculate the variance of the distance estimate.

    Returns:
        If `variance` is False, returns a tuple containing the Kimura distance and the transition /
            transversion ratio. If `variance` is True, returns a tuple containing the distance, the
            transition / transversion ratio, and the variance of the distance.
    """
    S, V, valid_columns = _IV_proportions(seq1, seq2)

    d, kappa = _K80_transform(S, V)

    if not variance:
        return d, kappa
    else:
        var_d = _K80_distance_var(S, V, valid_columns)
        return d, kappa, var_d


def _IV_proportions(seq1: str, seq2: str) -> tuple[float, float, int]:
    """Computes the proportions of transitions and transversions.

    Args:
        seq1: First aligned sequence.
        seq2: Second aligned sequence.

    Returns:
        The proportion of transitions.
        The proportion of transversions.
        The number of valid columns.

    Raises:
        ValueError: If the input sequences have different lengths.
    """
    purines = {"A", "a", "G", "g", "R", "r"}
    pyrimidines = {"C", "c", "T", "t", "U", "u", "Y", "y"}

    if len(seq1) != len(seq2):
        raise ValueError(f"unequal sequence lengths: {len(seq1)} and {len(seq2)}")

    transitions, transversions, valid_columns = 0, 0, 0

    for i in range(len(seq1)):
        if seq1[i] != "-" and seq2[i] != "-":
            valid_columns += 1

            if seq1[i] == seq2[i]:
                continue

            elif (seq1[i] in purines and seq2[i] in purines) or (
                seq1[i] in pyrimidines and seq2[i] in pyrimidines
            ):
                transitions += 1

            elif (seq1[i] in purines and seq2[i] in pyrimidines) or (
                seq1[i] in pyrimidines and seq2[i] in purines
            ):
                transversions += 1

    S = transitions / valid_columns if valid_columns > 0 else float("nan")
    V = transversions / valid_columns if valid_columns > 0 else float("nan")

    return S, V, valid_columns


def _K80_transform(S: float, V: float) -> tuple[float, float]:
    """Kimura 1980 distance and transition / transversion ratio.

    Args:
        S: Proportion of transitions.
        V: Proportion of transversions.

    Returns:
        The Kimura distance.
        The transition / transversion ratio.
    """
    a1 = 1.0 - 2.0 * S - V
    a2 = 1.0 - 2.0 * V

    if a1 <= 0.0 or a2 <= 0.0:
        return float("nan"), float("nan")
    else:
        d = -0.5 * np.log(a1) - 0.25 * np.log(a2)
        kappa = 2.0 * np.log(a1) / np.log(a2) - 1.0 if a2 < 1.0 else float("nan")

        return d, kappa


def _K80_distance_var(S: float, V: float, n: int) -> float:
    """Kimura 1980 distance variance.

    Args:
        S: Proportion of transitions.
        V: Proportion of transversions.
        n: Number of valid columns.

    Returns:
        The variance of the Kimura 1980 distance.
    """
    a1 = 1.0 - 2.0 * S - V
    a2 = 1.0 - 2.0 * V

    if a1 <= 0.0 or a2 <= 0.0:
        return float("nan")
    else:
        a = 1.0 / a1
        b = 0.5 * (1.0 / a1 + 1.0 / a2)

        return (a**2 * S + b**2 * V - (a * S + b * V) ** 2) / n
