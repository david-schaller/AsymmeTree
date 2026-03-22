"""Simulation of noisy distance matrices.

Distances derived from (real-life) gene or protein sequences are burdened with noise. Such data can
either be modeled by simulating sequences, or by disturbing the distances specified by a given tree
directly. This module implements functions for the latter alternative.
"""

from __future__ import annotations

import random
import numpy as np

from tralda.datastructures.Tree import Tree, TreeNode
from asymmetree.utils.phylogenetic_trees import distance_matrix


# --------------------------------------------------------------------------------------------------
#                                  RANDOM PERTURBATION NOISE
# --------------------------------------------------------------------------------------------------


def noisy_matrix(orig_matrix: np.ndarray, sd: float, metric_repair: str = "reject") -> np.ndarray:
    """Disturb a distance matrix (which must be a metric) with random noise.

    Args:
        orig_matrix: The 2-dimensional distance matrix.
        sd: Disturbance parameter. Serves as standard deviation of a normal distrubition with mean
            1.0.
        metric_repair: Strategy to ensure that the resulting matrix is still a metric and, in
            particular, satisfies the triangle inequality. The default is 'reject', see [1]. Other
            available options are 'general' and 'DOMR', see [2].

    Returns:
        The disturbed distance matrix.

    Raises:
        ValueError: If the argument 'metric_repair' is not one of 'reject', 'general', or 'DOMR'.

    References:
        .. [1] P. F. Stadler, M. Geiß, D. Schaller, A. López Sánchez, M. González Laffitte, D.
           Valdivia, M. Hellmuth, M. Hernández Rosales. From pairs of most similar sequences to
           phylogenetic best matches. In: Alg. Mol. Biol., 2020, 15, 5.
           doi: 10.1186/s13015-020-00165-2.
        .. [2] A. C. Gilbert, L. Jain. If it ain't broke, don't fix it: Sparse metric repair.
           In: 2017 55th Annual Allerton Conference on Communication, Control, and Computing
           (Allerton), pages 612-619, Monticello, IL, USA, October 2017. IEEE.
           ISBN 978-1-5386-3266-6. doi: 10.1109/ALLERTON.2017.8262793.
    """
    if metric_repair == "general":
        return _noise_general_metric_repair(orig_matrix, sd)
    elif metric_repair == "DOMR":
        return _noise_metric_repair_DOMR(orig_matrix, sd)
    elif metric_repair == "reject":
        return _noise_reject_method(orig_matrix, sd)

    raise ValueError(f"illegal argument '{metric_repair}'")


def _noise_reject_method(orig_matrix: np.ndarray, sd: float) -> np.ndarray:
    """Return a matrix D' with noise by accept/reject algorithm.

    Args:
        orig_matrix: The original distance matrix.
        sd: Standard deviation of the normal distribution used for noise.

    Returns:
        The disturbed distance matrix.
    """
    D = np.array(orig_matrix, copy=True)
    N = D.shape[0]

    success_count = 0
    stop_at = N * (N - 1) / 2

    while success_count < stop_at:
        i, j = np.random.randint(N, size=2)
        while i == j:
            j = np.random.randint(N)
        old_distance = D[i, j].copy()
        new_distance = 0.0
        while new_distance <= 0.0:
            new_distance = D[i, j] * np.random.normal(loc=1.0, scale=sd)
        D[i, j] = new_distance
        D[j, i] = new_distance

        if (
            D[i, j] > np.min(D[i, :] + D[:, j])
            or np.any(D[:, i] > D[:, j] + D[j, i])
            or np.any(D[:, j] > D[:, i] + D[i, j])
        ):
            D[i, j] = old_distance
            D[j, i] = old_distance
        else:
            success_count += 1

    return D


def _noise_metric_repair_DOMR(orig_matrix: np.ndarray, sd: float) -> np.ndarray:
    """Return a matrix D' with noise by metric repair (DOMR) algorithm.

    Args:
        orig_matrix: The original distance matrix.
        sd: Standard deviation of the normal distribution used for noise.

    Returns:
        The disturbed distance matrix.
    """
    D = np.array(orig_matrix, copy=True)
    N = D.shape[0]

    for i in range(N - 1):  # noise introduction
        for j in range(i + 1, N):
            new_distance = 0.0
            while new_distance <= 0.0:
                new_distance = D[i, j] * np.random.normal(loc=1.0, scale=sd)
            D[i, j] = new_distance
            D[j, i] = new_distance

    for k in range(N):  # metric repair: decrease
        for i in range(N):  # only metric repair (DOMR)
            for j in range(i):  # with Floyd-Warshall
                if D[i, j] >= D[i, k] + D[k, j]:
                    D[i, j] = D[i, k] + D[k, j]
                    D[j, i] = D[i, k] + D[k, j]

    return D


def _noise_general_metric_repair(orig_matrix: np.ndarray, sd: float) -> np.ndarray:
    """Return a matrix D' with noise by general metric repair algorithm.

    Args:
        orig_matrix: The original distance matrix.
        sd: Standard deviation of the normal distribution used for noise.

    Returns:
        The disturbed distance matrix.
    """
    D = np.array(orig_matrix, copy=True)
    N = D.shape[0]

    for i in range(N - 1):  # noise introduction
        for j in range(i + 1, N):
            new_distance = 0.0
            while new_distance <= 0.0:
                new_distance = D[i, j] * np.random.normal(loc=1.0, scale=sd)
            D[i, j] = new_distance
            D[j, i] = new_distance

    matrix_l, matrix_r = np.zeros_like(D), np.zeros_like(D)
    for k in range(N):  # metric repair: general metric repair
        for i in range(N):
            for j in range(i):
                if D[i, j] >= D[i, k] + D[k, j]:
                    matrix_l[i, j] += 1
                    matrix_r[i, k] += 1
                    matrix_r[j, k] += 1

    for k in range(N):  # metric repair: general metric repair
        for i in range(N):
            for j in range(i):
                if D[i, j] >= D[i, k] + D[k, j]:
                    if matrix_l[i, j] > max(matrix_r[i, k], matrix_r[j, k]):
                        D[i, j] = D[i, k] + D[k, j]
                        D[j, i] = D[i, k] + D[k, j]
                    elif matrix_r[i, k] > matrix_r[j, k]:
                        D[i, k] = D[i, j] - D[j, k]
                        D[k, i] = D[i, j] - D[j, k]
                    else:
                        D[j, k] = D[i, j] - D[i, k]
                        D[k, j] = D[i, j] - D[i, k]

    return D


# --------------------------------------------------------------------------------------------------
#                                    WRONG TOPOLOGY NOISE
# --------------------------------------------------------------------------------------------------


def convex_linear_comb(
    D1: np.ndarray,
    D2: np.ndarray,
    alpha: float = 0.05,
    first_only: bool = False,
) -> np.ndarray:
    """Convex linear combination of distance matrices.

    Returns the convex linear combinations of two distance matrices
        D1' = (1-alpha) * D1 + alpha * D2
        D2' = (1-alpha) * D2 + alpha * D1

    The motivation for this function is to simulate non-random noise from a wrong topology, which
    can be used to test the robustness of tree reconstruction methods.

    Args:
        D1: The first distance matrix.
        D2: The second distance matrix.
        alpha: Disturbance parameter, the default is 0.05.
        first_only: If True, only return the disturbed first matrix.

    Returns:
        The disturbed distance matrix or matrices.

    References:
        .. [1] P. F. Stadler, M. Geiß, D. Schaller, A. López Sánchez, M. González Laffitte, D.
           Valdivia, M. Hellmuth, M. Hernández Rosales. From pairs of most similar sequences to
           phylogenetic best matches. In: Alg. Mol. Biol., 2020, 15, 5.
           doi: 10.1186/s13015-020-00165-2.
    """
    if not first_only:
        if D1.shape == D2.shape:
            return ((1 - alpha) * D1 + alpha * D2, (1 - alpha) * D2 + alpha * D1)

        if D1.shape > D2.shape:
            D1, D2 = D2, D1  # now D1 is smaller
            swap = True
        else:
            swap = False

        N1, N2 = D1.shape[0], D2.shape[0]
        D1_alpha = (1 - alpha) * D1 + alpha * D2[:N1, :N1]
        D2_alpha = (1 - alpha) * D2

        for i in range(N2 // N1):
            for j in range(N2 // N1):
                D2_alpha[i * N1 : i * N1 + N1, j * N1 : j * N1 + N1] += alpha * D1
                if i == (N2 // N1) - 1 and i == j:
                    D2_alpha[i * N1 + N1 : N2, j * N1 + N1 : N2] += (
                        alpha * D1[: (N2 % N1), : (N2 % N1)]
                    )
                if i == (N2 // N1) - 1:
                    D2_alpha[i * N1 + N1 : N2, j * N1 : j * N1 + N1] += alpha * D1[: (N2 % N1), :]
                if j == (N2 // N1) - 1:
                    D2_alpha[i * N1 : i * N1 + N1, j * N1 + N1 : N2] += alpha * D1[:, : (N2 % N1)]

        if swap:
            return D2_alpha, D1_alpha
        else:
            return D1_alpha, D2_alpha

    else:
        if D1.shape == D2.shape:
            return (1 - alpha) * D1 + alpha * D2

        N1, N2 = D1.shape[0], D2.shape[0]
        if N1 < N2:
            return (1 - alpha) * D1 + alpha * D2[:N1, :N1]
        else:
            D1_alpha = (1 - alpha) * D1
            for i in range(N1 // N2):
                for j in range(N1 // N2):
                    D1_alpha[i * N2 : i * N2 + N2, j * N2 : j * N2 + N2] += alpha * D2
                    if i == (N1 // N2) - 1 and i == j:
                        D1_alpha[i * N2 + N2 : N1, j * N2 + N2 : N1] += (
                            alpha * D2[: (N1 % N2), : (N1 % N2)]
                        )
                    if i == (N1 // N2) - 1:
                        D1_alpha[i * N2 + N2 : N1, j * N2 : j * N2 + N2] += (
                            alpha * D2[: (N1 % N2), :]
                        )
                    if j == (N1 // N2) - 1:
                        D1_alpha[i * N2 : i * N2 + N2, j * N2 + N2 : N1] += (
                            alpha * D2[:, : (N1 % N2)]
                        )
            return D1_alpha


def wrong_topology_matrix(PGT: Tree) -> np.ndarray:
    """Return a wrong topology matrix by rearranging the edges of a binary tree.

    Args:
        PGT: The input binary tree, which is typically a pruned gene tree (PGT).

    Returns:
        The distance matrix of the wrong topology.
    """
    # do not include root
    distances = [v.dist for v in PGT.preorder()][1:]

    # pruned gene tree (PGT) is binary and not planted, hence |E| should be even
    if len(distances) % 2 != 0:
        raise ValueError("The input tree is not binary!")

    random.shuffle(distances)

    random_tree = Tree(TreeNode(label=0, dist=0.0))
    id_counter = 1
    current_leaves = [random_tree.root]

    while distances:
        v = current_leaves.pop(random.randint(0, len(current_leaves) - 1))
        dist1, dist2 = distances.pop(), distances.pop()
        new_child1 = TreeNode(label=id_counter, dist=dist1)
        new_child2 = TreeNode(label=id_counter + 1, dist=dist2)
        v.add_child(new_child1)
        v.add_child(new_child2)
        current_leaves.extend(v.children)
        id_counter += 1

    # implicit random bijection to original tree by random leaf order
    random_leaves = [leaf for leaf in random_tree.leaves()]
    random.shuffle(random_leaves)

    _, D = distance_matrix(random_tree, leaf_order=random_leaves)

    return D


# --------------------------------------------------------------------------------------------------
#                                         METRIC CHECK
# --------------------------------------------------------------------------------------------------


def check_metric(matrix: np.ndarray, verbose: bool = False) -> bool:
    """Check whether a given matrix is a metric.

    Args:
        matrix: The distance matrix to check.
        verbose: If True, print the reason for failure if the matrix is not a metric.

    Returns:
        True if the matrix is a metric, False otherwise.
    """
    N = matrix.shape[0]

    for i in range(N):
        if matrix[i, i] != 0.0:
            if verbose:
                print("Not all diagonal elements zero!")
            return False

    for i in range(N - 1):
        for j in range(i + 1, N):
            if matrix[i, j] != matrix[j, i]:
                if verbose:
                    print("Not symmetrical for", i, j)
                return False

            if abs(matrix[i, j] - np.min(matrix[i, :] + matrix[:, j])) > 1e-8:
                if verbose:
                    print(
                        "Violation of triangle inequality!",
                        i,
                        j,
                        np.min(matrix[i, :] + matrix[:, j]),
                        matrix[i, j],
                    )
                return False

    return True
