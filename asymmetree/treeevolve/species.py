"""Simulation of dated species trees."""

from __future__ import annotations

import random
import numpy as np

from tralda.datastructures import Tree
from tralda.datastructures import TreeNode


from asymmetree.utils.phylogenetic_trees import delete_losses_and_contract
from asymmetree.utils.phylogenetic_trees import distance_from_timing
from asymmetree.utils.phylogenetic_trees import remove_planted_root


def species_tree_n(
    n: int,
    model: str = "yule",
    innovation: bool = False,
    planted: bool = True,
    remove_extinct: bool = False,
    **kwargs,
) -> Tree:
    """Simulates a species tree S with n leaves.

    The tree is sampled under the Yule, BDP, or EBDP model conditioned on the number n of surviving
    species [1, 2, 3].

    Args:
        n: Number of leaves in the resulting tree that correspond to extant species.
        model: Simulation model to be applied, the default is 'yule', see [1]. Other available
            models are birth-death process ('BDP', [2]), and episodic birth-death process ('EBDP',
            [3]).
        innovation: If True, use the innovation model, see [4], to sample a lineage for the next
            speciation event. Only available for the Yule model. The default is False, in which case
            the lineage is chosen uniformly at random among the currently existing lineages.
        planted: Add a planted root that has the canonical root as its single neighbor.
        remove_extinct: Remove all branches that lead to extinctions, only relevant for some models.
        birth_rate: The birth rate for models such as 'yule' and 'BDP'.
        death_rate: The death rate for models such as 'BDP'.
        episodes: The episodes for the model 'EBDP'.
        contraction_probability: Probability that an inner edge is contracted. The default is 0.0,
            in which case the tree is binary. Only one of this parameter and
            'contraction_proportion' may be non-zero.
        contraction_proportion: The proportion of inner edges to be contracted. The default is 0.0,
            in which case the tree is binary. Only one of this parameter and
            'contraction_probability' may be non-zero.
        contraction_bias: Specifies whether shorter edges, i.e., with a smaller difference t of the
            time stamps, have a higher probability to be contracted. Only relevant if
            'contraction_proportion' > 0.0. The default is False, in which case all edges have the
            same probability to be contracted. The options 'inverse' and 'exponential' mean that an
            edge is sampled weighted by 1/t or e^(-t), respectively.
        bias_strength: Intensity factor for preferring shorter edges to be contracted. The default
            is 1.0.

    Raises:
        ValueError: If unknown or invalid parameter values are passed.

    Returns:
        The simulated species tree.

    References:
    ----------
        .. [1] G. U. Yule. A mathematical theory of evolution, based on the conclusions of Dr. J.
           C. Willis, F. R. S. In: Phil. Trans. R. Soc. Lond. B, 1924, 213, 21-87.
           doi:10.1098/rstb.1925.0002.
        .. [2] D. G. Kendall. On the Generalized "Birth-and-Death" Process. In: Ann. Math. Statist.
           1948, 19, 1-15. doi:10.1214/aoms/1177730285.
        .. [3] T. Stadler. Simulating trees with a fixed number of extant species. In: Syst. Biol.
           2011, 60, 676-684. doi:10.1093/sysbio/syr029.
        .. [4] S. Keller-Schmidt, K. Klemm. A model of macroevolution as a branching process based
           on innovations. In: Adv. Complex Syst.,2012, 15, 1250043. doi:10.1142/S0219525912500439.
    """
    # parameter checking
    if not isinstance(n, int) or n < 0:
        raise ValueError("n must be an int >=0")
    elif n == 0:
        return Tree(None)

    if not isinstance(model, str):
        raise ValueError("model must be of type 'str'")

    # choice of the main simulation algorithm
    if model.lower() == "yule":
        tree = _yule_n(n, kwargs.get("birth_rate"), innovation)
    elif model.upper() == "BDP":
        tree = _BDP_n(n, **kwargs)
    elif model.upper() == "EBDP":
        tree = _EBDP_n(n, **kwargs)
    else:
        raise ValueError(f"model '{model}' is not available")

    # remove extinct branches for models that include losses
    if remove_extinct and model.upper() in ("BDP", "EBDP"):
        delete_losses_and_contract(tree, inplace=True)

    # remove planted edge for models that are planted by construction
    if not planted and model.upper() in ("YULE", "BDP", "EBDP"):
        remove_planted_root(tree, inplace=True)

    # make tree non_binary by random contraction of edges
    nonbinary(tree, inplace=True, **kwargs)

    # assign the distance attribute to all vertices
    distance_from_timing(tree)

    return tree


def species_tree_age(age: float, model: str = "yule", innovation: bool = False, **kwargs) -> Tree:
    """Simulates a (planted) species tree S of the specified age.

    Args:
        age: Simulation time, i.e., the time span from the root of the tree to the leaves that
            correspond to extant species.
        model: Simulation model to be applied, the default is 'yule', see [1]. Other available
            models are birth-death process ('BDP', [2]), and episodic birth-death process ('EBDP',
            [3]).
        innovation: If True, use the innovation model, see [4], to sample a lineage for the next
            speciation event. The default is False, in which case the lineage is chosen uniformly at
            random among the currently existing lineages.
        birth_rate: The birth rate for models such as 'yule' and 'BDP'.
        death_rate: The death rate for models such as 'BDP'.
        episodes: The episodes for the model 'EBDP'.
        contraction_probability: Probability that an inner edge is contracted. The default is 0.0,
            in which case the tree is binary. Only one of this parameter and
            'contraction_proportion' may be non-zero.
        contraction_proportion: The proportion of inner edges to be contracted. The default is 0.0,
            in which case the tree is binary. Only one of this parameter and
            'contraction_probability' may be non-zero.
        contraction_bias: Specifies whether shorter edges, i.e., with a smaller difference t of the
            time stamps, have a higher probability to be contracted. Only relevant if
            'contraction_proportion' > 0.0. The default is False, in which case all edges have the
            same probability to be contracted. The options 'inverse' and 'exponential' mean that an
            edge is sampled weighted by 1/(a * t) or e^(-a * t), respectively, where a is a
            user-defined factor.
        bias_strength: Intensity factor for preferring shorter edges to be contracted. The default
            is 1.0.

    Raises:
        ValueError: If unknown or invalid parameter values are passed.

    Returns:
        The simulated species tree.

    References:
    ----------
        .. [1] G. U. Yule. A mathematical theory of evolution, based on the conclusions of Dr. J.
           C. Willis, F. R. S. In: Phil. Trans. R. Soc. Lond. B, 1924, 213, 21-87.
           doi:10.1098/rstb.1925.0002.
        .. [2] D. G. Kendall. On the Generalized "Birth-and-Death" Process. In: Ann. Math. Statist.
           1948, 19, 1-15. doi:10.1214/aoms/1177730285.
        .. [3] T. Stadler. Simulating trees with a fixed number of extant species. In: Syst. Biol.
           2011, 60, 676-684. doi:10.1093/sysbio/syr029.
        .. [4] S. Keller-Schmidt, K. Klemm. A model of macroevolution as a branching process based
           on innovations. In: Adv. Complex Syst.,2012, 15, 1250043. doi:10.1142/S0219525912500439.
    """
    # parameter checking
    if not isinstance(age, (float, int)) or age <= 0.0:
        raise ValueError("age must be a number >0")
    elif isinstance(age, int):
        age = float(age)

    if not isinstance(model, str):
        raise ValueError("model must be of type 'str'")

    # main simulation algorithm
    if model.lower() == "yule":
        tree = _yule_age(age, kwargs.get("birth_rate"), innovation)
    elif model.upper() == "BDP":
        tree = _BDP_age(age, innovation, **kwargs)
    elif model.upper() == "EBDP":
        tree = _EBDP_age(age, innovation, **kwargs)
    else:
        raise ValueError(f"model '{model}' is not available")

    # make tree non_binary by random contraction of edges
    nonbinary(tree, inplace=True, **kwargs)

    # assign the distance attribute to all vertices
    distance_from_timing(tree)

    return tree


def species_tree_n_age(
    n: int,
    age: float,
    model: str = "yule",
    innovation: bool = False,
    birth_rate: float = 1.0,
    death_rate: float = 0.0,
    **kwargs,
) -> Tree:
    """Simulate a (planted) species tree S with n leaves and of the specified age.

    The tree is sampled under the Yule or BDP model conditioned on the number n of surviving species
    and with a given age [1]. The resulting tree does not contain loss leaves even if the specified
    death rate is > 0 as a consequence of the tree sampling method.

    Args:
        n: Number of leaves in the resulting tree that correspond to extant species.
        age: Simulation time, i.e., the time span from the root of the tree to the leaves that
            correspond to extant species.
        model: Simulation model to be applied, the default is 'yule', see [2]. Alternatively the
            birth-death process ('BDP', [3]) is available.
        innovation: If True, use the innovation model, see [4], to sample a lineage for the next
            speciation event. The default is False, in which case the lineage is chosen uniformly at
            random among the currently existing lineages.
        birth_rate: The birth rate for models such as 'yule' and 'BDP'.
        death_rate: The death rate for models such as 'BDP'.
        episodes: The episodes for the model 'EBDP'.
        contraction_probability: Probability that an inner edge is contracted. The default is 0.0,
            in which case the tree is binary. Only one of this parameter and
            'contraction_proportion' may be non-zero.
        contraction_proportion: The proportion of inner edges to be contracted. The default is 0.0,
            in which case the tree is binary. Only one of this parameter and
            'contraction_probability' may be non-zero.
        contraction_bias: Specifies whether shorter edges, i.e., with a smaller difference t of the
            time stamps, have a higher probability to be contracted. Only relevant if
            'contraction_proportion' > 0.0. The default is False, in which case all edges have the
            same probability to be contracted. The options 'inverse' and 'exponential' mean that an
            edge is sampled weighted by 1/(a * t) or e^(-a * t), respectively, where a is a
            user-defined factor.
        bias_strength: Intensity factor for preferring shorter edges to be contracted. The default
            is 1.0.

    Raises:
        ValueError: If unknown or invalid parameter values are passed.

    Returns:
        The simulated species tree.

    References:
    ----------
        .. [1] G. U. Yule. A mathematical theory of evolution, based on the conclusions of Dr. J.
           C. Willis, F. R. S. In: Phil. Trans. R. Soc. Lond. B, 1924, 213, 21-87.
           doi:10.1098/rstb.1925.0002.
        .. [2] D. G. Kendall. On the Generalized "Birth-and-Death" Process. In: Ann. Math. Statist.
           1948, 19, 1-15. doi:10.1214/aoms/1177730285.
        .. [3] T. Stadler. Simulating trees with a fixed number of extant species. In: Syst. Biol.
           2011, 60, 676-684. doi:10.1093/sysbio/syr029.
        .. [4] S. Keller-Schmidt, K. Klemm. A model of macroevolution as a branching process based
           on innovations. In: Adv. Complex Syst.,2012, 15, 1250043. doi:10.1142/S0219525912500439.
    """
    # parameter checking
    if not isinstance(age, (float, int)) or age <= 0.0:
        raise ValueError("age must be a number >0")
    elif isinstance(age, int):
        age = float(age)

    if not isinstance(n, int) or n < 0:
        raise ValueError("n must be an int >=0")
    elif n == 0:
        return Tree(None)

    if not isinstance(model, str):
        raise ValueError("model must be of type 'str'")

    # main simulation algorithm
    if model.lower() == "yule":
        tree = _yule_n_age(n, age, birth_rate, innovation)
    elif model.upper() == "BDP":
        tree = _BDP_n_age(n, age, birth_rate, death_rate, innovation)
    else:
        raise ValueError(f"model '{model}' is not available")

    # make tree non_binary by random contraction of edges
    nonbinary(tree, inplace=True, **kwargs)

    # assign the distance attribute to all vertices
    distance_from_timing(tree)

    return tree


def nonbinary(
    tree: Tree,
    contraction_probability: float = 0.0,
    contraction_proportion: float = 0.0,
    contraction_bias: bool | str = False,
    bias_strength: float = 1.0,
    inplace: bool = False,
    **kwargs,
) -> Tree:
    """Introduce multifurcation into a tree by contraction of inner edges.

    Args:
        tree: The tree whose edges shall be contracted.
        contraction_probability: Probability that an inner edge is contracted; results in non-binary
            tree; default is 0.0. Only one of this parameter and 'contraction_proportion' may be
            non-zero.
        contraction_proportion: The proportion of inner edges to be contracted; results in
            non-binary tree; default is 0.0. Only one of this parameter and
            'contraction_probability' may be non-zero.
        contraction_bias: Specifies whether shorter edges, i.e., with a smaller difference t of the
            time stamps, have a higher probability to be contracted. Only relevant if
            'contraction_proportion' > 0.0. The default is False, in which case all edges have the
            same probability to be contracted. The options 'inverse' and 'exponential' mean that an
            edge is sampled weighted by 1/(a * t) or e^(-a * t), respectively, where a is a
            user-defined factor.
        bias_strength: Intensity factor for preferring shorter edges to be contracted.
        inplace: If True, the edges are contracted in the original tree instance; otherwise a copy
            of the tree in which edges are contracted is returned.

    Returns:
        The tree whose edges have been contracted according to the parameters.

    Raises:
        ValueError: If both contraction_probability and contraction_proportion are non-zero or do
            not lie in the interval [0, 1].
    """
    if contraction_probability != 0.0 and contraction_proportion != 0.0:
        raise ValueError(
            "only one of parameters 'contraction_probability' and 'contraction_proportion' may be "
            "non-zero"
        )

    if contraction_probability < 0.0 or contraction_probability > 1.0:
        raise ValueError("contraction probability must be in [0.0, 1.0]")

    if contraction_proportion < 0.0 or contraction_proportion > 1.0:
        raise ValueError("contraction proportion must be in [0.0, 1.0]")

    if not isinstance(bias_strength, (int, float)) or bias_strength <= 0.0:
        raise ValueError("factor for contraction bias must be > 0")

    if not inplace:
        tree = tree.copy()

    edges = None
    if contraction_probability > 0.0:
        edges = _select_edges_by_probability(
            tree, contraction_probability, exclude_planted_edge=True
        )
    elif contraction_proportion > 0.0:
        edges = _select_edges_by_proportion(
            tree, contraction_proportion, contraction_bias, bias_strength, exclude_planted_edge=True
        )

    if edges:
        tree.contract(edges)
        distance_from_timing(tree)

    return tree


# --------------------------------------------------------------------------------------------------
#                           auxiliary functions for tree manipulation
# --------------------------------------------------------------------------------------------------


def rescale(tree: Tree, height: float, inplace: bool = True) -> Tree:
    """Rescale the time stamps of a tree.

    The time stamps of all vertices are multiplied by a scaling factor such that the height of the
    tree is equal to the specified height. The height of the tree is defined as the time stamp of
    the root vertex.

    Args:
        tree: The tree whose time stamps shall be rescaled.
        height: The new height of the tree.
        inplace: If True, the time stamps are rescaled in the original tree instance; otherwise a
            copy of the tree with rescaled time stamps is returned.

    Returns:
        The tree with rescaled time stamps.
    """
    if not inplace:
        tree = tree.copy()

    old_height = tree.root.tstamp

    # not available for trees that only consist of a root
    if old_height <= 0.0:
        raise RuntimeError(f"cannot rescale tree of height {old_height}")

    scaling_factor = height / old_height

    for v in tree.preorder():
        v.tstamp *= scaling_factor

    return tree


def _select_edges_by_probability(
    tree: Tree,
    p: float,
    exclude_planted_edge: bool = True,
) -> list[tuple[TreeNode, TreeNode]]:
    """Select edges for contraction with a given probability.

    Args:
        tree: The tree whose edges shall be selected.
        p: The probability with which an edge is selected for contraction.
        exclude_planted_edge: If True, the edge between the planted root and the canonical root is
            not considered for selection.

    Returns:
        A list of edges that have been selected for contraction.
    """
    edges: list[tuple[TreeNode, TreeNode]] = []

    for u, v in tree.inner_edges():
        if exclude_planted_edge and (u is tree.root) and len(u.children) == 1:
            continue

        if random.random() < p:
            edges.append((u, v))

    return edges


def _select_edges_by_proportion(
    tree: Tree,
    p: float,
    weighting: str,
    weighting_factor: float,
    exclude_planted_edge: bool = True,
) -> list[tuple[TreeNode, TreeNode]]:
    """Select a given proportion of edges for contraction.

    Optionally, a bias towards shorter edges can be applied.

    Args:
        tree: The tree whose edges shall be selected.
        p: The proportion of edges to be selected for contraction.
        weighting: The weighting scheme for biased selection ("inverse" or "exponential").
        weighting_factor: The factor controlling the strength of the bias.
        exclude_planted_edge: If True, the edge between the planted root and the canonical root is
            not considered for selection.

    Returns:
        A list of edges that have been selected for contraction.
    """
    edges: list[tuple[TreeNode, TreeNode]] = [
        e
        for e in tree.inner_edges()
        if not (exclude_planted_edge and (e[0] is tree.root) and len(e[0].children) == 1)
    ]

    # number of edges to be sampled
    k = round(p * len(edges))

    if not weighting:
        return random.choices(edges, k=k)
    else:
        distances = [abs(u.tstamp - v.tstamp) for u, v in edges]
        if weighting == "inverse":
            weights = 1 / (weighting_factor * np.asarray(distances))
        elif weighting == "exponential":
            weights = np.exp(-weighting_factor * np.asarray(distances))
        else:
            raise ValueError(f"unknown mode for weighted sampling: {weighting}")

        return random.choices(edges, weights=weights, k=k)


def assign_losses(tree: Tree, proportion: float) -> None:
    """Randomly assigns a specified proportion of leaves as losses.

    Args:
        tree: The tree whose leaves shall be assigned as losses.
        proportion: The proportion of leaves to be assigned as losses.
    """
    for leaf in tree.random_leaves(proportion):
        leaf.event = "L"


# --------------------------------------------------------------------------------------------------
#                            forward simulation of species trees
# --------------------------------------------------------------------------------------------------


class _ForwardLineageSampler:
    """Sample lineages for speciation events (forward simulation).

    Optionally uses the innovation model [1] for sampling the lineage in which the next speciation
    occurs.

    References:
        .. [1] S. Keller-Schmidt, K. Klemm. A model of macroevolution as a branching process based
           on innovations. In: Adv. Complex Syst.,2012, 15, 1250043. doi:10.1142/S0219525912500439.
    """

    def __init__(self, innovation: bool) -> None:
        """Initialize the forward lineage sampler.

        Args:
            innovation: If True, use the innovation model for sampling the lineage in which the next
                speciation occurs; otherwise, the lineage is chosen uniformly at random among the
                currently existing lineages.
        """
        self.innovation = innovation

        root = TreeNode(label=0, event=None, tstamp=0.0)
        self.tree = Tree(root)

        # lineages (branch id, parent node, feature set index)
        self.lineages = [(1, root, 0)]
        self.node_counter = 2

        if self.innovation:
            self.feature_counter = 0
            self.feature_sets = [frozenset([])]

            # feature set --> lineage index
            self.set_to_species = {frozenset([]): 0}

    def speciation(self, t: float) -> TreeNode:
        """Simulate a speciation event at time t.

        A new lineage is sampled for the speciation event according to the innovation model if
        self.innovation is True; otherwise, the lineage is chosen uniformly at random among the
        currently existing lineages. The new lineage is added to the list of lineages and a
        speciation event is recorded in the tree.

        Args:
            t: The time at which the speciation event occurs.

        Returns:
            The node in the tree that corresponds to the speciation event.
        """
        if self.innovation:
            loss_candidates = set()  # species for which loss of a feature can trigger a speciation
            for s in self.feature_sets:
                for f in s:
                    if s - {f} not in self.set_to_species:
                        loss_candidates.add(s)

            if loss_candidates:  # speciation by loss of feature
                loss_candidates = list(loss_candidates)
                while True:
                    s = random.choice(loss_candidates)
                    f = np.random.randint(self.feature_counter)
                    new_s = s - {f}
                    if new_s not in self.set_to_species:
                        break

            else:  # speciation by gain of feature (innovation)
                s = random.choice(self.feature_sets)
                new_feature = self.feature_counter
                self.feature_counter += 1
                new_s = s.union([new_feature])

            new_set_index = len(self.feature_sets)
            self.feature_sets.append(new_s)
            i = self.set_to_species[s]
            self.set_to_species[new_s] = len(self.lineages)

        else:
            i = np.random.randint(len(self.lineages))
            new_set_index = 0

        lin_id, parent, set_index = self.lineages[i]
        spec_node = TreeNode(label=lin_id, event="S", tstamp=t)
        parent.add_child(spec_node)

        self.lineages[i] = (self.node_counter, spec_node, set_index)
        self.lineages.append((self.node_counter + 1, spec_node, new_set_index))
        self.node_counter += 2

        return spec_node

    def extinction(self, t: float) -> TreeNode:
        """Simulate an extinction event at time t.

        A lineage is sampled uniformly at random among the currently existing lineages and removed
        from the list of lineages. An extinction event is recorded in the tree.

        Args:
            t: The time at which the extinction event occurs.

        Returns:
            The node in the tree that corresponds to the extinction event.
        """
        return self.lineage_extinction(np.random.randint(len(self.lineages)), t)

    def mass_extinction(self, surviving_rate: float, t: float) -> None:
        """Simulate a mass extinction event at time t.

        A number of lineages are removed uniformly at random according to the surviving rate.

        Args:
            surviving_rate: The proportion of lineages that survive the mass extinction event.
            t: The time at which the mass extinction event occurs.
        """
        no_of_losses = round((1 - surviving_rate) * len(self.lineages))

        # indices in decreasing order so that indices in self.lineages remain valid upon removal
        chosen_losses = sorted(
            np.random.choice(len(self.lineages), replace=False, size=no_of_losses), reverse=True
        )

        for i in chosen_losses:
            self.lineage_extinction(i, t)

    def lineage_extinction(self, i: int, t: float) -> TreeNode:
        """Simulate the extinction of a lineage with index i at time t.

        Args:
            i: The index of the lineage to be extinct.
            t: The time at which the extinction event occurs.

        Returns:
            The node in the tree that corresponds to the extinction event.
        """
        lin_id, parent, j = self.lineages[i]
        loss_node = TreeNode(label=lin_id, event="L", tstamp=t)
        parent.add_child(loss_node)
        self.lineages.pop(i)

        if self.innovation:
            feature_set = self.feature_sets[j]
            self.feature_sets.pop(j)
            del self.set_to_species[feature_set]

        return loss_node

    def finalize(self, time: float) -> Tree:
        """Finalize the tree by adding speciation events for all remaining lineages.

        Args:
            time: The time at which the finalization occurs.

        Returns:
            The finalized tree with speciation events added for all remaining lineages and time
            stamps adjusted.
        """
        for lin_id, parent, _ in self.lineages:
            parent.add_child(TreeNode(label=lin_id, event="S", tstamp=time))

        for v in self.tree.preorder():
            v.tstamp = abs(v.tstamp - time)

        return self.tree


def _yule_n(n: int, birth_rate: float, innovation: bool) -> Tree:
    """Simulate a Yule tree with n leaves.

    Args:
        n: The number of leaves in the resulting tree.
        birth_rate: The birth rate for the Yule process.
        innovation: If True, use the innovation model for sampling the lineage in which the next
            speciation occurs; otherwise, the lineage is chosen uniformly at random among the
            currently existing lineages.

    Returns:
        The simulated tree with n leaves.
    """
    if birth_rate is None:
        birth_rate = 1.0
    elif birth_rate <= 0.0:
        raise ValueError("birth rate must be >0")

    fls = _ForwardLineageSampler(innovation)
    t = 0.0

    while len(fls.lineages) < n:
        rate = len(fls.lineages) * birth_rate
        t += np.random.exponential(1 / rate)
        fls.speciation(t)

    # add length for pendant lineages (cf. Hartmann et al. 2010)
    t += np.random.exponential(1 / rate)

    return fls.finalize(t)


def _yule_age(age: float, birth_rate: float, innovation: bool) -> Tree:
    """Simulate a Yule tree up to a given age.

    Args:
        age: The age of the resulting tree.
        birth_rate: The birth rate for the Yule process.
        innovation: If True, use the innovation model for sampling the lineage in which the next
            speciation occurs; otherwise, the lineage is chosen uniformly at random among the
            currently existing lineages.

    Returns:
        The simulated tree up to the given age.
    """
    if birth_rate is None:
        birth_rate = 1.0
    elif birth_rate <= 0.0:
        raise ValueError("birth rate must be >0")

    fls = _ForwardLineageSampler(innovation)
    t = 0.0

    while t < age:
        rate = len(fls.lineages) * birth_rate
        t += np.random.exponential(1 / rate)
        if t >= age:
            break
        fls.speciation(t)

    return fls.finalize(age)


def _yule_n_age(n: int, age: float, birth_rate: float, innovation: bool) -> Tree:
    """Simulate a Yule tree with n leaves and of the specified age.

    Args:
        n: The number of leaves in the resulting tree.
        age: The age of the resulting tree.
        birth_rate: The birth rate for the Yule process.
        innovation: If True, use the innovation model for sampling the lineage in which the next
            speciation occurs; otherwise, the lineage is chosen uniformly at random among the
            currently existing lineages.

    Returns:
        The simulated tree with n leaves and of the specified age.
    """
    return _BDP_n_age(n, age, birth_rate, 0.0, innovation)


def _BDP_age(
    age: float,
    innovation: bool,
    birth_rate: float = None,
    death_rate: float = None,
    **kwargs,
) -> Tree:
    """Simulate a birth-death tree up to a given age.

    Args:
        age: The age of the resulting tree.
        innovation: If True, use the innovation model for sampling the lineage in which the next
            speciation occurs; otherwise, the lineage is chosen uniformly at random among the
            currently existing lineages.
        birth_rate: The birth rate for the birth-death process.
        death_rate: The death rate for the birth-death process.

    Returns:
        The simulated tree up to the given age.
    """
    # ignore potentially supplied 'episodes' argument in kwargs
    episodes = _EBDP_age_check_episodes(birth_rate=birth_rate, death_rate=death_rate)

    return _EBDP_age_forward(age, episodes, innovation)


def _BDP_n_age(
    n: int,
    age: float,
    birth_rate: float,
    death_rate: float,
    innovation: bool,
) -> Tree:
    """Simulate a birth-death tree with n leaves and of the specified age.

    The tree is sampled under the birth-death process conditioned on the number n of surviving
    species and with a given age [1]. The resulting tree does not contain loss leaves even if
    the specified death rate is > 0 as a consequence of the tree sampling method.

    Args:
        n: The number of leaves in the resulting tree.
        age: The age of the resulting tree.
        birth_rate: The birth rate for the birth-death process.
        death_rate: The death rate for the birth-death process.
        innovation: If True, use the innovation model for sampling the lineage in which the next
            speciation occurs; otherwise, the lineage is chosen uniformly at random among the
            currently existing lineages.

    Returns:
        The simulated tree with n leaves and of the specified age.

    Raises:
        ValueError: If birth_rate is not > 0 or if death_rate is not >= 0.

    References:
        .. [1] T. Stadler. Simulating trees with a fixed number of extant species. In: Syst. Biol.
           2011, 60, 676-684. doi:10.1093/sysbio/syr029.
    """
    if birth_rate is None:
        birth_rate = 1.0
    elif birth_rate <= 0.0:
        raise ValueError("birth rate must be >0")

    if death_rate is None:
        death_rate = 0.0
    elif death_rate < 0.0:
        raise ValueError("death rate must be >=0")

    fls = _ForwardLineageSampler(innovation)
    inner_vertices = []

    while len(fls.lineages) < n:
        inner_vertices.append(fls.speciation(0.0))
    fls.finalize(0.0)

    tree = fls.tree
    tree.root.tstamp = age

    # the following code is adapted from TreeSim (Stadler 2011, see [1])

    spec_times = []
    rho = 1.0  # proportion of sampled extant species (for future extension of the function)
    random_uniform = np.random.random(len(inner_vertices))
    for r in random_uniform:
        lamb1 = rho * birth_rate
        mu1 = death_rate - birth_rate * (1 - rho)

        if birth_rate > death_rate:
            t = (
                1
                / (lamb1 - mu1)
                * np.log(
                    (
                        lamb1
                        - mu1 * np.exp((-lamb1 + mu1) * age)
                        - mu1 * (1 - np.exp((-lamb1 + mu1) * age)) * r
                    )
                    / (
                        lamb1
                        - mu1 * np.exp((-lamb1 + mu1) * age)
                        - lamb1 * (1 - np.exp((-lamb1 + mu1) * age)) * r
                    )
                )
            )
        else:
            t = -((age * r) / (-1 - birth_rate * rho * age + birth_rate * rho * age * r))

        spec_times.append(t)

    spec_times.sort(reverse=True)

    for v, t in zip(inner_vertices, spec_times):
        v.tstamp = t

    return tree


def _EBDP_age(
    age: float,
    innovation: bool,
    birth_rate: float | None = None,
    death_rate: float | None = None,
    episodes: list[tuple[float, float, float, float]] | None = None,
    **kwargs,
) -> Tree:
    """Simulate an episodic birth-death tree up to a given age.

    Args:
        age: The age of the resulting tree.
        innovation: If True, use the innovation model for sampling the lineage in which the next
            speciation occurs; otherwise, the lineage is chosen uniformly at random among the
            currently existing lineages.
        birth_rate: The birth rate for the episodic birth-death process.
        death_rate: The death rate for the episodic birth-death process.
        episodes: The episodes for the episodic birth-death process.

    Returns:
        The simulated tree up to the given age.
    """
    episodes = _EBDP_age_check_episodes(
        birth_rate=birth_rate, death_rate=death_rate, episodes=episodes
    )

    return _EBDP_age_forward(age, episodes, innovation)


def _EBDP_age_check_episodes(
    birth_rate: float | None = None,
    death_rate: float | None = None,
    episodes: list[tuple[float, float, float, float]] | None = None,
    **kwargs,
) -> list[tuple[float, float, float, float]]:
    """Check the validity of the episodes for the episodic birth-death process.

    The episodes parameter is preferred over birth_rate and death_rate parameters. If episodes is
    not provided, birth_rate and death_rate will be used to create a single episode.

    Args:
        birth_rate: The birth rate for the episodic birth-death process.
        death_rate: The death rate for the episodic birth-death process.
        episodes: The episodes for the episodic birth-death process.

    Returns:
        The validated episodes for the episodic birth-death process.

    Raises:
        ValueError: If the episodes are not valid, or if birth_rate or death_rate are invalid when
            episodes is not provided.
    """
    # episodes parameter is preferred
    if episodes is not None:
        if len(episodes) == 0:
            raise ValueError("list of episodes must not be empty")

        for i in range(len(episodes)):
            if len(episodes[i]) != 4:
                raise ValueError(
                    "all episodes must contain 4 values: birth rate, death rate, proportion of "
                    "survivors, time stamp (from recent time as 0)"
                )

            birth_rate, death_rate, rho, t = episodes[i]
            if i == 0 and t != 0.0:
                raise ValueError("first episode must be at t=0.0")
            elif i > 0 and episodes[i - 1][3] >= t:
                print(episodes[i - 1][3], t)
                raise ValueError("episodes must be in correct temporal order")

            if birth_rate < 0.0 or death_rate < 0.0:
                raise ValueError("birth and death rate must be >=0 in all episodes")

            if rho <= 0.0 or rho > 1.0:
                raise ValueError("proportion of survivors must be in (0.0, 1.0]")

        return episodes

    elif birth_rate is not None:
        if birth_rate < 0.0 or (death_rate is not None and death_rate < 0.0):
            raise ValueError("birth and death rate must be >=0")

        if death_rate is None:
            # default death rate = 0.0
            return [(birth_rate, 0.0, 1.0, 0.0)]
        else:
            return [(birth_rate, death_rate, 1.0, 0.0)]

    else:
        if death_rate:
            raise ValueError("birth rate (>=0) must be specified if death rate is supplied")

        # default birth rate = 1.0 and death rate = 0.0
        return [(1.0, 0.0, 1.0, 0.0)]


def _EBDP_age_forward(
    age: float,
    episodes: list[tuple[float, float, float, float]],
    innovation: float,
) -> Tree:
    """Episodic birth-death process (EBDP), forward algorithm with maximum age.

    Simulate a tree under the episodic birth-death process up to a given age.

    Args:
        age: The age of the resulting tree.
        episodes: The episodes for the episodic birth-death process.
        innovation: If True, use the innovation model for sampling the lineage in which the next
            speciation occurs; otherwise, the lineage is chosen uniformly at random among the
            currently existing lineages.

    Returns:
        The simulated tree up to the given age.
    """
    fls = _ForwardLineageSampler(innovation)
    t = 0.0  # current time (forward simulation, is reversed in the end)
    i = 0  # current episode

    # may lead to extinction of the single branch at time t=0
    fls.mass_extinction(episodes[i][2], episodes[i][3])

    while t < age:
        birth_rate, death_rate, *_ = episodes[i]

        rate = len(fls.lineages) * (birth_rate + death_rate)
        waiting_time = np.random.exponential(1 / rate) if rate > 0.0 else float("inf")

        if i + 1 < len(episodes) and t + waiting_time >= episodes[i + 1][3]:
            fls.mass_extinction(episodes[i + 1][2], episodes[i + 1][3])
            t = episodes[i + 1][3]
            i += 1

        elif t + waiting_time >= age:
            break

        else:
            t += waiting_time

            if birth_rate > np.random.uniform(low=0.0, high=birth_rate + death_rate):
                fls.speciation(t)
            else:
                fls.extinction(t)

    return fls.finalize(age)


# --------------------------------------------------------------------------------------------------
#             backward simulation of species trees conditioned on n extant species
# --------------------------------------------------------------------------------------------------


def _BDP_n(
    n: int,
    birth_rate: float | None = None,
    death_rate: float | None = None,
    **kwargs,
) -> Tree:
    """Simulate a birth-death tree with n leaves under the backward algorithm.

    Args:
        n: The number of leaves in the resulting tree.
        birth_rate: The birth rate for the birth-death process.
        death_rate: The death rate for the birth-death process.

    Returns:
        The simulated tree with n leaves.
    """
    # ignore potentially supplied 'episodes' argument
    episodes = _EBDP_check_episodes(birth_rate=birth_rate, death_rate=death_rate)

    return _EBDP_backward(n, episodes)


def _EBDP_n(
    n: int,
    birth_rate: float | None = None,
    death_rate: float | None = None,
    episodes: list[tuple[float, float, float, float]] | None = None,
    **kwargs,
) -> Tree:
    """Simulate an episodic birth-death tree with n leaves under the backward algorithm.

    Args:
        n: The number of leaves in the resulting tree.
        birth_rate: The birth rate for the episodic birth-death process.
        death_rate: The death rate for the episodic birth-death process.
        episodes: The episodes for the episodic birth-death process.

    Returns:
        The simulated tree with n leaves.
    """
    episodes = _EBDP_check_episodes(birth_rate=birth_rate, death_rate=death_rate, episodes=episodes)

    return _EBDP_backward(n, episodes)


def _EBDP_check_episodes(
    birth_rate: float | None = None,
    death_rate: float | None = None,
    episodes: list[tuple[float, float, float, float]] | None = None,
    **kwargs,
) -> list[tuple[float, float, float, float]]:
    """Check the validity of the episodes for the episodic birth-death process.

    The episodes parameter is preferred over birth_rate and death_rate parameters. If episodes is
    not provided, birth_rate and death_rate will be used to create a single episode.

    Args:
        birth_rate: The birth rate for the episodic birth-death process.
        death_rate: The death rate for the episodic birth-death process.
        episodes: The episodes for the episodic birth-death process.

    Returns:
        The validated episodes for the episodic birth-death process.

    Raises:
        ValueError: If the episodes are not valid, or if birth_rate or death_rate are invalid when
            episodes is not provided.
    """
    # episodes parameter is preferred
    if episodes is not None:
        if len(episodes) == 0:
            raise ValueError("list of episodes must not be empty")

        for i in range(len(episodes)):
            if len(episodes[i]) != 4:
                raise ValueError(
                    "all episodes must contain 4 values: birth rate, "
                    "death rate, proportion of survivors, time stamp "
                    "(from recent time as 0)"
                )

            birth_rate, death_rate, rho, t = episodes[i]
            if i == 0 and t != 0.0:
                raise ValueError("first episode must be at t=0.0")
            elif i > 0 and episodes[i - 1][3] >= t:
                print(episodes[i - 1][3], t)
                raise ValueError("episodes must be in correct temporal order")

            if birth_rate <= 0.0 or birth_rate < death_rate:
                raise ValueError("birth rate must be >0 and >=death rate in all episodes")

            if rho <= 0.0 or rho > 1.0:
                raise ValueError("proportion of survivors must be in (0.0, 1.0]")

        return episodes

    elif birth_rate is not None:
        if birth_rate <= 0.0 or (death_rate and birth_rate < death_rate):
            raise ValueError("birth rate must be >0 and >=death rate")

        if death_rate and death_rate < 0.0:
            raise ValueError("death rate must be >=0")

        if death_rate is None:
            return [(birth_rate, 0.0, 1.0, 0.0)]
        else:
            return [(birth_rate, death_rate, 1.0, 0.0)]

    else:
        if death_rate:
            raise ValueError("birth rate (>0) must be specified if death rate is supplied")

        return [(1.0, 0.0, 1.0, 0.0)]


def _EBDP_backward(
    n: int,
    episodes: list[tuple[float, float, float, float]],
    max_tries: int = 500,
) -> Tree:
    """Episodic birth-death process (EBDP).

    Simulate a tree under the episodic birth-death process conditioned on the number n of surviving
    species. The tree is sampled using the backward algorithm [1]. The resulting tree does not
    contain loss leaves even if the specified death rates are > 0 as a consequence of the tree
    sampling method.

    The algorithm is not guaranteed to return a tree after a finite number of simulations. The
    parameter `max_tries` is used to limit the number of simulation attempts.

    Args:
        n: The number of leaves in the resulting tree.
        episodes: The episodes for the episodic birth-death process.
        max_tries: The maximum number of simulations to be performed in order to return a tree
            that satisfies the conditioning on n extant species; default is 500.

    Returns:
        The simulated tree with n leaves.

    Raises:
        ValueError: If the episodes are not valid.
        RuntimeError: If a tree with n leaves could not be returned after max_tries simulations.

    References:
        .. [1] T. Stadler. Simulating trees with a fixed number of extant species. In: Syst. Biol.
           2011, 60, 676-684. doi:10.1093/sysbio/syr029.
    """
    birth_inv_sum = sum([1 / episodes[i][0] for i in range(len(episodes))])

    for _ in range(max_tries):
        tree = None
        t = 0.0
        i = 0

        branches = [TreeNode(label=j, event="S", tstamp=t) for j in range(n)]
        id_counter = n

        while branches:
            birth_i, death_i, rho_i, t_i = episodes[i]

            losses_to_add = round(len(branches) / rho_i) - len(branches)
            for _ in range(losses_to_add):
                branches.append(TreeNode(label=id_counter, event="L", tstamp=t))
            id_counter += losses_to_add

            while branches:
                w = np.random.exponential(1 / ((birth_i + death_i) * len(branches)))

                if i + 1 < len(episodes) and t + w > episodes[i + 1][3]:
                    t = episodes[i + 1][3]
                    i += 1
                    break

                else:
                    t += w

                    if birth_i > np.random.uniform(low=0.0, high=birth_i + death_i):
                        # speciation event drawn
                        spec_node = TreeNode(label=id_counter, event="S", tstamp=t)
                        id_counter += 1
                        if len(branches) > 1:
                            k, m = np.random.choice(len(branches), 2, replace=False)
                            if k > m:
                                k, m = m, k
                            spec_node.add_child(branches[k])
                            spec_node.add_child(branches[m])
                            branches[k] = spec_node
                            branches.pop(m)
                        else:
                            spec_node.add_child(branches[0])
                            tree = Tree(spec_node)
                            # planted root is not a speciation event
                            spec_node.event = None
                            branches.clear()
                    else:
                        # extinction event drawn
                        branches.append(TreeNode(label=id_counter, event="L", tstamp=t))
                        id_counter += 1

        # return tree with the following probability
        if np.random.random() < (1 / birth_i) / birth_inv_sum:
            return tree

    raise RuntimeError(
        f"could not return a tree with {n} leaves after {max_tries} simulations under the "
        "episodic birth-death process with the specified parameters"
    )
