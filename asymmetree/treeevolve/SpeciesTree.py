# -*- coding: utf-8 -*-

"""
Simulation of dated species trees.
"""

import random, sys
import numpy as np

from tralda.datastructures.Tree import Tree, TreeNode

from asymmetree.tools.PhyloTreeTools import (delete_losses_and_contract,
                                             remove_planted_root,
                                             distance_from_timing)


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                         USER INTERFACE FUNCTION
# --------------------------------------------------------------------------

def species_tree_n(n, 
                   model='yule',
                   innovation=False,
                   planted=True,
                   remove_extinct=False,
                   **kwargs):
    """Simulates a species tree S with n leaves.
    
    Parameters
    ----------
    n : int
        Number of leaves in the resulting tree that correspond to extant
        species.
    model : str, optional
        Simulation model to be applied, the default is 'yule', see [1].
        Other available model are birth-death process ('BDP', [2]), and
        episodic birth-death process ('EBDP', [3]).
    innovation : bool, optional
        If True, use the innovation model, see [4], to sample a lineage for
        the next speciation event. Only available for the Yule model.
        The default is False, in which case the lineage is chosen uniformly at
        random among the currently existing lineages.
    planted : bool, optional
        Add a planted root that has the canonical root as its single neighbor; 
        default is True.
    remove_extinct : bool, optional
        Remove all branches that lead to extinctions, only relevant for some
        models; default is False.
    birth_rate : float, optional
        The birth rate for models such as 'yule' and 'BDP'.
    death_rate : float, optional
        The death rate for models such as 'BDP'.
    episodes : list, optional
        The episodes for the model 'EBDP'.
    contraction_probability : float, optional
        Probability that an inner edge is contracted. The default is 0.0, in
        which case the tree is binary. Only one of this parameter and
        'contraction_proportion' may be non-zero.
    contraction_proportion : float, optional
        The proportion of inner edges to be contracted. The default is 0.0, in
        which case the tree is binary. Only one of this parameter and 
        'contraction_probabilty' may be non-zero.
    contraction_bias : str, optional
        Specifies whether shorter edges, i.e., with a smaller difference t of
        the time stamps, have a higher probability to be contracted. Only
        relevant if 'contraction_proportion' > 0.0. The default is False,
        in which case all edges have the same probability to be contracted.
        The options 'inverse' and 'exponential' mean that an edge is sampled
        weighted by 1/t or e^(-t), respectively.
    bias_strength : float, optional
        Intensity factor for preferring shorter edges to be contracted. The
        default is 1.0.
        
    Raises
    ------
    ValueError
        If unknown parameter values are passed.
        
    Returns
    -------
    Tree
        The simulated species tree.
        
    References
    ----------
    .. [1] G. U. Yule.
       A mathematical theory of evolution, based on the conclusions of Dr. J. 
       C.Willis, F. R. S. 
       In: Phil. Trans. R. Soc. Lond. B, 1924, 213, 21–87.
       doi:10.1098/rstb.1925.0002.
    .. [2] D. G. Kendall.
       On the Generalized "Birth-and-Death" Process.
       In: Ann. Math. Statist. 1948, 19, 1–15.
       doi:10.1214/aoms/1177730285.
    .. [3] T. Stadler.
       Simulating trees with a fixed number of extant species.
       In: Syst. Biol. 2011, 60, 676–684.
       doi:10.1093/sysbio/syr029.
    .. [4] S. Keller-Schmidt, K. Klemm.
       A model of macroevolution as a branching process based on innovations.
       In: Adv. Complex Syst.,2012, 15, 1250043.
       doi:10.1142/S0219525912500439.
    """
    
    # parameter checking
    if not isinstance(n, int) or n < 0:
        raise ValueError('n must be an int >=0')
    elif n == 0:
        return Tree(None)
    
    if not isinstance(model, str):
        raise ValueError("model must be of type 'str'")
    
    # choice of the main simulation algorithm
    if model.lower() == 'yule':
        tree = _yule_n(n, kwargs.get('birth_rate'), innovation)
    elif model.upper() == 'BDP':
        tree = _BDP_n(n, **kwargs)
    elif model.upper() == 'EBDP':
        tree = _EBDP_n(n, **kwargs)
    else:
        raise ValueError(f"model '{model}' is not available")
        
    # remove extinct branches for models that include losses
    if remove_extinct and model.upper() in ('BDP', 'EBDP'):
        delete_losses_and_contract(tree, inplace=True)
        
    # remove planted edge for models that are planted by construction
    if not planted and model.upper() in ('YULE', 'BDP', 'EBDP'):
        remove_planted_root(tree, inplace=True)
    
    # make tree non_binary by random contraction of edges
    nonbinary(tree, inplace=True, **kwargs)
    
    # assign the distance attribute to all vertices
    distance_from_timing(tree)
        
    return tree


def species_tree_age(age, 
                     model='yule', 
                     innovation=False,
                     **kwargs):
    """Simulates a (planted) species tree S of the specified age.
    
    Parameters
    ----------
    age : float
        Simulation time, i.e., the time span from the root of the tree to the
        leaves that correspond to extant species.
    model : str, optional
        Simulation model to be applied, the default is 'yule', see [1].
        Other available model are birth-death process ('BDP', [2]), and
        episodic birth-death process ('EBDP', [3]).
    innovation : bool, optional
        If True, use the innovation model, see [4], to sample a lineage for
        the next speciation event. The default is False, in which case the
        lineage is chosen uniformly at random among the currently existing
        lineages.
    birth_rate : float, optional
        The birth rate for models such as 'yule' and 'BDP'.
    death_rate : float, optional
        The death rate for models such as 'BDP'.
    episodes : list, optional
        The episodes for the model 'EBDP'.
    contraction_probability : float, optional
        Probability that an inner edge is contracted. The default is 0.0, in
        which case the tree is binary. Only one of this parameter and
        'contraction_proportion' may be non-zero.
    contraction_proportion : float, optional
        The proportion of inner edges to be contracted. The default is 0.0, in
        which case the tree is binary. Only one of this parameter and 
        'contraction_probabilty' may be non-zero.
    contraction_bias : str, optional
        Specifies whether shorter edges, i.e., with a smaller difference t of
        the time stamps, have a higher probability to be contracted. Only
        relevant if 'contraction_proportion' > 0.0. The default is False,
        in which case all edges have the same probability to be contracted.
        The options 'inverse' and 'exponential' mean that an edge is sampled
        weighted by 1/(a * t) or e^(-a * t), respectively, where a is a
        user-defined factor.
    bias_strength : float, optional
        Intensity factor for preferring shorter edges to be contracted. The
        default is 1.0.
        
    Raises
    ------
    ValueError
        If unknown or invalid parameter values are passed.
        
    Returns
    -------
    Tree
        The simulated species tree.
        
    References
    ----------
    .. [1] G. U. Yule.
       A mathematical theory of evolution, based on the conclusions of Dr. J. 
       C.Willis, F. R. S. 
       In: Phil. Trans. R. Soc. Lond. B, 1924, 213, 21–87.
       doi:10.1098/rstb.1925.0002.
    .. [2] D. G. Kendall.
       On the Generalized "Birth-and-Death" Process.
       In: Ann. Math. Statist. 1948, 19, 1–15.
       doi:10.1214/aoms/1177730285.
    .. [3] T. Stadler.
       Simulating trees with a fixed number of extant species.
       In: Syst. Biol. 2011, 60, 676–684.
       doi:10.1093/sysbio/syr029.
    .. [4] S. Keller-Schmidt, K. Klemm.
       A model of macroevolution as a branching process based on innovations.
       In: Adv. Complex Syst.,2012, 15, 1250043.
       doi:10.1142/S0219525912500439.
    """
    
    # parameter checking
    if not isinstance(age, (float, int)) or age <= 0.0:
        raise ValueError('age must be a number >0')
    elif isinstance(age, int):
        age = float(age)
        
    if not isinstance(model, str):
        raise ValueError("model must be of type 'str'")
    
    # main simulation algorithm
    if model.lower() == 'yule':
        tree = _yule_age(age, kwargs.get('birth_rate'), innovation)
    elif model.upper() == 'BDP':
        tree = _BDP_age(age, innovation, **kwargs)
    elif model.upper() == 'EBDP':
        tree = _EBDP_age(age, innovation, **kwargs)
    else:
        raise ValueError(f"model '{model}' is not available")
        
    # make tree non_binary by random contraction of edges
    nonbinary(tree, inplace=True, **kwargs)
    
    # assign the distance attribute to all vertices
    distance_from_timing(tree)
        
    return tree


def species_tree_n_age(n, age, 
                       model='yule',
                       innovation=False,
                       birth_rate=1.0,
                       death_rate=0.0,
                       **kwargs):
    """Simulate a (planted) species tree S with n leaves and of the specified
    age.
    
    The tree is sampled under the Yule or BDP model conditioned on the number
    n of surviving species and with a given age [1].
    The resulting tree does not contain loss leaves even if the specified death
    rate is > 0 as a consequence of the tree sampling method.
    
    Parameters
    ----------
    n : int
        Number of leaves in the resulting tree that correspond to extant
        species.
    age : float
        Simulation time, i.e., the time span from the root of the tree to the
        leaves that correspond to extant species.
    model : str, optional
        Simulation model to be applied, the default is 'yule', see [2].
        Alternatively the birth-death process ('BDP', [3]) is available.
    innovation : bool, optional
        If True, use the innovation model, see [4], to sample a lineage for
        the next speciation event. The default is False, in which case the
        lineage is chosen uniformly at random among the currently existing
        lineages.
    birth_rate : float, optional
        The birth rate for models 'yule' and 'BDP'. The default is 1.0.
    death_rate : float, optional
        The death rate for model 'BDP'. The deafult is 0.0.
    contraction_probability : float, optional
        Probability that an inner edge is contracted. The default is 0.0, in
        which case the tree is binary. Only one of this parameter and
        'contraction_proportion' may be non-zero.
    contraction_proportion : float, optional
        The proportion of inner edges to be contracted. The default is 0.0, in
        which case the tree is binary. Only one of this parameter and 
        'contraction_probabilty' may be non-zero.
    contraction_bias : str, optional
        Specifies whether shorter edges, i.e., with a smaller difference t of
        the time stamps, have a higher probability to be contracted. Only
        relevant if 'contraction_proportion' > 0.0. The default is False,
        in which case all edges have the same probability to be contracted.
        The options 'inverse' and 'exponential' mean that an edge is sampled
        weighted by 1/(a * t) or e^(-a * t), respectively, where a is a
        user-defined factor.
    bias_strength : float, optional
        Intensity factor for preferring shorter edges to be contracted. The
        default is 1.0.
        
    Raises
    ------
    ValueError
        If unknown or invalid parameter values are passed.
        
    Returns
    -------
    Tree
        The simulated species tree.
        
    References
    ----------
    .. [1] T. Stadler.
       Simulating trees with a fixed number of extant species.
       In: Syst. Biol. 2011, 60, 676–684.
       doi:10.1093/sysbio/syr029.
    .. [2] G. U. Yule.
       A mathematical theory of evolution, based on the conclusions of Dr. J. 
       C.Willis, F. R. S. 
       In: Phil. Trans. R. Soc. Lond. B, 1924, 213, 21–87.
       doi:10.1098/rstb.1925.0002.
    .. [3] D. G. Kendall.
       On the Generalized "Birth-and-Death" Process.
       In: Ann. Math. Statist. 1948, 19, 1–15.
       doi:10.1214/aoms/1177730285.
    .. [4] S. Keller-Schmidt, K. Klemm.
       A model of macroevolution as a branching process based on innovations.
       In: Adv. Complex Syst.,2012, 15, 1250043.
       doi:10.1142/S0219525912500439.
    """
    
    # parameter checking
    if not isinstance(age, (float, int)) or age <= 0.0:
        raise ValueError('age must be a number >0')
    elif isinstance(age, int):
        age = float(age)
        
    if not isinstance(n, int) or n < 0:
        raise ValueError('n must be an int >=0')
    elif n == 0:
        return Tree(None)
        
    if not isinstance(model, str):
        raise ValueError("model must be of type 'str'")
    
    # main simulation algorithm
    if model.lower() == 'yule':
        tree = _yule_n_age(n, age, kwargs.get('birth_rate'), innovation)
    elif model.upper() == 'BDP':
        tree = _BDP_n_age(n, age,
                          kwargs.get('birth_rate'), kwargs.get('death_rate'),
                          innovation)
    else:
        raise ValueError(f"model '{model}' is not available")
        
    # make tree non_binary by random contraction of edges
    nonbinary(tree, inplace=True, **kwargs)
    
    # assign the distance attribute to all vertices
    distance_from_timing(tree)
        
    return tree


def nonbinary(tree,
              contraction_probability=0.0,
              contraction_proportion=0.0,
              contraction_bias=False,
              bias_strength=1.0,
              inplace=False,
              **kwargs):
    """Introduce multifurcation into a tree by contraction of inner edges.
    
    Parameters
    ----------
    tree : Tree
        The tree whose edges shall be contracted.
    contraction_probability : float, optional
        Probability that an inner edge is contracted; results in non-binary 
        tree; default is 0.0. Only one of this parameter and 
        'contraction_proportion' may be non-zero.
    contraction_proportion : float, optional
        The proportion of inner edges to be contracted; results in non-binary 
        tree; default is 0.0. Only one of this parameter and 
        'contraction_probabilty' may be non-zero.
    contraction_bias : str, optional
        Specifies whether shorter edges, i.e., with a smaller difference t of
        the time stamps, have a higher probability to be contracted. Only
        relevant if 'contraction_proportion' > 0.0. The default is False,
        in which case all edges have the same probability to be contracted.
        The options 'inverse' and 'exponential' mean that an edge is sampled
        weighted by 1/(a * t) or e^(-a * t), respectively, where a is a
        user-defined factor.
    bias_strength : float, optional
        Intensity factor for preferring shorter edges to be contracted. The
        default is 1.0.
    inplace : bool, optional
        If True, the edges are contracted in the original tree instance;
        otherwise a copy of the tree in which edges are contracted is returned.
    
    Returns
    -------
    Tree
        The tree whose edges have been contracted according to the parameters.
    
    Raises
    ------
    ValueError
        If both contraction_probability and contraction_proportion are
        non-zero or do not lie in the interval [0, 1].
    """
    
    if contraction_probability != 0.0 and contraction_proportion != 0.0:
        raise ValueError("only one of parameters 'contraction_probability' "\
                         "and 'contraction_proportion' may be non-zero")
    
    if contraction_probability < 0.0 or contraction_probability > 1.0:
        raise ValueError('contraction probability must be in [0.0, 1.0]')
    
    if contraction_proportion < 0.0 or contraction_proportion > 1.0:
        raise ValueError('contraction proportion must be in [0.0, 1.0]')
    
    if not isinstance(bias_strength, (int, float)) or \
        bias_strength <= 0.0:
        raise ValueError('factor for contraction bias must be > 0')
    
    if not inplace:
        tree = tree.copy()
    
    edges = False
    if contraction_probability > 0.0:
         edges = _select_edges_by_probability(tree, contraction_probability,
                                              exclude_planted_edge=True)
         
    elif contraction_proportion > 0.0:
        edges = _select_edges_by_proportion(tree, contraction_proportion,
                                            contraction_bias, bias_strength,
                                            exclude_planted_edge=True)
    if edges:
        tree.contract(edges)
        distance_from_timing(tree)
    
    return tree


# --------------------------------------------------------------------------
#                         AUXILIARY FUNCTIONS 
# --------------------------------------------------------------------------


def _rescale(tree, height, inplace=True):
    
    if not inplace:
        tree = tree.copy()
    
    old_height = tree.root.tstamp
    
    # not available for trees that only consist of a root
    if old_height <= 0.0:
        raise RuntimeError("cannot rescale tree of "\
                           "height '{}'".format(old_height))
        
    scaling_factor = height / old_height
    
    for v in tree.preorder():
        v.tstamp *= scaling_factor
    
    return tree
            
            
def _select_edges_by_probability(tree, p, exclude_planted_edge=True):
    
    edges = []
    
    for u, v in tree.inner_edges():
        
        if exclude_planted_edge and (u is tree.root) and len(u.children) == 1:
            continue
        
        if random.random() < p:
            edges.append((u,v))
    
    return edges


def _select_edges_by_proportion(tree, p, weighting, weighting_factor,
                                exclude_planted_edge=True):
    
    edges = [e for e in tree.inner_edges()
             if not (exclude_planted_edge and (e[0] is tree.root) and 
                     len(e[0].children) == 1)]
    
    # number of edges to be sampled
    k = round(p * len(edges))
    
    if not weighting:
        return random.choices(edges, k=k)
    else:
        distances = [abs(u.tstamp-v.tstamp) for u, v in edges]
        if weighting == 'inverse':
            weights = 1 / (weighting_factor * np.asarray(distances))
        elif weighting == 'exponential':
            weights = np.exp(-weighting_factor * np.asarray(distances))
        else:
            raise ValueError(f'unknown mode for weighted sampling: {weighting}')
        
        return random.choices(edges, weights=weights, k=k)



def assign_losses(tree, proportion):
    """Randomly assigns a specified proportion of leaves as losses."""
    
    for leaf in tree.random_leaves(proportion):
        leaf.event = 'L'


# --------------------------------------------------------------------------
#                             FORWARD MODELS
# --------------------------------------------------------------------------

class _ForwardLineageSampler:
    """Sample lineages for speciation events (forward simulation).
    
    Optionally uses the innovation model [1] for sampling the lineage in which
    the next speciation occurs.
    
    References
    ----------
    .. [1] S. Keller-Schmidt, K. Klemm.
       A model of macroevolution as a branching process based on innovations.
       In: Adv. Complex Syst.,2012, 15, 1250043.
       doi:10.1142/S0219525912500439.
    """
    
    def __init__(self, innovation):
        
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
    
    
    def speciation(self, t):
        
        if self.innovation:
            
            loss_candidates = set()     # species for which loss of a feature 
                                        # can trigger a speciation
            for s in self.feature_sets:
                for f in s:
                    if s - {f} not in self.set_to_species:
                        loss_candidates.add(s)
            
            if loss_candidates: # speciation by oss of feature
            
                loss_candidates = list(loss_candidates)
                while True:
                    s = random.choice(loss_candidates)
                    f = np.random.randint(self.feature_counter)
                    new_s = s - {f}
                    if new_s not in self.set_to_species:
                        break
                    
            else:               # speciation by gain of feature (innovation)
                
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
        spec_node = TreeNode(label=lin_id, event='S', tstamp=t)
        parent.add_child(spec_node)
        
        self.lineages[i] = (self.node_counter, spec_node, set_index)
        self.lineages.append((self.node_counter+1, spec_node, new_set_index))
        self.node_counter += 2
        
        return spec_node
    
    
    def extinction(self, t):
        
        return self.lineage_extinction(np.random.randint(len(self.lineages)), t)
        
    
    def mass_extinction(self, surviving_rate, t):
        
        no_of_losses = round((1-surviving_rate) * len(self.lineages))
        
        # indices in decreasing order so that indices in self.lineages remain
        # valid upon removal
        chosen_losses = sorted(np.random.choice(len(self.lineages), 
                                                replace=False,
                                                size=no_of_losses),
                               reverse=True)
        
        for i in chosen_losses:
            self.lineage_extinction(i, t)
            
    
    def lineage_extinction(self, i, t):
        
        lin_id, parent, j = self.lineages[i]
        loss_node = TreeNode(label=lin_id, event='L', tstamp=t)
        parent.add_child(loss_node)
        self.lineages.pop(i)
        
        if self.innovation:
            feature_set = self.feature_sets[j]
            self.feature_sets.pop(j)
            del self.set_to_species[feature_set]
        
        return loss_node
    
    
    def finalize(self, time):
        
        for lin_id, parent, _ in self.lineages:
            parent.add_child( TreeNode(label=lin_id, event='S', tstamp=time) )
        
        for v in self.tree.preorder():
            v.tstamp = abs(v.tstamp - time)
        
        return self.tree


def _yule_n(n, birth_rate, innovation):
    
    if birth_rate is None:
        birth_rate = 1.0
    elif birth_rate <= 0.0:
        raise ValueError('birth rate must be >0')
    
    fls = _ForwardLineageSampler(innovation)
    t = 0.0
    
    while len(fls.lineages) < n:
        rate = len(fls.lineages) * birth_rate
        t += np.random.exponential(1/rate)
        fls.speciation(t)
        
    # add length for pendant lineages (cf. Hartmann et al. 2010)
    t += np.random.exponential(1/rate)
    
    return fls.finalize(t)


def _yule_age(age, birth_rate, innovation):
    
    if birth_rate is None:
        birth_rate = 1.0
    elif birth_rate <= 0.0:
        raise ValueError("birth rate must be >0")
    
    fls = _ForwardLineageSampler(innovation)
    t = 0.0
    
    while t < age:
        rate = len(fls.lineages) * birth_rate
        t += np.random.exponential(1/rate)
        if t >= age: break
        fls.speciation(t)
    
    return fls.finalize(age)


def _yule_n_age(n, age, birth_rate, innovation):
    
    return _BDP_n_age(n, age, birth_rate, 0.0, innovation)
    

def _BDP_age(age, innovation, **kwargs):
    
    # remove potentially supplied 'episodes' argument
    episodes = _EBDP_age_check_episodes(birth_rate = kwargs.get('birth_rate'),
                                        death_rate = kwargs.get('death_rate'))
    
    return _EBDP_age_forward(age, episodes, innovation)


def _BDP_n_age(n, age, birth_rate, death_rate, innovation):
    
    if birth_rate is None:
        birth_rate = 1.0
    elif birth_rate <= 0.0:
        raise ValueError('birth rate must be >0')
    if death_rate is None:
        death_rate = 0.0
    elif death_rate < 0.0:
        raise ValueError('death rate must be >=0')
    
    fls = _ForwardLineageSampler(innovation)
    inner_vertices = []
    
    while len(fls.lineages) < n:
        inner_vertices.append(fls.speciation(0.0))
    fls.finalize(0.0)
    
    tree = fls.tree
    tree.root.tstamp = age
    
    # the following code is adapated from TreeSim (Stadler 2011)
    
    spec_times = []
    rho = 1.0           # proportion of sampled extant species (for future
                        # extension of the function)
    random_uniform = np.random.random(len(inner_vertices))
    for r in random_uniform:
        
        lamb1 = rho * birth_rate
        mu1   = death_rate - birth_rate * (1 - rho)
        
        if birth_rate > death_rate:
            t = 1 / (lamb1 - mu1) * \
                np.log((lamb1-mu1 * np.exp((-lamb1+mu1) * age) - \
                        mu1*(1-np.exp((-lamb1+mu1)*age)) * r ) /
                       (lamb1 - mu1 * np.exp((-lamb1+mu1) * age) - 
                        lamb1 * (1-np.exp((-lamb1+mu1) * age)) * r)
                      )
        else:
            t = - ((age * r) / 
                   (-1 - birth_rate * rho * age + birth_rate * rho* age * r ))
            
        spec_times.append(t)
        
    spec_times.sort(reverse=True)
    
    for v, t in zip(inner_vertices, spec_times):
        v.tstamp = t
    
    return tree


def _EBDP_age(age, innovation, **kwargs):
    
    episodes = _EBDP_age_check_episodes(**kwargs)
    
    return _EBDP_age_forward(age, episodes, innovation)


def _EBDP_age_check_episodes(**kwargs):
    
    birth_rate = kwargs.get('birth_rate')
    death_rate = kwargs.get('death_rate')
    episodes = kwargs.get('episodes')
    
    # episodes parameter is prefered
    if episodes is not None:
        
        if len(episodes) == 0:
            raise ValueError("list of episodes must not be empty")
        
        for i in range(len(episodes)):
            
            if len(episodes[i]) != 4:
                raise ValueError("all episodes must contain 4 values: birth rate, "\
                                 "death rate, proportion of survivors, time stamp "\
                                 "(from recent time as 0)")
            
            birth_rate, death_rate, rho, t = episodes[i]
            if i == 0 and t != 0.0:
                raise ValueError("first episode must be at t=0.0")
            elif i > 0 and episodes[i-1][3] >= t:
                print(episodes[i-1][3], t)
                raise ValueError("episodes must be in correct temporal order")
            
            if birth_rate < 0.0 or death_rate < 0.0:
                raise ValueError("birth and death rate must be >=0 "\
                                 "in all episodes")
                
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
            raise ValueError("birth rate (>=0) must be specified if death rate "\
                             "is supplied")
        
        # default birth rate = 1.0 and death rate = 0.0
        return [(1.0, 0.0, 1.0, 0.0)]
    

def _EBDP_age_forward(age, episodes, innovation):
    """Episodic birth–death process (EBDP), forward algorithm with max. age."""
    
    fls = _ForwardLineageSampler(innovation)
    t = 0.0  # current time (forward simulation, is reversed in the end)
    i = 0    # current episode
    
    # may lead to extinction of the single branch at time t=0
    fls.mass_extinction(episodes[i][2], episodes[i][3])
    
    while t < age:
        birth_rate, death_rate, *_ = episodes[i]
        
        rate = len(fls.lineages) * (birth_rate + death_rate)
        waiting_time = np.random.exponential(1/rate) if rate > 0.0 else float('inf')
        
        if i+1 < len(episodes) and t + waiting_time >= episodes[i+1][3]:
            fls.mass_extinction(episodes[i+1][2], episodes[i+1][3])
            t = episodes[i+1][3]
            i += 1
        
        elif t + waiting_time >= age:
            break
        
        else:
            t += waiting_time
            
            if birth_rate > np.random.uniform(low=0.0, 
                                              high=birth_rate+death_rate):
                fls.speciation(t)
            else:
                fls.extinction(t)
    
    return fls.finalize(age)


# --------------------------------------------------------------------------
#                            BACKWARD MODELS
# --------------------------------------------------------------------------


def _BDP_n(n, **kwargs):
    
    # remove potentially supplied 'episodes' argument
    episodes = _EBDP_check_episodes(birth_rate = kwargs.get('birth_rate'),
                                    death_rate = kwargs.get('death_rate'))
    
    return _EBDP_backward(n, episodes)


def _EBDP_n(n, **kwargs):
    
    episodes = _EBDP_check_episodes(**kwargs)
    
    return _EBDP_backward(n, episodes)


def _EBDP_check_episodes(**kwargs):
    
    birth_rate = kwargs.get('birth_rate')
    death_rate = kwargs.get('death_rate')
    episodes = kwargs.get('episodes')
    
    # episodes parameter is prefered
    if episodes is not None:
        
        if len(episodes) == 0:
            raise ValueError("list of episodes must not be empty")
        
        for i in range(len(episodes)):
            
            if len(episodes[i]) != 4:
                raise ValueError("all episodes must contain 4 values: birth rate, "\
                                 "death rate, proportion of survivors, time stamp "\
                                 "(from recent time as 0)")
            
            birth_rate, death_rate, rho, t = episodes[i]
            if i == 0 and t != 0.0:
                raise ValueError("first episode must be at t=0.0")
            elif i > 0 and episodes[i-1][3] >= t:
                print(episodes[i-1][3], t)
                raise ValueError("episodes must be in correct temporal order")
            
            if birth_rate <= 0.0 or birth_rate < death_rate:
                raise ValueError("birth rate must be >0 and >=death rate "\
                                 "in all episodes")
                
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
            raise ValueError("birth rate (>0) must be specified if death rate "\
                             "is supplied")
            
        return [(1.0, 0.0, 1.0, 0.0)]


def _EBDP_backward(n, episodes, max_tries=500):
    """Episodic birth–death process (EBDP).
    
    References
    ----------
    .. [1] T. Stadler.
       Simulating trees with a fixed number of extant species.
       In: Syst. Biol. 2011, 60, 676–684.
       doi:10.1093/sysbio/syr029.
    """
    
    birth_inv_sum = sum([1/episodes[i][0] for i in range(len(episodes))])
    
    for _ in range(max_tries):
        
        tree = None
        t = 0.0
        i = 0
        
        branches = [TreeNode(label=j, event='S', tstamp=t)
                    for j in range(n)]
        id_counter = n
        
        while branches:
            birth_i, death_i, rho_i, t_i = episodes[i]
            
            losses_to_add = round(len(branches) / rho_i) - len(branches)
            for j in range(losses_to_add):
                branches.append( TreeNode(label=id_counter,
                                          event='L', tstamp=t) )
            id_counter += losses_to_add
            
            while branches:
                w = np.random.exponential( 1 / ((birth_i + death_i) * 
                                                len(branches)) )
                
                if i+1 < len(episodes) and t + w > episodes[i+1][3]:
                    t = episodes[i+1][3]
                    i += 1
                    break
                
                else:
                    t += w
                    
                    if birth_i > np.random.uniform(low=0.0, 
                                                   high=birth_i+death_i):
                        # speciation event drawn
                        spec_node = TreeNode(label=id_counter, event='S',
                                             tstamp=t)
                        id_counter += 1
                        if len(branches) > 1:
                            k, l = np.random.choice(len(branches), 2,
                                                    replace=False)
                            if k > l:
                                k, l = l, k
                            spec_node.add_child(branches[k])
                            spec_node.add_child(branches[l])
                            branches[k] = spec_node
                            branches.pop(l)
                        else:
                            spec_node.add_child(branches[0])
                            tree = Tree(spec_node)
                            # planted root is not a speciation event
                            spec_node.event = None
                            branches.clear()
                    else:
                        # extinction event drawn
                        branches.append( TreeNode(label=id_counter,
                                                  event='L', tstamp=t) )
                        id_counter += 1
        
        # return tree with the following probability
        if np.random.random() < (1 / birth_i) / birth_inv_sum:
            
            return tree
        
    print(f'Could not return a tree after {max_tries} simulations',
          file=sys.stderr)
