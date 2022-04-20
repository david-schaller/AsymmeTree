# -*- coding: utf-8 -*-

"""
Simulation of dated species trees.
"""

import random, sys
import numpy as np

from tralda.datastructures.Tree import Tree, TreeNode

from asymmetree.tools.PhyloTreeTools import (delete_losses_and_contract,
                                             remove_planted_root)


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                         USER INTERFACE FUNCTION
# --------------------------------------------------------------------------

def simulate_species_tree(N, 
                          model='innovation',
                          planted=True,
                          remove_extinct=False,
                          rescale_to_height=False,
                          **kwargs):
    """Simulates a species tree S with N leaves.
    
    Parameters
    ----------
    N : int
        Number of leaves in the resulting tree that correspond to extant
        species.
    model : str, optional
        Simulation model to be applied, the default is 'innovation', see [1].
        Other available model are 'yule' [2], birth-death process ('BDP', [3]),
        and episodic birth-death process ('EBDP', [4]).
    planted : bool, optional
        Add a planted root that has the canonical root as its single neighbor; 
        default is True.
    remove_extinct : bool, optional
        Remove all branches that lead to extinctions, only relevant for some
        models; default is False.
    rescale_to_height : bool or float, optional
        Determines the final distance from the root to the (surviving) leaves,
        default is False, i.e., model-dependent.
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
    .. [1] S. Keller-Schmidt, K. Klemm.
       A model of macroevolution as a branching process based on innovations.
       In: Adv. Complex Syst.,2012, 15, 1250043.
       doi:10.1142/S0219525912500439.
    .. [2] G. U. Yule.
       A mathematical theory of evolution, based on the conclusions of Dr. J. 
       C.Willis, F. R. S. 
       In: Phil. Trans. R. Soc. Lond. B, 1924, 213, 21–87.
       doi:10.1098/rstb.1925.0002.
    .. [3] D. G. Kendall.
       On the Generalized "Birth-and-Death" Process.
       In: Ann. Math. Statist. 1948, 19, 1–15.
       doi:10.1214/aoms/1177730285.
    .. [4] T. Stadler.
       Simulating trees with a fixed number of extant species.
       In: Syst. Biol. 2011, 60, 676–684.
       doi:10.1093/sysbio/syr029.
    """
    
    # parameter checking
    if not isinstance(N, int) or N < 0:
        raise ValueError('N must be an int >=0')
    elif N == 0:
        return Tree(None)
    
    if not isinstance(model, str):
        raise ValueError("model must be of type 'str'")
        
    if (rescale_to_height is not False and
        (not isinstance(rescale_to_height, (int, float)) or
         rescale_to_height < 0.0)):
        raise ValueError('height must be a number >=0')
    elif rescale_to_height is not False and N == 1 and not planted:
        raise ValueError('rescaling is not applicable to unplanted trees '\
                         'with only one leaf')
    
    # choice of the main simulation algorithm
    if model.lower() in ('innovation', 'innovations'):
        tree = _innovation_model(N, planted)
    elif model.lower() == 'yule':
        tree = _yule(N, kwargs.get('birth_rate'))
    elif model.upper() == 'BDP':
        tree = _BDP(N, **kwargs)
    elif model.upper() == 'EBDP':
        tree = _EBDP(N, **kwargs)
    else:
        raise ValueError("model '{}' is not available".format(model))
        
    # remove extinct branches for models that include losses
    if remove_extinct and model.upper() in ('BDP', 'EBDP'):
        delete_losses_and_contract(tree, inplace=True)
        
    # remove planted edge for models that are planted by construction
    if not planted and model.upper() in ('YULE', 'BDP', 'EBDP'):
        remove_planted_root(tree, inplace=True)
    
    # make tree non_binary by random contraction of edges
    nonbinary(tree, inplace=True, **kwargs)
    
    # rescale to specified height
    if rescale_to_height is not False:
        _rescale(tree, rescale_to_height, inplace=True)
        
    return tree


def simulate_species_tree_age(age, model='yule', **kwargs):
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
        tree = _yule_age(age, kwargs.get('birth_rate'))
    elif model.upper() == 'BDP':
        tree = _BDP_age(age, **kwargs)
    elif model.upper() == 'EBDP':
        tree = _EBDP_age(age, **kwargs)
    else:
        raise ValueError("model '{}' is not available".format(model))
        
    # make tree non_binary by random contraction of edges
    nonbinary(tree, inplace=True, **kwargs)
        
    return tree


def nonbinary(tree,
              contraction_probability=0.0,
              contraction_proportion=0.0,
              contraction_bias=False,
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
        weighted by 1/t or e^(-t), respectively.
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
    
    if not inplace:
        tree = tree.copy()
    
    if contraction_probability > 0.0:
         edges = _select_edges_by_probability(tree, contraction_probability,
                                              exclude_planted_edge=True)
         tree.contract(edges)
         
    elif contraction_proportion > 0.0:
        edges = _select_edges_by_proportion(tree, contraction_proportion,
                                            contraction_bias,
                                            exclude_planted_edge=True)
        tree.contract(edges)
    
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
        v.dist   *= scaling_factor
    
    return tree
            
            
def _select_edges_by_probability(tree, p, exclude_planted_edge=True):
    
    edges = []
    
    for u, v in tree.inner_edges():
        
        if exclude_planted_edge and (u is tree.root) and len(u.children) == 1:
            continue
        
        if random.random() < p:
            edges.append((u,v))
    
    return edges


def _select_edges_by_proportion(tree, p, weighting, exclude_planted_edge=True):
    
    edges = [e for e in tree.inner_edges() if not
             (exclude_planted_edge and (e[0] is tree.root) and len(e[0].children) == 1)]
    
    # number of edges to be sampled
    k = round(p * len(edges))
    
    if not weighting:
        return random.choices(edges, k=k)
    else:
        distances = [abs(u.tstamp-v.tstamp) for u, v in edges]
        if weighting == 'inverse':
            weights = 1 / np.asarray(distances)
        elif weighting == 'exponential':
            weights = np.exp(-np.asarray(distances))
        else:
            raise ValueError(f'unknown mode for weighted sampling: {weighting}')
        
        return random.choices(edges, weights=weights, k=k)



def assign_losses(tree, proportion):
    """Randomly assigns a specified proportion of leaves as losses."""
    
    for leaf in tree.random_leaves(proportion):
        leaf.event = 'L'


# --------------------------------------------------------------------------
#                     SPECIES TREE MODEL FUNCTIONS 
# --------------------------------------------------------------------------


def _innovation_model(N, planted, ultrametric=True):
    """Builds a species tree S with N leaves with the innovation model.
    
    Parameters
    ----------
    N : int
        Number of extant species in the resulting tree.
    planted : bool
        Add a planted root that has the canonical root as its single neighbor.
    ultrametric : bool, optional
        If True, make tree ultrametric and rescale it to height 1.0, else all
        edges have length 1.0; the default is True.
    
    Returns
    -------
    Tree
        The simulated species tree.
    
    References
    ----------
    .. [1] S. Keller-Schmidt, K. Klemm.
       A model of macroevolution as a branching process based on innovations.
       In: Adv. Complex Syst.,2012, 15, 1250043.
       doi:10.1142/S0219525912500439.
    """
    
    tree = Tree(TreeNode(label=0, event='S'))
    tree.number_of_species = N
    node_counter = 1
    
    # planted tree (root is an implicit outgroup with outdegree = 1)
    if planted:
        root = TreeNode(label=1, event='S')
        tree.root.add_child(root)
        node_counter += 1
    else:
        root = tree.root
    
    features = [0]                  # set of available features
    species = {(0,): root}          # extant species
    
    while len(species) < N:
        
        loss_candidates = set()     # species for which loss of a feature
        for s in species.keys():    # can trigger a speciation
            for i in range(0, len(s)):
                    if s[:i] + s[i+1:] not in species:
                        loss_candidates.add(s)
        
        if not loss_candidates:     # INNOVATION EVENT
            s = random.choice(list(species))
            new_feature = len(features)
            
            new_s = s + (new_feature,)
            
            child1 = TreeNode(label=node_counter, event='S')
            species[s].add_child(child1)
            child2 = TreeNode(label=node_counter+1, event='S')
            species[s].add_child(child2)
            
            node_counter += 2
            
            species[s] = child1
            species[new_s] = child2
            features.append(new_feature)
            
        else:
            s = random.choice(list(loss_candidates))
            
            if len(s) > 1:
                feature_index = random.randint(0, len(s)-1)
            else:
                feature_index = 0
            
            new_s = s[:feature_index] + s[feature_index+1:]
            
            if new_s not in species:    # LOSS EVENT
                
                child1 = TreeNode(label=node_counter, event='S')
                species[s].add_child(child1)
                child2 = TreeNode(label=node_counter+1, event='S')
                species[s].add_child(child2)

                node_counter += 2
                
                species[s] = child1
                species[new_s] = child2
                
    if ultrametric:         
        simulate_timing(tree)
        distance_from_timing(tree)
    
    return tree


def simulate_timing(tree):
    """Simulates a timing for the tree.
    
    It is t(root) = 1 and t(x) = 0 for x in L(S)."""
    
    for v in tree.preorder():
        if not v.children:
            v.tstamp = 0.0
        elif not v.parent:
            v.tstamp = 1.0
        else:                               # random walk to a leaf
            pos = v                         # current position
            length = 0                      # path length |P|
            while pos.children:
                length += 1
                pos = pos.children[np.random.randint(len(pos.children))]
            v.tstamp = v.parent.tstamp * (1 - 2 * np.random.uniform() / (length+1))
            

def distance_from_timing(tree):
    """Adjusts all distances according to the time stamp difference."""
    
    if tree.root:
        tree.root.dist = 0.0
    
    for u, v in tree.edges():
        v.dist = abs(u.tstamp - v.tstamp)


def _reverse_time_stamps(tree):
    
    max_depth = 0.0
    for v in tree.preorder():
        max_depth = max(max_depth, v.tstamp)
            
    for v in tree.preorder():
        v.tstamp = abs(v.tstamp - max_depth)


def _yule(N, birth_rate):
    
    if birth_rate is None:
        birth_rate = 1.0
    elif birth_rate <= 0.0:
        raise ValueError("birth rate must be >0")
    
    tree = Tree(TreeNode(label=0, event='S', dist=0.0, tstamp=0.0))
    tree.number_of_species = N
    
    branches = [(1, tree.root)]
    forward_time = 0.0
    node_counter = 1
    
    while len(branches) < N:
        
        rate = len(branches) * birth_rate
        forward_time += np.random.exponential(1/rate)
        
        i = np.random.randint(len(branches))
        branch_id, parent = branches[i]
        spec_node = TreeNode(label=branch_id, event='S',
                             dist=forward_time-parent.tstamp,
                             tstamp=forward_time)
        parent.add_child(spec_node)
        branches[i] = (node_counter, spec_node)
        branches.append((node_counter+1, spec_node))
        node_counter += 2
        
    # add length for pendant branches (cf. Hartmann et al. 2010)
    forward_time += np.random.exponential(1/rate)
    
    # finalize the branches
    for branch_id, parent in branches:
        parent.add_child( TreeNode(label=branch_id, event='S',
                                   dist=forward_time-parent.tstamp,
                                   tstamp=forward_time) )
    
    _reverse_time_stamps(tree)
    
    return tree


def _yule_age(age, birth_rate):
    
    if birth_rate is None:
        birth_rate = 1.0
    elif birth_rate <= 0.0:
        raise ValueError("birth rate must be >0")
    
    tree = Tree(TreeNode(label=0, event='S', dist=0.0, tstamp=0.0))
    
    branches = [(1, tree.root)]
    forward_time = 0.0
    node_counter = 1
    
    while forward_time < age:
        
        rate = len(branches) * birth_rate
        forward_time += np.random.exponential(1/rate)
        
        # do not branch if age is already reached
        if forward_time >= age:
            break
        
        i = np.random.randint(len(branches))
        branch_id, parent = branches[i]
        spec_node = TreeNode(label=branch_id, event='S',
                             dist=forward_time-parent.tstamp,
                             tstamp=forward_time)
        parent.add_child(spec_node)
        branches[i] = (node_counter, spec_node)
        branches.append((node_counter+1, spec_node))
        node_counter += 2
    
    # finalize the branches
    for branch_id, parent in branches:
        parent.add_child( TreeNode(label=branch_id, event='S',
                                   dist=age-parent.tstamp,
                                   tstamp=age) )
    
    # reverse such that t(root) = age and t(leaves) = 0.0
    _reverse_time_stamps(tree)
    
    return tree


def _BDP(N, **kwargs):
    
    # remove potentially supplied 'episodes' argument
    episodes = _EBDP_check_episodes(birth_rate = kwargs.get('birth_rate'),
                                    death_rate = kwargs.get('death_rate'))
    
    return _EBDP_backward(N, episodes)


def _EBDP(N, **kwargs):
    
    episodes = _EBDP_check_episodes(**kwargs)
    
    return _EBDP_backward(N, episodes)


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


def _EBDP_backward(N, episodes, max_tries=500):
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
        
        branches = [TreeNode(label=j, event='S', dist=0.0, tstamp=t)
                    for j in range(N)]
        id_counter = N
        
        while branches:
            birth_i, death_i, rho_i, t_i = episodes[i]
            
            losses_to_add = round(len(branches) / rho_i) - len(branches)
            for j in range(losses_to_add):
                branches.append( TreeNode(label=id_counter,
                                          event='L',
                                          dist=0.0, tstamp=t) )
            id_counter += losses_to_add
            
            while branches:
                w = np.random.exponential( 1 / ((birth_i + death_i) * len(branches)) )
                
                if i+1 < len(episodes) and t + w > episodes[i+1][3]:
                    t = episodes[i+1][3]
                    i += 1
                    break
                
                else:
                    t += w
                    
                    if birth_i > np.random.uniform(low=0.0, high=birth_i+death_i):
                        # speciation event drawn
                        spec_node = TreeNode(label=id_counter, event='S',
                                             dist=0.0, tstamp=t)
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
                            branches.clear()
                    else:
                        # extinction event drawn
                        branches.append( TreeNode(label=id_counter,
                                                  event='L',
                                                  dist=0.0, tstamp=t) )
                        id_counter += 1
        
        # return tree with the following probability
        if np.random.random() < (1 / birth_i) / birth_inv_sum:
            
            for v in tree.preorder():
                if v.parent:
                    v.dist = v.parent.tstamp - v.tstamp
            
            return tree
        
    print("Could not return a tree after {} simulations!".format(max_tries),
          file=sys.stderr)
   

def _BDP_age(age, **kwargs):
    
    # remove potentially supplied 'episodes' argument
    episodes = _EBDP_age_check_episodes(birth_rate = kwargs.get('birth_rate'),
                                        death_rate = kwargs.get('death_rate'))
    
    return _EBDP_age_forward(age, episodes)


def _EBDP_age(age, **kwargs):
    
    episodes = _EBDP_age_check_episodes(**kwargs)
    
    return _EBDP_age_forward(age, episodes)


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

   
def _EBDP_mass_extinction(branches, surviving_rate, t):
    
    no_of_losses = round((1-surviving_rate) * len(branches))
    chosen_losses = sorted(np.random.choice(len(branches), replace=False,
                                            size=no_of_losses),
                           reverse=True)
    for j in chosen_losses:
        branch_id, parent = branches[j]
        loss_node = TreeNode(label=branch_id,
                             event='L',
                             dist=t-parent.tstamp,
                             tstamp=t)
        parent.add_child(loss_node)
        branches.pop(j)
    

def _EBDP_age_forward(age, episodes):
    """Episodic birth–death process (EBDP), forward algorithm with max. age."""
    
    tree = Tree(TreeNode(label=0, event='S', dist=0.0, tstamp=0.0))
    
    branches = [(1, tree.root)]
    forward_time = 0.0
    node_counter = 1
    i = 0               # current episode
    
    # may lead to extinction of the single branch at time t=0
    _EBDP_mass_extinction(branches, episodes[i][2], episodes[i][3])
    
    while forward_time < age:
        birth_rate, death_rate, *_ = episodes[i]
        
        rate = len(branches) * (birth_rate + death_rate)
        waiting_time = np.random.exponential(1/rate) if rate > 0.0 else float('inf')
        
        if i+1 < len(episodes) and forward_time + waiting_time >= episodes[i+1][3]:
            _EBDP_mass_extinction(branches, episodes[i+1][2], episodes[i+1][3])
            forward_time = episodes[i+1][3]
            i += 1
        
        elif forward_time + waiting_time >= age:
            break
        
        else:
            forward_time += waiting_time
            
            j = np.random.randint(len(branches))
            branch_id, parent = branches[j]
            
            if birth_rate > np.random.uniform(low=0.0, high=birth_rate+death_rate):
                # speciation event drawn
                spec_node = TreeNode(label=branch_id, event='S',
                                     dist=forward_time-parent.tstamp,
                                     tstamp=forward_time)
                parent.add_child(spec_node)
                branches[j] = (node_counter, spec_node)
                branches.append((node_counter+1, spec_node))
                node_counter += 2
            else:
                # extinction event drawn
                loss_node = TreeNode(label=branch_id,
                                     event='L',
                                     dist=forward_time-parent.tstamp,
                                     tstamp=forward_time)
                parent.add_child(loss_node)
                branches.pop(j)
                
    # finalize the (surviving) branches
    for branch_id, parent in branches:
        parent.add_child( TreeNode(label=branch_id, event='S',
                                   dist=age-parent.tstamp,
                                   tstamp=age) )
    
    # reverse such that t(root) = age and t(surviving leaves) = 0.0
    for v in tree.preorder():
        v.tstamp = age - v.tstamp
    
    return tree