# -*- coding: utf-8 -*-

import random, sys
import numpy as np

from asymmetree.tools.PhyloTree import PhyloTree, PhyloTreeNode


def _innovations_model(N, planted, ultrametric=True):
    """Builds a species tree S with N leaves with the innovations model.
    
    Keyword arguments:
        ultrametric - it True make tree ultrametric and rescale it to
            depth 1.0, else all edges have length 1.0; default is True.
    """
    
    tree = PhyloTree(PhyloTreeNode(0, label='0'))
    tree.number_of_species = N
    
    if planted:                     # planted tree (root is an implicit
                                    # outgroup with outdegree = 1)
        root = PhyloTreeNode(1, label="1")
        tree.root.add_child(root)
        node_counter = 2
    else:
        root = tree.root
        node_counter = 1
    
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
            
            child1 = PhyloTreeNode(node_counter, label=str(node_counter))
            species[s].add_child(child1)
            child2 = PhyloTreeNode(node_counter+1, label=str(node_counter+1))
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
                
                child1 = PhyloTreeNode(node_counter, label=str(node_counter))
                species[s].add_child(child1)
                child2 = PhyloTreeNode(node_counter+1, label=str(node_counter+1))
                species[s].add_child(child2)

                node_counter += 2
                
                species[s] = child1
                species[new_s] = child2
                
    if ultrametric:         
        _make_ultrametric(tree)
    
    return tree


def _make_ultrametric(tree):
    """Makes a given species tree S ultrametric.
    
    It is t(root) = 1 and t(x) = 0 for x in L(S)."""
    
    for v in tree.preorder():
        if not v.parent:
            v.dist = 0.0
            v.tstamp = 1.0
        elif not v.children:
            v.dist = v.parent.tstamp
            v.tstamp = 0.0
        else:                               # random walk to a leaf
            pos = v                         # current position
            length = 0                      # path length |P|
            while pos.children:
                length += 1
                pos = pos.children[np.random.randint(len(pos.children))]
            v.dist = (v.parent.tstamp) * 2 * np.random.uniform() / (length+1)
            v.tstamp = v.parent.tstamp - v.dist


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
        raise ValueError("Birth rate must be >0!")
    
    tree = PhyloTree(PhyloTreeNode(0, label='0', dist=0.0, tstamp=0.0))
    tree.number_of_species = N
    
    branches = [(1, tree.root)]
    forward_time = 0.0
    node_counter = 1
    
    while len(branches) < N:
        
        rate = len(branches) * birth_rate
        forward_time += np.random.exponential(1/rate)
        
        i = np.random.randint(len(branches))
        branch_id, parent = branches[i]
        spec_node = PhyloTreeNode(branch_id, label=str(branch_id),
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
        parent.add_child( PhyloTreeNode(branch_id, label=str(branch_id),
                                        dist=forward_time-parent.tstamp,
                                        tstamp=forward_time) )
        
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
        
        for episode in episodes:
            if len(episode) != 4:
                raise ValueError("All episodes must contain 4 values: birth rate, "\
                                 "death rate, proportion of survivors, time stamp "\
                                 "(from recent time as 0)!")
            
            birth_rate, death_rate, rho, t = episode
            
            if birth_rate <= 0.0 or birth_rate < death_rate:
                raise ValueError("Birth rate must be >0 and >=death rate "\
                                 "in all episodes!")
                
            if rho <= 0.0 or rho > 1.0:
                raise ValueError("Proportion of survivors must be in (0.0, 1.0]!")
        
        return episodes
    
    elif birth_rate is not None:
        
        if birth_rate <= 0.0 or (death_rate and birth_rate < death_rate):
            raise ValueError("Birth rate must be >0 and >=death rate!")
            
        if death_rate and death_rate < 0.0:
            raise ValueError("Death rate must be >=0!")
        
        if death_rate is None:
            return [(birth_rate, 0.0, 1.0, 0.0)]
        else:
            return [(birth_rate, death_rate, 1.0, 0.0)]
        
    else:
        if death_rate:
            raise ValueError("Birth rate (>0) must be specified if death rate "\
                             "is supplied!")
            
        return [(birth_rate, 0.0, 1.0, 0.0)]


def _EBDP_backward(N, episodes, max_tries=500):
    """Episodic birth–death process (EBDP), backward algorithm by Stadler 2001."""
    
    birth_inv_sum = sum([1/episodes[i][0] for i in range(len(episodes))])
    
    for _ in range(max_tries):
        
        tree = None
        t = 0.0
        i = 0
        
        branches = [PhyloTreeNode(j, str(j), dist=0.0, tstamp=t)
                    for j in range(N)]
        id_counter = N
        
        while branches:
            birth_i, death_i, rho_i, t_i = episodes[i]
            
            losses_to_add = round(len(branches) / rho_i) - len(branches)
            for j in range(losses_to_add):
                branches.append( PhyloTreeNode(id_counter+j, '*', dist=0.0, tstamp=t) )
            id_counter += losses_to_add
            
            while branches:
                w = np.random.exponential( 1 / ((birth_i + death_i) * len(branches)) )
                
                if i+1 < len(episodes) and t + w < episodes[i+1][3]:
                    
                    i += 1
                    t = episodes[i+1][3]
                    break
                
                else:
                    
                    t += w
                    
                    if birth_i > np.random.uniform(low=0.0, high=birth_i+death_i):
                        # speciation event drawn
                        spec_node = PhyloTreeNode(id_counter, label=str(id_counter),
                                                  dist=0.0, tstamp=t)
                        id_counter += 1
                        if len(branches) > 1:
                            k, l = np.random.choice(len(branches), 2, replace=False)
                            if k > l:
                                k, l = l, k
                            spec_node.add_child(branches[k])
                            spec_node.add_child(branches[l])
                            branches[k] = spec_node
                            branches.pop(l)
                        else:
                            spec_node.add_child(branches[0])
                            tree = PhyloTree(spec_node)
                            branches.clear()
                    else:
                        # extinction event drawn
                        branches.append( PhyloTreeNode(id_counter, '*', dist=0.0, tstamp=t) )
                        id_counter += 1
        
        # return tree with the following probability
        if np.random.random() < (1 / birth_i) / birth_inv_sum:
            
            for v in tree.preorder():
                if v.parent:
                    v.dist = v.parent.tstamp - v.tstamp
            
            return tree
        
    print("Could not return a tree after {} simulations!".format(max_tries),
          file=sys.stderr)
    


if __name__ == '__main__':
    
    T = _yule(10, 1.0)
    print(T.to_newick())
    
    print('------------------------')
#    T2 = _EBDP(10, episodes=[(1.0, 0.3, 0.8, 0.0)])
    T2 = _EBDP(10, birth_rate=1, death_rate=0.5)
    print(T2.to_newick())
    