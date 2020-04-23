# -*- coding: utf-8 -*-

import random, sys
import numpy as np

from asymmetree.tools.PhyloTree import PhyloTree, PhyloTreeNode


def _innovations_model(N, planted=True):
    """Builds a species tree S with N leaves with the innovations model."""
    
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
    
    return tree


def _reverse_time_stamps(tree):
    
    max_depth = 0.0
    for v in tree.preorder():
        max_depth = max(max_depth, v.tstamp)
            
    for v in tree.preorder():
        v.tstamp = abs(v.tstamp - max_depth)   


def _yule(N, birth_rate):
    
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


def _EBDP_backward(N, parameters, max_tries=500):
    
    birth_inv_sum = sum([1/parameters[i][0] for i in range(len(parameters))])
    
    for _ in range(max_tries):
        
        tree = None
        t = 0.0
        i = 0
        
        branches = [PhyloTreeNode(j, str(j), dist=0.0, tstamp=t)
                    for j in range(N)]
        id_counter = N
        
        while branches:
            birth_i, death_i, rho_i, t_i = parameters[i]
            
            losses_to_add = round(len(branches) / rho_i) - len(branches)
            for j in range(losses_to_add):
                branches.append( PhyloTreeNode(id_counter+j, '*', dist=0.0, tstamp=t) )
            id_counter += losses_to_add
            
            while branches:
                w = np.random.exponential( 1 / ((birth_i + death_i) * len(branches)) )
                
                if i+1 < len(parameters) and t + w < parameters[i+1][3]:
                    
                    i += 1
                    t = parameters[i+1][3]
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
    T2 = _EBDP_backward(10, [(1.0, 0.3, 0.8, 0.0)])
    print(T2.to_newick())
    