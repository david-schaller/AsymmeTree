# -*- coding: utf-8 -*-

"""
Tree Simulator.

Simulate species and gene trees.

Methods in this module:
    - build_species_tree
    - make_ultrametric
    - build_gene_tree
    - observable_tree
"""

import random
from collections import deque

import numpy as np

from asymmetree.tools.PhyloTree import PhyloTree, PhyloTreeNode


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2020, David Schaller"


# --------------------------------------------------------------------------
#                    CONSTRUCTION OF THE SPECIES TREE
#
#                       with the innovations model
# --------------------------------------------------------------------------

def build_species_tree(N, planted=True, model="innovations",
                       non_binary=0.0):
    """Builds a species tree S with N leaves.
    
    Keyword parameters:
        planted -- add a planted root that has the canonical root as its
                   single neighbor; default is True
        model -- simulation model to be applied; default is 'innovations'
        non_binary -- probality that an inner edge is contracted;
                      results in non-binary tree; default is 0.0
    """
    
    if isinstance(model, str) and model.lower() in ('innovation', 'innovations'):
        tree = _innovations_model(N, planted=planted)
    else:
        raise ValueError(f"Model '{model}' is not available!")
        
    if non_binary > 0.0:
         edges = _select_edges_for_contraction(tree,
                                               min(non_binary, 1.0),
                                               exclude_planted_edge=True)
         tree.contract(edges)
         
    make_ultrametric(tree)
        
    return tree
        

def _innovations_model(N, planted=True):
    """Builds a species tree S with N leaves with the innovations model."""
    
    tree = PhyloTree(PhyloTreeNode(0, label="0"))
    tree.number_of_species = N
    
    if planted:                     # planted tree (root is an implicit
                                    # outgroup with outdegree = 1)
        root = PhyloTreeNode(1, label="1")
        tree.root.add_child(root)
        node_counter = 2
    else:
        root = tree.root
        node_counter = 1
    
    t = 1
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
            t += 1
            
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
                t += 1
    
    return tree


def make_ultrametric(tree):
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
            
            
def _select_edges_for_contraction(tree, p, exclude_planted_edge=True):
    
    edges = []
    
    for u, v in tree.inner_edges():
        
        if exclude_planted_edge and (u is tree.root) and len(u.children) == 1:
            continue
        
        if random.random() < p:
            edges.append((u,v))
    
    return edges


# --------------------------------------------------------------------------
#                    CONSTRUCTION OF THE GENE TREE
#
# --------------------------------------------------------------------------
            
def build_gene_tree(S, DLH_rates):
    
    return _gillespie_direct(S, DLH_rates)

    
def _get_tstamp(t, total_rate):
    
    if total_rate <= 0.0:
        return -1
    else:
        return t - np.random.exponential(1/total_rate)


def _get_branch_and_type(DLH_rates, total_rate, branch_rates, i_to_b, ES_to_b):
    
    r = np.random.uniform() * total_rate
    index, current_sum = 0, 0
    
    for rate in branch_rates:
        if r <= current_sum + rate:
            break
        current_sum += rate
        index += 1
    
    branch = i_to_b[index]
    loss_factor = 1 if len(ES_to_b[(branch[2],branch[3])]) > 1 else 0
    
    if r <= current_sum + DLH_rates[0]:
        event_type = "D"
    elif r <= current_sum + DLH_rates[0] + loss_factor * DLH_rates[1]:
        event_type = "L"
    else:
        event_type = "H"
    
#    rate_sum = sum(branch_rates)
#    if total_rate != rate_sum:
#        raise KeyboardInterrupt
    
    return branch, event_type


def _coexisting_edges(sorted_edges, tstamp, exclude_edge=None):
    """Return list of edges for the given timestamp."""
    
    valid_edges = []
    for edge in sorted_edges:
        if edge[0].tstamp <= tstamp:
            break
        elif edge[1].tstamp < tstamp and edge != exclude_edge:
            valid_edges.append(edge)
    return valid_edges


def _gillespie_direct(S, DLH_rates):
    
    # --------------- initialization -----------------
    d, l, h = DLH_rates
    S_edges = S.sorted_edges()                      # edges of the species tree
                                                    # sort (u,v) by tstamp of u
    speciations = deque(S.sorted_nodes())           # queue for speciation events
    
    total_rate = 0                                  # total event rate (all branches)
    branch_rates = []                               # array for branch rates
    
    ES_to_b = {e: [] for e in S_edges}              # maps E(S) --> existing branches
    i_to_b = {}                                     # maps array index --> branch
    b_to_i = {}                                     # maps branch --> array index
    
    if len(S.root.children) > 1:
        T = PhyloTree(PhyloTreeNode(0, label="S",   # root is a speciation event
                                    color=S.root.ID, 
                                    dist=0.0, tstamp=1.0))
    else:                                           # planted tree
        T = PhyloTree(PhyloTreeNode(0, color=S.root.ID,
                                    dist=0.0, tstamp=1.0))
    speciations.popleft()
    
    id_counter = 1
    
    for S_v in S.root.children:
        new_branch = (id_counter, T.root, S.root, S_v, 0)
        ES_to_b[(S.root, S_v)].append(new_branch)
        index = len(branch_rates)
        branch_rates.append(d+h)
        total_rate += d+h
        i_to_b[index] = new_branch
        b_to_i[new_branch] = index
        id_counter += 1
        
    t = T.root.tstamp                               # start time = 1.0
    
    while speciations:
        
        event_tstamp = _get_tstamp(t, total_rate)
        
        # ----------------- SPECIATION -----------------
        if event_tstamp <= speciations[0].tstamp:
            S_v = speciations.popleft()
            S_u = S_v.parent
            
            for branch in ES_to_b[(S_u, S_v)]:
                b_id, b_parent, _, _, b_transferred = branch
                spec_node = PhyloTreeNode(b_id, label="S",
                                     color=S_v.ID, tstamp=S_v.tstamp,
                                     dist=abs(S_v.tstamp-b_parent.tstamp),
                                     transferred=b_transferred)
                b_parent.add_child(spec_node)
                if not S_v.children:
                    spec_node.label = str(spec_node.ID)
                for S_w in S_v.children:
                    new_branch = (id_counter, spec_node, S_v, S_w, 0)
                    ES_to_b[(S_v, S_w)].append(new_branch)
                    if S_w is S_v.children[0]:
                        index = b_to_i[branch]
                    else:
                        index = len(branch_rates)
                        branch_rates.append(branch_rates[b_to_i[branch]])
                        total_rate += branch_rates[-1]
                    i_to_b[index] = new_branch
                    b_to_i[new_branch] = index
                    id_counter += 1
                del b_to_i[branch]
                    
                
            t = S_v.tstamp
            
        # ----------------- D / L / H ------------------
        else:
            branch, event_type = _get_branch_and_type(DLH_rates, total_rate,
                                                      branch_rates, i_to_b, ES_to_b)
            b_id, b_parent, S_u, S_v, b_transferred = branch
            
            # ----------------- DUPLICATION -----------------
            if event_type == "D":
                dupl_node = PhyloTreeNode(b_id, label="D",
                                          color=(S_u.ID,S_v.ID),
                                          tstamp=event_tstamp,
                                          dist=abs(event_tstamp-b_parent.tstamp),
                                          transferred=b_transferred)
                b_parent.add_child(dupl_node)
                ES_to_b[(S_u, S_v)].remove(branch)
                for i in range(2):
                    new_branch = (id_counter, dupl_node, S_u, S_v, 0)
                    ES_to_b[(S_u, S_v)].append(new_branch)
                    if i == 0:
                        index = b_to_i[branch]
                        branch_rates[index] = d + l + h
                    else:
                        index = len(branch_rates)
                        branch_rates.append(d+l+h)
                    i_to_b[index] = new_branch
                    b_to_i[new_branch] = index
                    id_counter += 1
                del b_to_i[branch]
                
                if len(ES_to_b[(S_u, S_v)]) == 2:
                    total_rate += (d + l + h + l)
                else:
                    total_rate += (d + l + h)
            
            # -------------------- LOSS ---------------------
            elif event_type == "L":
                loss_node = PhyloTreeNode(b_id, label="*",
                                          color=(S_u.ID,S_v.ID),
                                          tstamp=event_tstamp,
                                          dist=abs(event_tstamp-b_parent.tstamp),
                                          transferred=b_transferred)
                b_parent.add_child(loss_node)
                ES_to_b[(S_u, S_v)].remove(branch)
                branch_rates[b_to_i[branch]] = 0
                del b_to_i[branch]
                
                if len(ES_to_b[(S_u, S_v)]) == 1:
                    total_rate -= (d + l + h + l)
                    branch_rates[b_to_i[ES_to_b[(S_u, S_v)][0]]] -= l
                else:
                    total_rate -= (d + l + h)
            
            # -------------------- HGT ----------------------
            elif event_type == "H":
                valid_edges = _coexisting_edges(S_edges, event_tstamp,
                                                exclude_edge=(S_u, S_v))
                if valid_edges:
                    trans_edge = random.choice(valid_edges)
                    hgt_node = PhyloTreeNode(b_id, label="H",
                                             color=(S_u.ID,S_v.ID),
                                             tstamp=event_tstamp,
                                             dist=abs(event_tstamp-b_parent.tstamp),
                                             transferred=b_transferred)
                    b_parent.add_child(hgt_node)
                    ES_to_b[(S_u, S_v)].remove(branch)
                    
                    new_branch = (id_counter, hgt_node, S_u, S_v, 0)        # original branch
                    ES_to_b[(S_u, S_v)].append(new_branch)
                    index = b_to_i[branch]
                    i_to_b[index] = new_branch
                    b_to_i[new_branch] = index
                    id_counter += 1
                    
                    trans_branch = (id_counter, hgt_node, *trans_edge, 1)   # receiving branch
                    ES_to_b[trans_edge].append(trans_branch)
                    index = len(branch_rates)
                    branch_rates.append(d+l+h)
                    i_to_b[index] = trans_branch
                    b_to_i[trans_branch] = index
                    id_counter += 1
                    
                    del b_to_i[branch]
                    
                    if len(ES_to_b[trans_edge]) == 2:
                        total_rate += (d + l + h + l)
                        branch_rates[b_to_i[ES_to_b[trans_edge][0]]] += l
                    else:
                        total_rate += (d + l + h)
            
            t = event_tstamp

#    # check that there is no extinction in any species
#    VS_to_VT = {l.ID: [] for l in S.preorder() if not l.children}
#    for v in T.preorder():
#        if not v.children and not v.label == "*":
#            VS_to_VT[v.color].append(v.ID)
#    for leaf_list in VS_to_VT.values():
#        if not leaf_list:
#            raise KeyboardInterrupt
    
    return T


# --------------------------------------------------------------------------
#                 CONSTRUCTION OF THE OBSERVABLE TREE
#
# --------------------------------------------------------------------------
    
def observable_tree(tree):
    obs_tree = tree.copy()
    
    loss_nodes = []
    for node in obs_tree.postorder():
        if not node.children and node.label == "*":
            loss_nodes.append(node)
    
    # traverse from loss node to root delete if degree <= 1
    for loss_node in loss_nodes:
        current = obs_tree.delete_and_reconnect(loss_node,
                                                add_distances=True,
                                                keep_transferred=True)
        
        while len(current.children) < 2 and current.parent:
            current = obs_tree.delete_and_reconnect(current,
                                                    add_distances=True,
                                                    keep_transferred=True)
    
    # delete the root if the tree is planted
    if len(obs_tree.root.children) == 1:
        obs_tree.delete_and_reconnect(obs_tree.root.children[0],
                                      add_distances=False,
                                      keep_transferred=False)
    
    return obs_tree