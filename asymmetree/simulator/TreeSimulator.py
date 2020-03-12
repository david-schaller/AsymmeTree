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
        non_binary -- probability that an inner edge is contracted;
                      results in non-binary tree; default is 0.0
    """
    
    if isinstance(model, str) and model.lower() in ('innovation', 'innovations'):
        tree = _innovations_model(N, planted=planted)
    else:
        raise ValueError("Model '{}' is not available!".format(model))
        
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


class Branch:
    
    def __init__(self, ID, parent, S_u, S_v, transferred):
        
        self.ID = ID
        self.parent = parent
        self.S_u = S_u
        self.S_v = S_v
        self.transferred = transferred
        

class GeneTreeSimulator:
    
    def __init__(self, S):
        
        self.S = S                                      # species tree
        self.sorted_speciations = S.sorted_nodes()      # list of speciations sorted by time stamp
        self.sorted_edges = S.sorted_edges()            # edges of the species tree
                                                        # sort (u,v) by tstamp of u
    
    
    def simulate(self, DLH_rates,
                 dupl_polytomy=0.0):
        """Simulate a gene tree along the specified species tree.
        
        Simulate a gene tree along the specified species tree where the parameter
        'DLH_rates' contains the rates for duplications, losses and HGTs.
        
        Keyword arguments:
            dupl_polytomy - allows non-binary duplication events by specifying
                            the lambda parameter for a poisson distribution
                            (copy number = drawn number + 2); default is 0.0
        """
        
        self.DLH_rates = DLH_rates
        self.rate_sum = sum(DLH_rates)
        self.d = DLH_rates[0]
        self.l = DLH_rates[1]
        self.h = DLH_rates[2]
        
        self.dupl_polytomy = dupl_polytomy
        
        self._reset()
        
        return self._gillespie_direct()
    
    
    def _reset(self):
        
        self.spec_queue = deque(self.sorted_speciations)    # queue for speciation events
        self.id_counter = 0
        
        self.total_rate = 0                                 # total event rate (all branches)
        self.branch_rates = []                              # array for branch rates
        
        self.ES_to_b = {e: [] for e in self.sorted_edges}   # maps E(S) --> existing branches
        self.i_to_b = {}                                    # maps array index --> branch
        self.b_to_i = {}                                    # maps branch --> array index

    
    def _get_tstamp(self, t):
        
        if self.total_rate <= 0.0:
            return -1
        else:
            return t - np.random.exponential(1/self.total_rate)
    
    
    def _get_branch_and_type(self):
        
        r = np.random.uniform() * self.total_rate
        index, current_sum = 0, 0
        
        for rate in self.branch_rates:
            if r <= current_sum + rate:
                break
            current_sum += rate
            index += 1
        
        branch = self.i_to_b[index]
        loss_factor = 1 if len(self.ES_to_b[(branch.S_u,branch.S_v)]) > 1 else 0
        
        if r <= current_sum + self.DLH_rates[0]:
            event_type = "D"
        elif r <= current_sum + self.DLH_rates[0] + loss_factor * self.DLH_rates[1]:
            event_type = "L"
        else:
            event_type = "H"
        
#        assert np.isclose(self.total_rate, sum(self.branch_rates)), "Sum of rates and total rate are not equal."
        
        return branch, event_type
    
    
    def _get_copy_number(self):
        
        if self.dupl_polytomy <= 0.0:
            return 2
        else:
            return 2 + np.random.poisson(lam=self.dupl_polytomy)


    def _coexisting_edges(self, tstamp, exclude_edge=None):
        """Return list of edges for the given timestamp."""
        
        valid_edges = []
        for edge in self.sorted_edges:
            if edge[0].tstamp <= tstamp:
                break
            elif edge[1].tstamp < tstamp and edge != exclude_edge:
                valid_edges.append(edge)
        return valid_edges


    def _gillespie_direct(self):
        
        T = self._initiatialize_tree()
        t = T.root.tstamp                   # start time = 1.0
        
        while self.spec_queue:
            
            event_tstamp = self._get_tstamp(t)
            next_spec_tstamp = self.spec_queue[0].tstamp
            
            # speciation
            if event_tstamp <= next_spec_tstamp:
                self._speciation()      
                t = next_spec_tstamp
                
            else:
                branch, event_type = self._get_branch_and_type()
                
                # duplication
                if event_type == "D":
                    self._duplication(event_tstamp, branch, event_type)
                
                # loss
                elif event_type == "L":
                    self._loss(event_tstamp, branch, event_type)
                    
                # HGT    
                elif event_type == "H":
                    self._hgt(event_tstamp, branch, event_type)
                
                t = event_tstamp

        return T
    
    
    def _initiatialize_tree(self):
        
        if len(self.S.root.children) > 1:
            # root is a speciation event
            T = PhyloTree(PhyloTreeNode(self.id_counter, label="S",
                                        color=self.S.root.ID, 
                                        dist=0.0,
                                        tstamp=self.S.root.tstamp))
        else:                    
            # planted species tree
            T = PhyloTree(PhyloTreeNode(self.id_counter,
                                        color=self.S.root.ID,
                                        dist=0.0,
                                        tstamp=self.S.root.tstamp))
        self.id_counter += 1
        self.spec_queue.popleft()
                
        for S_v in self.S.root.children:
            new_branch = Branch(self.id_counter, T.root, self.S.root, S_v, 0)
            self.ES_to_b[(self.S.root, S_v)].append(new_branch)
            index = len(self.branch_rates)
            self.branch_rates.append(self.d + self.h)
            self.total_rate += self.d + self.h
            self.i_to_b[index] = new_branch
            self.b_to_i[new_branch] = index
            self.id_counter += 1
            
        return T
            
    
    def _speciation(self):
        
        S_v = self.spec_queue.popleft()
        S_u = S_v.parent
        
        for branch in self.ES_to_b[(S_u, S_v)]:
            spec_node = PhyloTreeNode(branch.ID, label="S",
                                 color=S_v.ID, tstamp=S_v.tstamp,
                                 dist=abs(S_v.tstamp-branch.parent.tstamp),
                                 transferred=branch.transferred)
            branch.parent.add_child(spec_node)
            if not S_v.children:
                spec_node.label = str(spec_node.ID)
            for S_w in S_v.children:
                new_branch = Branch(self.id_counter, spec_node, S_v, S_w, 0)
                self.ES_to_b[(S_v, S_w)].append(new_branch)
                if S_w is S_v.children[0]:
                    index = self.b_to_i[branch]
                else:
                    index = len(self.branch_rates)
                    self.branch_rates.append(self.branch_rates[self.b_to_i[branch]])
                    self.total_rate += self.branch_rates[index]
                self.i_to_b[index] = new_branch
                self.b_to_i[new_branch] = index
                self.id_counter += 1
            del self.b_to_i[branch]
            
    
    def _duplication(self, event_tstamp, branch, event_type):
        
        S_u, S_v = branch.S_u, branch.S_v
        
        dupl_node = PhyloTreeNode(branch.ID, label="D",
                                  color=(S_u.ID,S_v.ID),
                                  tstamp=event_tstamp,
                                  dist=abs(event_tstamp-branch.parent.tstamp),
                                  transferred=branch.transferred)
        branch.parent.add_child(dupl_node)
        self.ES_to_b[(S_u, S_v)].remove(branch)
        
        copy_number = self._get_copy_number()
        
        for i in range(copy_number):
            new_branch = Branch(self.id_counter, dupl_node, S_u, S_v, 0)
            self.ES_to_b[(S_u, S_v)].append(new_branch)
            if i == 0:
                index = self.b_to_i[branch]
                self.branch_rates[index] = self.rate_sum
            else:
                index = len(self.branch_rates)
                self.branch_rates.append(self.rate_sum)
            self.i_to_b[index] = new_branch
            self.b_to_i[new_branch] = index
            self.id_counter += 1
        del self.b_to_i[branch]
        
        self.total_rate += (copy_number - 1) * (self.rate_sum)
        if len(self.ES_to_b[(S_u, S_v)]) == copy_number:
            self.total_rate += self.l
            
            
    def _loss(self, event_tstamp, branch, event_type):
        
        S_u, S_v = branch.S_u, branch.S_v
        
        loss_node = PhyloTreeNode(branch.ID, label="*",
                                  color=(S_u.ID,S_v.ID),
                                  tstamp=event_tstamp,
                                  dist=abs(event_tstamp-branch.parent.tstamp),
                                  transferred=branch.transferred)
        branch.parent.add_child(loss_node)
        self.ES_to_b[(S_u, S_v)].remove(branch)
        self.branch_rates[self.b_to_i[branch]] = 0
        del self.b_to_i[branch]
        
        if len(self.ES_to_b[(S_u, S_v)]) == 1:
            self.total_rate -= (self.rate_sum + self.l)
            self.branch_rates[self.b_to_i[self.ES_to_b[(S_u, S_v)][0]]] -= self.l
        else:
            self.total_rate -= self.rate_sum
    
    
    def _hgt(self, event_tstamp, branch, event_type):
        
        S_u, S_v = branch.S_u, branch.S_v
        
        valid_edges = self._coexisting_edges(event_tstamp,
                                             exclude_edge=(S_u, S_v))
        if valid_edges:
            trans_edge = random.choice(valid_edges)
            hgt_node = PhyloTreeNode(branch.ID, label="H",
                                     color=(S_u.ID,S_v.ID),
                                     tstamp=event_tstamp,
                                     dist=abs(event_tstamp-branch.parent.tstamp),
                                     transferred=branch.transferred)
            branch.parent.add_child(hgt_node)
            self.ES_to_b[(S_u, S_v)].remove(branch)
            
            # original branch
            new_branch = Branch(self.id_counter, hgt_node, S_u, S_v, 0)        
            self.ES_to_b[(S_u, S_v)].append(new_branch)
            index = self.b_to_i[branch]
            self.i_to_b[index] = new_branch
            self.b_to_i[new_branch] = index
            self.id_counter += 1
            
            # receiving branch
            trans_branch = Branch(self.id_counter, hgt_node, *trans_edge, 1)
            self.ES_to_b[trans_edge].append(trans_branch)
            index = len(self.branch_rates)
            self.branch_rates.append(self.rate_sum)
            self.i_to_b[index] = trans_branch
            self.b_to_i[trans_branch] = index
            self.id_counter += 1
            
            del self.b_to_i[branch]
            
            if len(self.ES_to_b[trans_edge]) == 2:
                self.total_rate += (self.rate_sum + self.l)
                self.branch_rates[self.b_to_i[self.ES_to_b[trans_edge][0]]] += self.l
            else:
                self.total_rate += self.rate_sum
                
    
    def _assert_no_extinction(self, T):
        """Returns False if gene family is extinct in some species."""
        
        VS_to_VT = {l.ID: [] for l in self.S.preorder() if not l.children}
        
        for v in T.preorder():
            if not v.children and not v.label == "*":
                VS_to_VT[v.color].append(v.ID)
                
        for leaf_list in VS_to_VT.values():
            if not leaf_list:
                return False
        
        return True
        

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