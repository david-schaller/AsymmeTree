# -*- coding: utf-8 -*-

"""
Tree Simulator.

Simulate species and gene trees.

Methods and classes in this module:
    - simulate_species_tree
    - make_ultrametric
    - GeneTreeSimulator
    - observable_tree
"""

import random
from collections import deque

import numpy as np

from asymmetree.tools.PhyloTree import PhyloTree, PhyloTreeNode
import asymmetree.simulator._SpeciesTreeModels as stm


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2020, David Schaller"


# --------------------------------------------------------------------------
#                    CONSTRUCTION OF THE SPECIES TREE
#
#                       with the innovations model
# --------------------------------------------------------------------------

def simulate_species_tree(N, planted=True, model='innovations',
                          non_binary=0.0):
    """Simulates a species tree S with N leaves.
    
    Keyword parameters:
        planted -- add a planted root that has the canonical root as its
                   single neighbor; default is True
        model -- simulation model to be applied; default is 'innovations'
        non_binary -- probability that an inner edge is contracted;
                      results in non-binary tree; default is 0.0
    """
    
    if isinstance(model, str) and model.lower() in ('innovation', 'innovations'):
        tree = stm._innovations_model(N, planted=planted)
    else:
        raise ValueError("Model '{}' is not available!".format(model))
        
    if non_binary > 0.0:
         edges = _select_edges_for_contraction(tree,
                                               min(non_binary, 1.0),
                                               exclude_planted_edge=True)
         tree.contract(edges)
         
    make_ultrametric(tree)
        
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
    
    __slots__ = ('ID', 'array_id', 'rate',
                 'parent', 'S_u', 'S_v', 'transferred')
    
    def __init__(self, ID, array_id, rate,
                 parent, S_u, S_v, transferred):
        
        self.ID = ID                                    # unique branch id
        self.array_id = array_id                        # index in rate array
        self.rate = rate                                # total event rate
        
        self.parent = parent                            # parent node
        self.S_u = S_u                                  # edge (S_u, S_v) into which
        self.S_v = S_v                                  # the branch is embedded
        self.transferred = transferred                  # 1 if HGT edge, 0 otherwise
        

class GeneTreeSimulator:
    
    def __init__(self, S):
        
        self.S = S                                      # species tree
        self.sorted_speciations = S.sorted_nodes()      # list of speciations sorted by time stamp
        self.sorted_edges = S.sorted_edges()            # edges of the species tree
                                                        # sort (u,v) by tstamp of u
    
    
    def simulate(self, DLH_rates,
                 dupl_polytomy=0.0,
                 prohibit_extinction='per_species'):
        """Simulate a gene tree along the specified species tree.
        
        Simulate a gene tree along the specified species tree where the parameter
        'DLH_rates' contains the rates for duplications, losses and HGTs.
        
        Keyword arguments:
            dupl_polytomy - allows non-binary duplication events by specifying
                the lambda parameter for a poisson distribution (copy number =
                drawn number + 2); default is 0.0
            prohibit_extinction - avoid the extinction of all members in any
                species ('per_species'), of the complete gene family
                ('per_family'), or no constraints (False); default is 
                'per_species'.
        """
        
        self.DLH_rates = DLH_rates
        self.rate_sum = sum(DLH_rates)
        self.d = DLH_rates[0]
        self.l = DLH_rates[1]
        self.h = DLH_rates[2]
        
        self._prohibit_extinction = prohibit_extinction
        
        self._dupl_polytomy = dupl_polytomy
        
        self._reset()
        
        return self._gillespie_simulation()
    
    
    def _reset(self):
        
        self.spec_queue = deque(self.sorted_speciations)    # queue for speciation events
        self.id_counter = 0
        
        self.total_surviving = 0                            # counter for surviving genes
        self.total_rate = 0                                 # total event rate (all branches)
        self.branches = []
        
        self.ES_to_b = {e: [] for e in self.sorted_edges}   # maps E(S) --> existing branches

    
    def _get_tstamp(self, t):
        
        if self._prohibit_extinction == 'per_family' and self.total_surviving == 1:
            rate = self.d + self.l
        else:
            rate = self.total_rate
        
        
        if rate <= 0.0:
            return -1
        else:
            return t - np.random.exponential(1/rate)
    
    
    def _get_branch_and_type(self):
        
        if self._prohibit_extinction == 'per_family' and self.total_surviving == 1:
            
            branch = None
            for b in self.branches:
                if b.rate > 0.0:
                    branch = b
                    break
            
            if np.random.uniform(high=self.d+self.h) <= self.d:
                event_type = 'D'
            else:
                event_type = 'H'
        
        else:
            r = np.random.uniform(high=self.total_rate)
            current_sum = 0.0
            
            for i in range(len(self.branches)):
                if r <= current_sum + self.branches[i].rate:
                    break
                current_sum += self.branches[i].rate
            
            branch = self.branches[i]
            
            loss_factor = 1
            if (self._prohibit_extinction == 'per_species' and
                len(self.ES_to_b[(branch.S_u,branch.S_v)]) <= 1):
                loss_factor = 0
            
            if r <= current_sum + self.d:
                event_type = 'D'
            elif r <= current_sum + self.d + loss_factor * self.l:
                event_type = 'L'
            else:
                event_type = 'H'
        
#        assert np.isclose(self.total_rate, sum([b.rate for b in self.branches])), "Sum of rates and total rate are not equal."
        
        return branch, event_type
    
    
    def _get_copy_number(self):
        
        if self._dupl_polytomy <= 0.0:
            return 2
        else:
            return 2 + np.random.poisson(lam=self._dupl_polytomy)


    def _coexisting_edges(self, tstamp, exclude_edge=None):
        """Return list of edges for the given timestamp."""
        
        valid_edges = []
        
        for edge in self.sorted_edges:
            
            if edge[0].tstamp <= tstamp:
                break
            elif edge[1].tstamp < tstamp and edge != exclude_edge:
                valid_edges.append(edge)
                
        return valid_edges


    def _gillespie_simulation(self):
        
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
                if event_type == 'D':
                    self._duplication(event_tstamp, branch)
                
                # loss
                elif event_type == 'L':
                    self._loss(event_tstamp, branch)
                    
                # HGT    
                elif event_type == 'H':
                    self._hgt(event_tstamp, branch)
                
                t = event_tstamp

        return T
    
    
    def _initiatialize_tree(self):
        
        if len(self.S.root.children) > 1:
            # root is a speciation event
            root = PhyloTreeNode(self.id_counter, label='S', color=self.S.root.ID, 
                                 dist=0.0, tstamp=self.S.root.tstamp)
        else:                    
            # planted species tree
            root = PhyloTreeNode(self.id_counter, color=self.S.root.ID,
                                 dist=0.0, tstamp=self.S.root.tstamp)
            
        T = PhyloTree(root)
        self.id_counter += 1
        self.spec_queue.popleft()
        
        if self._prohibit_extinction == 'per_species':
            rate = self.d + self.h
        else:
            rate = self.rate_sum
                
        for S_v in self.S.root.children:
            
            array_id = len(self.branches)
            new_branch = Branch(self.id_counter, array_id, rate,
                                T.root, self.S.root, S_v, 0)
            self.ES_to_b[(self.S.root, S_v)].append(new_branch)
            self.branches.append(new_branch)
            self.id_counter += 1
            
        self.total_rate += len(self.S.root.children) * rate   
        self.total_surviving = len(self.S.root.children)
            
        return T
            
    
    def _speciation(self):
        
        S_v = self.spec_queue.popleft()
        S_u = S_v.parent
        
        for branch in self.ES_to_b[(S_u, S_v)]:
            
            spec_node = PhyloTreeNode(branch.ID, label='S',
                                      color=S_v.ID, tstamp=S_v.tstamp,
                                      dist=abs(S_v.tstamp-branch.parent.tstamp),
                                      transferred=branch.transferred)
            branch.parent.add_child(spec_node)
            
            if not S_v.children:
                spec_node.label = str(spec_node.ID)
                
            for S_w in S_v.children:
                
                if S_w is S_v.children[0]:
                    array_id = branch.array_id
                    new_branch = Branch(self.id_counter, array_id, branch.rate,
                                        spec_node, S_v, S_w, 0)
                    self.branches[array_id] = new_branch
                else:
                    array_id = len(self.branches)
                    new_branch = Branch(self.id_counter, array_id, branch.rate,
                                        spec_node, S_v, S_w, 0)
                    self.branches.append(new_branch)
                
                self.ES_to_b[(S_v, S_w)].append(new_branch)
                self.id_counter += 1
                
            self.total_rate += (len(S_v.children) - 1) * branch.rate
            self.total_surviving += len(S_v.children) - 1
            
    
    def _duplication(self, event_tstamp, branch):
        
        S_u, S_v = branch.S_u, branch.S_v
        
        dupl_node = PhyloTreeNode(branch.ID, label='D',
                                  color=(S_u.ID,S_v.ID),
                                  tstamp=event_tstamp,
                                  dist=abs(event_tstamp-branch.parent.tstamp),
                                  transferred=branch.transferred)
        branch.parent.add_child(dupl_node)
        self.ES_to_b[(S_u, S_v)].remove(branch)
        
        copy_number = self._get_copy_number()
        
        for i in range(copy_number):
            
            if i == 0:
                array_id = branch.array_id
                new_branch = Branch(self.id_counter, array_id, self.rate_sum,
                                    dupl_node, S_u, S_v, 0)
                self.branches[array_id] = new_branch
            else:
                array_id = len(self.branches)
                new_branch = Branch(self.id_counter, array_id, self.rate_sum,
                                    dupl_node, S_u, S_v, 0)
                self.branches.append(new_branch)
            
            self.ES_to_b[(S_u, S_v)].append(new_branch)
            self.id_counter += 1
        
        self.total_surviving += (copy_number - 1)
        self.total_rate += (copy_number - 1) * (self.rate_sum)
        
        if (self._prohibit_extinction == 'per_species' and
            len(self.ES_to_b[(S_u, S_v)]) == copy_number):
            self.total_rate += self.l
            
            
    def _loss(self, event_tstamp, branch):
        
        S_u, S_v = branch.S_u, branch.S_v
        
        loss_node = PhyloTreeNode(branch.ID, label='*',
                                  color=(S_u.ID,S_v.ID),
                                  tstamp=event_tstamp,
                                  dist=abs(event_tstamp-branch.parent.tstamp),
                                  transferred=branch.transferred)
        branch.parent.add_child(loss_node)
        self.ES_to_b[(S_u, S_v)].remove(branch)
        branch.rate = 0.0
        
        self.total_surviving -= 1
        self.total_rate -= self.rate_sum
        
        if (self._prohibit_extinction == 'per_species' and
            len(self.ES_to_b[(S_u, S_v)]) == 1):
            self.total_rate -= self.l
            self.ES_to_b[(S_u, S_v)][0].rate -= self.l
            
    
    
    def _hgt(self, event_tstamp, branch):
        
        S_u, S_v = branch.S_u, branch.S_v
        
        valid_edges = self._coexisting_edges(event_tstamp,
                                             exclude_edge=(S_u, S_v))
        if valid_edges:
            trans_edge = random.choice(valid_edges)
            hgt_node = PhyloTreeNode(branch.ID, label='H',
                                     color=(S_u.ID,S_v.ID),
                                     tstamp=event_tstamp,
                                     dist=abs(event_tstamp-branch.parent.tstamp),
                                     transferred=branch.transferred)
            branch.parent.add_child(hgt_node)
            self.ES_to_b[(S_u, S_v)].remove(branch)
            
            # original branch
            array_id = branch.array_id
            new_branch = Branch(self.id_counter, array_id, branch.rate,
                                hgt_node, S_u, S_v, 0) 
            self.branches[array_id] = new_branch
            self.ES_to_b[(S_u, S_v)].append(new_branch)
            self.id_counter += 1
            
            # receiving branch
            array_id = len(self.branches)
            trans_branch = Branch(self.id_counter, array_id, self.rate_sum,
                                  hgt_node, *trans_edge, 1)
            self.branches.append(trans_branch)
            self.ES_to_b[trans_edge].append(trans_branch)
            self.id_counter += 1
            
            self.total_surviving += 1
            self.total_rate += self.rate_sum
            
            if (self._prohibit_extinction == 'per_species' and
                len(self.ES_to_b[trans_edge]) == 2):
                self.total_rate += self.l
                self.ES_to_b[trans_edge][0].rate += self.l 
                
    
    def _assert_no_extinction(self, T):
        """Returns False if gene family is extinct in some species."""
        
        VS_to_VT = {l.ID: [] for l in self.S.preorder() if not l.children}
        
        for v in T.preorder():
            if not v.children and not v.label == '*':
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
        if not node.children and node.label == '*':
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
        
        new_root = obs_tree.root.children[0]
        new_root.detach()
        obs_tree.root = new_root
        new_root.dist = 0.0
        
    elif not obs_tree.root.children and not obs_tree.root.label:
        # no surviving genes --> return empty tree
        obs_tree.root = None
    
    return obs_tree