# -*- coding: utf-8 -*-

"""
Gene Tree Simulator.

Simulate dated gene trees.
"""

import random, warnings
from collections import deque

import numpy as np

from tralda.datastructures.Tree import Tree, TreeNode

from asymmetree.tools.PhyloTreeTools import (sorted_nodes,
                                             sorted_edges,
                                             delete_losses_and_contract,
                                             remove_planted_root)


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                    CONSTRUCTION OF THE GENE TREE
# --------------------------------------------------------------------------
    

def simulate_dated_gene_tree(S, **kwargs):
    
    gene_tree_simulator = GeneTreeSimulator(S)
    return gene_tree_simulator.simulate(**kwargs)


class _Branch:
    
    __slots__ = ('label', 'array_id', 'rate',
                 'parent', 'S_u', 'S_v', 'transferred')
    
    def __init__(self, label, array_id, rate,
                 parent, S_u, S_v, transferred):
        
        self.label = label                          # unique branch label
        self.array_id = array_id                    # index in rate array
        self.rate = rate                            # total event rate
        
        self.parent = parent                        # parent node
        self.S_u = S_u                              # edge (S_u, S_v) into which
        self.S_v = S_v                              # the branch is embedded
        self.transferred = transferred              # 1 if HGT edge, 0 otherwise
        

class GeneTreeSimulator:
    
    def __init__(self, S):
        
        if not isinstance(S, Tree) or not S.root:
            raise TypeError("'S' must be a non-empty tree of type 'Tree'")
        
        self.S = S                                  # species tree
        self.sorted_speciations = sorted_nodes(S)   # list of speciations sorted by time stamp
        self.sorted_edges = sorted_edges(S)         # edges of the species tree
                                                    # sort (u,v) by tstamp of u
        self._analyze_secies_tree()
    
    
    def _analyze_secies_tree(self):
        
        self.S_subtree_survivors = {}
        
        for u in self.S.postorder():
            
            if not u.children:
                if u.event == 'L':
                    self.S_subtree_survivors[u] = 0
                else:
                    self.S_subtree_survivors[u] = 1
            
            else:
                self.S_subtree_survivors[u] = sum(self.S_subtree_survivors[v]
                                                  for v in u.children)
        
        if not self.S_subtree_survivors[self.S.root]:
            warnings.warn('species tree has no non-loss leaves',
                          category=warnings.UserWarning)
    
    
    def simulate(self,
                 dupl_rate=0.0,
                 loss_rate=0.0,
                 hgt_rate=0.0,
                 dupl_polytomy=0.0,
                 prohibit_extinction='per_species',
                 replace_prob=0.0,
                 **kwargs):
        """Simulate a gene tree along the specified species tree.
        
        Keyword arguments:
            dupl_rate -- duplication rate, default is 0.0
            loss_rate -- loss rate, default is 0.0
            hgt_rate -- horizontal gene transfer rate, default is 0.0
            dupl_polytomy -- allows non-binary duplication events by specifying
                the lambda parameter for a poisson distribution (copy number =
                drawn number + 2); default is 0.0
            prohibit_extinction -- avoid the extinction of all members in any
                species ('per_species'), of the complete gene family
                ('per_family'), or no constraints (False); default is 
                'per_species'.
            replace_prob -- replacing HGT events, probability by which one
                random homolog in the receiving branch of the receiving branch
                gets lost immediately after the transfer
        """
        
        self.rate_sum = dupl_rate + loss_rate + hgt_rate
        self.d = dupl_rate
        self.l = loss_rate
        self.h = hgt_rate
        
        self._prohibit_extinction = prohibit_extinction
        
        self._dupl_polytomy = dupl_polytomy
        
        self._replace_prob = replace_prob
        
        self._reset()
        
        return self._run()
    
    
    def _reset(self):
        
        # queue for speciation events
        self.spec_queue = deque(self.sorted_speciations)
        self.id_counter = 0
        
        # counter for surviving genes
        self.total_surviving = 0
        
        # total event rate (all branches)
        self.total_rate = 0
        
        # keep track of surving branches that are in species branches with
        # at least 1 surviving species leaf
        self.surv_non_loss_lineages = set()
        
        self.branches = []
        
        # maps E(S) --> existing branches
        self.ES_to_b = {e: [] for e in self.sorted_edges}

    
    def _get_tstamp(self, t):
        
        rate = self.total_rate
        if (self._prohibit_extinction == 'per_family' and
            len(self.surv_non_loss_lineages) == 1):
            rate -= self.l
        
        if rate <= 0.0:
            return -1
        else:
            return t - np.random.exponential(1/rate)
    
    
    def _get_branch_and_type(self):
        
        if (self._prohibit_extinction == 'per_family' and
            len(self.surv_non_loss_lineages) == 1):
            
            # get the single branch that is in a non-loss species branch
            special_branch = next(iter(self.surv_non_loss_lineages))
            temp_rate = special_branch.rate
            
            r = np.random.uniform(high=self.total_rate-self.l)
            current_sum = 0.0
            
            for i in range(len(self.branches)):
                if r <= current_sum + self.branches[i].rate:
                    break
                current_sum += self.branches[i].rate
            
            branch = self.branches[i]
            
            loss_factor = 0 if branch is special_branch else 1
            
            if r <= current_sum + self.d:
                event_type = 'D'
            elif r <= current_sum + self.d + loss_factor * self.l:
                event_type = 'L'
            else:
                event_type = 'H'
            
            # finally reset branch rate
            special_branch.rate = temp_rate
        
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
        
        # assert np.isclose(self.total_rate, sum([b.rate for b in self.branches])), "Sum of rates and total rate are not equal."
        
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


    def _run(self):
        
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
                    trans_edge, _ = self._hgt(event_tstamp, branch)
                    
                    if (self._replace_prob > 0.0            and
                        trans_edge                          and
                        len(self.ES_to_b[trans_edge]) > 1   and
                        np.random.random() < self._replace_prob):
                        
                        # choose among all except the last added branch
                        i = np.random.randint(0, 
                                      high=len(self.ES_to_b[trans_edge])-1)
                        self._loss(event_tstamp,
                                   self.ES_to_b[trans_edge][i])
                
                t = event_tstamp

        return T
    
    
    def _initiatialize_tree(self):
        
        if len(self.S.root.children) > 1:
            # root is a speciation event
            root = TreeNode(label=0, event='S',
                            color=self.S.root.label, 
                            dist=0.0, tstamp=self.S.root.tstamp)
        else:                    
            # planted species tree
            root = TreeNode(label=0, event=None,
                            color=self.S.root.label,
                            dist=0.0, tstamp=self.S.root.tstamp)
            
        T = Tree(root)
        self.id_counter += 1
        self.spec_queue.popleft()
        
        if self._prohibit_extinction == 'per_species':
            rate = self.d + self.h
        else:
            rate = self.rate_sum
                
        for S_v in self.S.root.children:
            
            array_id = len(self.branches)
            new_branch = _Branch(self.id_counter, array_id, rate,
                                 T.root, self.S.root, S_v, 0)
            self.ES_to_b[(self.S.root, S_v)].append(new_branch)
            self.branches.append(new_branch)
            self.id_counter += 1
            if self.S_subtree_survivors[S_v]:
                self.surv_non_loss_lineages.add(new_branch)
                
            self.total_rate += rate   
            self.total_surviving = 1
            
        return T
            
    
    def _speciation(self):
        
        # also handles loss and leaf nodes of the species tree
        
        S_v = self.spec_queue.popleft()
        S_u = S_v.parent
        
        # copy since we modify this list
        branches = self.ES_to_b[(S_u, S_v)].copy()
        
        for branch in branches:
            
            spec_node = TreeNode(label=branch.label,
                                 event='S',
                                 color=S_v.label, tstamp=S_v.tstamp,
                                 dist=abs(S_v.tstamp-branch.parent.tstamp),
                                 transferred=branch.transferred)
            branch.parent.add_child(spec_node)
            
            for S_w in S_v.children:
                
                if S_w is S_v.children[0]:
                    array_id = branch.array_id
                    new_branch = _Branch(self.id_counter, array_id, branch.rate,
                                         spec_node, S_v, S_w, 0)
                    self.branches[array_id] = new_branch
                else:
                    array_id = len(self.branches)
                    new_branch = _Branch(self.id_counter, array_id, branch.rate,
                                         spec_node, S_v, S_w, 0)
                    self.branches.append(new_branch)
                
                self.ES_to_b[(S_v, S_w)].append(new_branch)
                self.id_counter += 1
                if self.S_subtree_survivors[S_w]:
                    self.surv_non_loss_lineages.add(new_branch)
            
            # losses and (extant) leaves
            if not S_v.children:
                
                self.ES_to_b[(S_u, S_v)].remove(branch)
                self.total_rate -= branch.rate
                branch.rate = 0.0
                
                if S_v.label == 'L':
                    spec_node.event = 'L'
                    self.total_surviving -= 1
                    self.surv_non_loss_lineages.discard(branch)
            
            else:
                self.total_surviving += len(S_v.children) - 1
                self.total_rate += (len(S_v.children) - 1) * branch.rate
                self.surv_non_loss_lineages.discard(branch)
            
    
    def _duplication(self, event_tstamp, branch):
        
        S_u, S_v = branch.S_u, branch.S_v
        
        dupl_node = TreeNode(label=branch.label,
                             event='D',
                             color=(S_u.label,S_v.label),
                             tstamp=event_tstamp,
                             dist=abs(event_tstamp-branch.parent.tstamp),
                             transferred=branch.transferred)
        branch.parent.add_child(dupl_node)
        self.ES_to_b[(S_u, S_v)].remove(branch)
        self.surv_non_loss_lineages.discard(branch)
        
        copy_number = self._get_copy_number()
        
        for i in range(copy_number):
            
            if i == 0:
                array_id = branch.array_id
                new_branch = _Branch(self.id_counter, array_id, self.rate_sum,
                                     dupl_node, S_u, S_v, 0)
                self.branches[array_id] = new_branch
            else:
                array_id = len(self.branches)
                new_branch = _Branch(self.id_counter, array_id, self.rate_sum,
                                     dupl_node, S_u, S_v, 0)
                self.branches.append(new_branch)
            
            self.ES_to_b[(S_u, S_v)].append(new_branch)
            self.id_counter += 1
            if self.S_subtree_survivors[S_v]:
                self.surv_non_loss_lineages.add(new_branch)
        
        self.total_surviving += (copy_number - 1)
        self.total_rate += (copy_number - 1) * (self.rate_sum)
        
        if (self._prohibit_extinction == 'per_species' and
            len(self.ES_to_b[(S_u, S_v)]) == copy_number):
            self.total_rate += self.l
            
            
    def _loss(self, event_tstamp, branch):
        
        S_u, S_v = branch.S_u, branch.S_v
        
        loss_node = TreeNode(label=branch.label,
                             event='L',
                             color=(S_u.label,S_v.label),
                             tstamp=event_tstamp,
                             dist=abs(event_tstamp-branch.parent.tstamp),
                             transferred=branch.transferred)
        branch.parent.add_child(loss_node)
        self.ES_to_b[(S_u, S_v)].remove(branch)
        self.surv_non_loss_lineages.discard(branch)
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
            hgt_node = TreeNode(label=branch.label,
                                event='H',
                                color=(S_u.label,S_v.label),
                                tstamp=event_tstamp,
                                dist=abs(event_tstamp-branch.parent.tstamp),
                                transferred=branch.transferred)
            branch.parent.add_child(hgt_node)
            self.ES_to_b[(S_u, S_v)].remove(branch)
            self.surv_non_loss_lineages.discard(branch)
            
            # original branch
            array_id = branch.array_id
            new_branch = _Branch(self.id_counter, array_id, branch.rate,
                                 hgt_node, S_u, S_v, 0) 
            self.branches[array_id] = new_branch
            self.ES_to_b[(S_u, S_v)].append(new_branch)
            self.id_counter += 1
            if self.S_subtree_survivors[S_v]:
                self.surv_non_loss_lineages.add(new_branch)
            
            # receiving branch
            array_id = len(self.branches)
            trans_branch = _Branch(self.id_counter, array_id, self.rate_sum,
                                   hgt_node, *trans_edge, 1)
            self.branches.append(trans_branch)
            self.ES_to_b[trans_edge].append(trans_branch)
            self.id_counter += 1
            if self.S_subtree_survivors[trans_edge[1]]:
                self.surv_non_loss_lineages.add(trans_branch)
            
            self.total_surviving += 1
            self.total_rate += self.rate_sum
            
            if (self._prohibit_extinction == 'per_species' and
                len(self.ES_to_b[trans_edge]) == 2):
                self.total_rate += self.l
                self.ES_to_b[trans_edge][0].rate += self.l
            
            return trans_edge, new_branch
        
        else:
            return False, None
                
    
    def _assert_no_extinction(self, T):
        """Returns False if gene family is extinct in some species."""
        
        VS_to_VT = {l.label: [] for l in self.S.preorder() if not l.children and
                                                                  l.event != 'L'}
        
        for v in T.preorder():
            if not v.children and v.event != 'L':
                VS_to_VT[v.color].append(v.label)
                
        for leaf_list in VS_to_VT.values():
            if not leaf_list:
                return False
        
        return True
        

# --------------------------------------------------------------------------
#                 CONSTRUCTION OF THE OBSERVABLE TREE
# --------------------------------------------------------------------------

def observable_tree(tree):
    
    obs_tree = delete_losses_and_contract(tree, inplace=False)
    
    remove_planted_root(obs_tree, inplace=True)
    
    return obs_tree