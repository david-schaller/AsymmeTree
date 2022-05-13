# -*- coding: utf-8 -*-

"""
Simulation of dated gene trees.
"""

import random, warnings
from collections import deque
from dataclasses import dataclass
from warnings import warn

import numpy as np

from tralda.datastructures.Tree import Tree, TreeNode, LCA

from asymmetree.tools.PhyloTreeTools import (sorted_nodes,
                                             sorted_edges,
                                             delete_losses_and_contract,
                                             remove_planted_root,
                                             distance_from_timing)


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                    CONSTRUCTION OF THE GENE TREE
# --------------------------------------------------------------------------
    

def dated_gene_tree(S, **kwargs):
    """Simulate a gene tree along the specified species tree.
    
    Parameters
    ----------
    dupl_rate : float, optional
        Duplication rate, default is 0.0.
    loss_rate : float, optional
        Loss rate, default is 0.0.
    hgt_rate : float, optional
        Horizontal gene transfer rate, default is 0.0.
    gc_rate : float, optional
        Gene conversion rate, default is 0.0.
    dupl_polytomy : float, optional
        Allows non-binary duplication events by specifying the lambda
        parameter for a poisson distribution (copy number =
        drawn number + 2); default is 0.0.
    prohibit_extinction : str or bool, optional 
        Avoid the extinction of all members in any species ('per_species'),
        of the complete gene family ('per_family'), or no constraints
        (False); default is 'per_species'.
    replace_prob : float, optional
        Enables replacing HGT events, probability by which one random
        homolog in the receiving branch of the receiving branch gets lost
        immediately after the transfer.
    additive_transfer_distance_bias : str or bool, optional
        Specifies whether closer related species have a higher probability
        to be the recipient species in an additive HGT event. The default
        is False, in which case the recipient species is chosen at random
        among the co-existing species. The options 'inverse' and
        'exponential' mean that a species branch is sampled weighted by
        1/(a * t) or e^(-(a * t)), resp., where t is the elapsed time between
        the last common ancestor of the two species branches and the time of
        the event, see [1], and a is a user-defined factor.
    replacing_transfer_distance_bias : str or bool, optional
        Specifies whether closer related gene branches have a higher
        probability to be replaced in a replacing HGT event. The default
        is False, in which case the replaced gene is chosen at random
        among the co-existing gene branches. The options 'inverse' and
        'exponential' mean that a gene branch is sampled weighted by
        1/(a * t) or e^(-(a * t)), resp., where t is the elapsed time between
        the last common ancestor of the two gene branches and the time of the
        event, see [1], and a is a user-defined factor.
    transfer_distance_bias : str or bool, optional
        Set a common bias mode for additive and replacing HGT, see
        description of parameters 'additive_transfer_distance_bias' and
        'replacing_transfer_distance_bias'. If the latter are no set to
        the default (False), then these optioned are prioritized.
    transfer_distance_bias_strength : float, optional
        Intensity of the transfer distance bias (factor a) for additive and
        replacing HGT. The default is 1.0.
    gc_distance_bias : str or bool, optional
        Specifies whether closer related gene branches have a higher
        probability to be replaced in a gene conversion event. The default
        is False, in which case the replaced gene is chosen at random
        among the paralogs in the respective species lineage. The options
        'inverse' and 'exponential' mean that a paralog is sampled weighted
        by 1/(a * t) or e^(-(a * t)), resp., where t is the elapsed time
        between the last common ancestor of the two gene branches and the
        time of the event, see [1], and a is a user-defined factor.
    gc_distance_bias_strength : float, optional
        Intensity of the distance bias (factor a) for gene conversion.
        The default is 1.0.
    
    Returns
    -------
    Tree
        The simulated gene tree.
    
    References
    ----------
    .. [1] S. Kundu, M. S. Bansal.
       SaGePhy: an improved phylogenetic simulation framework for gene and 
       subgene evolution.
       In: Bioinformatics, 35(18), 2019, 3496–3498.
       doi:10.1093/bioinformatics/btz081.
    """
    
    gene_tree_simulator = GeneTreeSimulator(S)
    return gene_tree_simulator.simulate(**kwargs)


@dataclass
class _Branch:
    
    label:          int         # unique branch label
    list_id:        int         # index in list of branches
    parent:         TreeNode    # parent node
    S_edge:         TreeNode    # v of species tree edge (u, v) into which
                                # the branch is embedded
    transferred:    TreeNode    # 1 if HGT edge, 0 otherwise
    
    def __hash__(self):
        
        return hash(self.label)
        

class GeneTreeSimulator:
    
    def __init__(self, S):
        
        if not isinstance(S, Tree) or not S.root:
            raise TypeError("'S' must be a non-empty tree of type 'Tree'")
        
        self.S = S                                  # species tree
        self.lca_S = LCA(S)                         # lca data struct. for the species tree
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
                 gc_rate=0.0,
                 dupl_polytomy=0.0,
                 prohibit_extinction='per_species',
                 replace_prob=0.0,
                 additive_transfer_distance_bias=False,
                 replacing_transfer_distance_bias=False,
                 transfer_distance_bias=False,
                 transfer_distance_bias_strength=1.0,
                 gc_distance_bias=False,
                 gc_distance_bias_strength=1.0,
                 **kwargs):
        """Simulate a gene tree along the specified species tree.
        
        Parameters
        ----------
        dupl_rate : float, optional
            Duplication rate, default is 0.0.
        loss_rate : float, optional
            Loss rate, default is 0.0.
        hgt_rate : float, optional
            Horizontal gene transfer rate, default is 0.0.
        gc_rate : float, optional
            Gene conversion rate, default is 0.0.
        dupl_polytomy : float, optional
            Allows non-binary duplication events by specifying the lambda
            parameter for a poisson distribution (copy number =
            drawn number + 2); default is 0.0.
        prohibit_extinction : str or bool, optional 
            Avoid the extinction of all members in any species ('per_species'),
            of the complete gene family ('per_family'), or no constraints
            (False); default is 'per_species'.
        replace_prob : float, optional
            Enables replacing HGT events, probability by which one random
            homolog in the receiving branch of the receiving branch gets lost
            immediately after the transfer.
        additive_transfer_distance_bias : str or bool, optional
            Specifies whether closer related species have a higher probability
            to be the recipient species in an additive HGT event. The default
            is False, in which case the recipient species is chosen at random
            among the co-existing species. The options 'inverse' and
            'exponential' mean that a species branch is sampled weighted by
            1/(a * t) or e^(-(a * t)), resp., where t is the elapsed time
            between the last common ancestor of the two species branches and
            the time of the event, see [1], and a is a user-defined factor.
        replacing_transfer_distance_bias : str or bool, optional
            Specifies whether closer related gene branches have a higher
            probability to be replaced in a replacing HGT event. The default
            is False, in which case the replaced gene is chosen at random
            among the co-existing gene branches. The options 'inverse' and
            'exponential' mean that a gene branch is sampled weighted by
            1/(a * t) or e^(-(a * t)), resp., where t is the elapsed time
            between the last common ancestor of the two gene branches and the
            time of the event, see [1], and a is a user-defined factor.
        transfer_distance_bias : str or bool, optional
            Set a common bias mode for additive and replacing HGT, see
            description of parameters 'additive_transfer_distance_bias' and
            'replacing_transfer_distance_bias'. If the latter are no set to
            the default (False), then these optioned are prioritized.
        transfer_distance_bias_strength : float, optional
            Intensity of the transfer distance bias (factor a) for additive and
            replacing HGT. The default is 1.0.
        gc_distance_bias : str or bool, optional
            Specifies whether closer related gene branches have a higher
            probability to be replaced in a gene conversion event. The default
            is False, in which case the replaced gene is chosen at random
            among the paralogs in the respective species lineage. The options
            'inverse' and 'exponential' mean that a paralog is sampled weighted
            by 1/(a * t) or e^(-(a * t)), resp., where t is the elapsed time
            between the last common ancestor of the two gene branches and the
            time of the event, see [1], and a is a user-defined factor.
        gc_distance_bias_strength : float, optional
            Intensity of the distance bias (factor a) for gene conversion.
            The default is 1.0.
        
        Returns
        -------
        Tree
            The simulated gene tree.
        
        References
        ----------
        .. [1] S. Kundu, M. S. Bansal.
           SaGePhy: an improved phylogenetic simulation framework for gene and 
           subgene evolution.
           In: Bioinformatics, 35(18), 2019, 3496–3498.
           doi:10.1093/bioinformatics/btz081.
        """
        
        self.rate_sum = 0.0
        for rate in (dupl_rate, loss_rate, hgt_rate, gc_rate):
            self.rate_sum += rate
            if rate < 0.0:
                raise ValueError(f'negative rate {rate}')
        
        self.d = dupl_rate
        self.l = loss_rate
        self.h = hgt_rate
        self.gc = gc_rate
        
        if prohibit_extinction not in (False, 'per_family', 'per_species'):
            raise ValueError(f"unknown mode prohibit_extinction attribute: "\
                             f"{prohibit_extinction}")
        self._prohibit_extinction = prohibit_extinction
        
        self._dupl_polytomy = dupl_polytomy
        
        self._replace_prob = replace_prob
        
        for m in (additive_transfer_distance_bias,
                  replacing_transfer_distance_bias,
                  transfer_distance_bias,
                  gc_distance_bias):
            if m not in (False, 'inverse', 'exponential'):
                raise ValueError(f"unknown mode for transfer distance bias: "\
                                 f"{m}")
        
        # mode for transfer distance bias for additive and replacing HGT
        self._additive_transfer_distance_bias = additive_transfer_distance_bias
        self._replacing_transfer_distance_bias = replacing_transfer_distance_bias
        if transfer_distance_bias:
            if not additive_transfer_distance_bias:
                self._additive_transfer_distance_bias = transfer_distance_bias
            if not replacing_transfer_distance_bias:
                self._replacing_transfer_distance_bias = transfer_distance_bias
        
        # intensity factor for transfer distance bias
        if not isinstance(transfer_distance_bias_strength, (int, float)) or \
            transfer_distance_bias_strength <= 0.0:
            raise ValueError('factor for transfer distance bias must be > 0')
        self._transfer_distance_bias_strength = transfer_distance_bias_strength
        
        self._gc_distance_bias = gc_distance_bias
        if not isinstance(gc_distance_bias_strength, (int, float)) or \
            gc_distance_bias_strength <= 0.0:
            raise ValueError('factor for g.c. distance bias must be > 0')
        self._gc_distance_bias_strength = gc_distance_bias_strength
        
        return self._run()
    
    
    def _new_branch(self, parent, S_edge, transferred):
        
        new_branch = _Branch(self.id_counter, len(self.branches),
                             parent, S_edge, transferred)
        self.branches.append(new_branch)
        
        self.ES_to_b[S_edge].append(new_branch)
        self.id_counter += 1
        
        if self.S_subtree_survivors[S_edge]:
            self.surv_non_loss_lineages.add(new_branch)
        
        return new_branch
    
    
    def _remove_branch(self, branch):
        
        # efficient removal from list
        if branch.list_id != len(self.branches) - 1:
            branch2 = self.branches[-1]
            self.branches[branch.list_id] = branch2
            branch2.list_id = branch.list_id
        self.branches.pop()
        
        self.ES_to_b[branch.S_edge].remove(branch)
        self.surv_non_loss_lineages.discard(branch)
    
    
    def _initiatialize_tree(self):
        
        root = TreeNode(label=0, 
                        event='S' if len(self.S.root.children) > 1 else None,
                        reconc=self.S.root.label, 
                        tstamp=self.S.root.tstamp)
            
        T = Tree(root)
        self.id_counter += 1
        self.spec_queue.popleft()
                
        for S_edge in self.S.root.children:
            self._new_branch(T.root, S_edge, 0)
            
        return T
    
    
    def _get_branch_and_type(self):
        
        total_rate = len(self.branches) * self.rate_sum
        r = np.random.uniform(high=total_rate)
        i = int(len(self.branches) * r / total_rate)
        r2 = r - i * self.rate_sum
        
        if r2 <= self.d:
            event_type = 'D'
        elif r2 <= self.d + self.l:
            event_type = 'L'
        elif r2 <= self.d + self.l + self.h:
            event_type = 'H'
        else:
            event_type = 'GC'
        
        return self.branches[i], event_type


    def _run(self):
        
        # queue for speciation events
        self.spec_queue = deque(self.sorted_speciations)
        self.id_counter = 0
        
        # keep track of surving branches that are in species branches with
        # at least 1 surviving species leaf
        self.surv_non_loss_lineages = set()
        
        self.branches = []
        
        # maps species tree branches to existing gene branches
        self.ES_to_b = {S_edge: [] for _, S_edge in self.sorted_edges}
        
        self.T = self._initiatialize_tree()
        t = self.T.root.tstamp
        
        while self.spec_queue:
            
            total_rate = len(self.branches) * self.rate_sum
            event_tstamp = t - np.random.exponential(1/total_rate) \
                           if total_rate > 0.0 else -1
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
                
                # gene conversion
                elif event_type == 'GC':
                    self._gene_conversion(event_tstamp, branch)
                
                t = event_tstamp
        
        # add the 'dist' attribute to the nodes
        distance_from_timing(self.T)

        return self.T
            
    
    def _speciation(self):
        
        # also handles loss and leaf nodes of the species tree
        
        S_edge = self.spec_queue.popleft()
        
        # copy since we modify this list
        for branch in self.ES_to_b[S_edge].copy():
            
            spec_node = TreeNode(label=branch.label,
                                 event='S',
                                 reconc=S_edge.label, tstamp=S_edge.tstamp,
                                 transferred=branch.transferred)
            branch.parent.add_child(spec_node)
            self._remove_branch(branch)
            
            for S_w in S_edge.children:
                self._new_branch(spec_node, S_w, 0)
            
            # loss leaves if it was a species extinction event
            if (not S_edge.children) and S_edge.label == 'L':
                spec_node.event = 'L'
            
    
    def _duplication(self, event_tstamp, branch):
        
        S_edge = branch.S_edge
        
        dupl_node = TreeNode(label=branch.label,
                             event='D',
                             reconc=(S_edge.parent.label, S_edge.label),
                             tstamp=event_tstamp,
                             transferred=branch.transferred)
        branch.parent.add_child(dupl_node)
        self._remove_branch(branch)
        
        copy_number = 2 if self._dupl_polytomy <= 0.0 else \
                      2 + np.random.poisson(lam=self._dupl_polytomy)
        
        for i in range(copy_number):
            self._new_branch(dupl_node, S_edge, 0)
            
            
    def _loss(self, event_tstamp, branch):
        
        # not executing the loss event if extinction is prohibited is
        # equivalent to temporarily setting the loss rate to zero in the
        # respective branches
        if (self._prohibit_extinction == 'per_family' and
            len(self.surv_non_loss_lineages) == 1 and
            next(iter(self.surv_non_loss_lineages)) is branch):
            return
        
        if (self._prohibit_extinction == 'per_species' and
            len(self.ES_to_b[branch.S_edge]) <= 1):
            return
        
        S_edge = branch.S_edge
        
        loss_node = TreeNode(label=branch.label,
                             event='L',
                             reconc=(S_edge.parent.label, S_edge.label),
                             tstamp=event_tstamp,
                             transferred=branch.transferred)
        branch.parent.add_child(loss_node)
        self._remove_branch(branch)
            
    
    def _coexisting_species_edges(self, tstamp, exclude_edge=None):
        """Return list of edges for the given timestamp."""
        
        valid_species = []
        
        for S_u, S_v in self.sorted_edges:
            
            if S_u.tstamp <= tstamp:
                break
            elif S_v.tstamp < tstamp and S_v != exclude_edge:
                valid_species.append(S_v)
                
        return valid_species
    
    
    def _sample_recipient(self, event_tstamp, branch):
        
        S_edge = branch.S_edge
        trans_edge, replaced_gene_branch = None, None
        
        valid_species = self._coexisting_species_edges(event_tstamp,
                                                       exclude_edge=S_edge)
        
        if not valid_species:
            return None, None
        
        # ---- ADDITIVE HGT EVENT ----
        if np.random.random() >= self._replace_prob:
            
            # ---- no transfer distance bias ---
            if not self._additive_transfer_distance_bias:
                trans_edge = random.choice(valid_species)
            # ---- transfer distance bias ---
            else:
                distances = [(self.lca_S(S_edge, e).tstamp-event_tstamp)
                             for e in valid_species]
                a = self._transfer_distance_bias_strength
                if self._additive_transfer_distance_bias == 'inverse':
                    weights = 1 / (a * np.asarray(distances))
                elif self._additive_transfer_distance_bias == 'exponential':
                    weights = np.exp(-a * np.asarray(distances))
                
                trans_edge = random.choices(valid_species, weights=weights)[0]
            
        # ---- REPLACING HGT EVENT ----
        else:
            valid_genes = [b for e in valid_species for b in self.ES_to_b[e]]
            
            if not valid_genes:
                return None, None
            
            # ---- no transfer distance bias ---
            if not self._replacing_transfer_distance_bias:
                replaced_gene_branch = random.choice(valid_genes)
            # ---- transfer distance bias ---
            else:
                lca_T = LCA(self.T)
                distances = [(lca_T(branch.parent, b.parent).tstamp-event_tstamp)
                             for b in valid_genes]
                a = self._transfer_distance_bias_strength
                if self._replacing_transfer_distance_bias == 'inverse':
                    weights = 1 / (a * np.asarray(distances))
                elif self._replacing_transfer_distance_bias == 'exponential':
                    weights = np.exp(-a * np.asarray(distances))
                
                replaced_gene_branch = random.choices(valid_genes,
                                                      weights=weights)[0]
            
            trans_edge = replaced_gene_branch.S_edge
        
        return trans_edge, replaced_gene_branch
    
    
    def _hgt(self, event_tstamp, branch):
        
        S_edge = branch.S_edge
        
        trans_edge, replaced_gene_branch = self._sample_recipient(event_tstamp,
                                                                  branch)
        
        if trans_edge:
            hgt_node = TreeNode(label=branch.label,
                                event='H',
                                reconc=(S_edge.parent.label, S_edge.label),
                                tstamp=event_tstamp,
                                transferred=branch.transferred)
            branch.parent.add_child(hgt_node)
            self._remove_branch(branch)
            
            # original branch
            self._new_branch(hgt_node, S_edge, 0)
            
            # receiving branch
            self._new_branch(hgt_node, trans_edge, 1)
        
        # replacing HGT leads to loss in the recipient species
        if replaced_gene_branch:
            self._loss(event_tstamp, replaced_gene_branch)
            # save replaced gene information in HGT node
            hgt_node.replaced_gene = replaced_gene_branch.label
    
    
    def _sample_replaced_gene_gc(self, event_tstamp, branch):
        
        candidates = [b for b in self.ES_to_b[branch.S_edge] if b is not branch]
        
        # no gene conversion possible if there is currently no other branch in
        # the species lineage
        if not candidates:
            return
        
        if not self._gc_distance_bias:
            return random.choice(candidates)
        else:
            lca_T = LCA(self.T)
            distances = [(lca_T(branch.parent, b.parent).tstamp-event_tstamp)
                         for b in candidates]
            a = self._gc_distance_bias_strength
            if self._gc_distance_bias == 'inverse':
                weights = 1 / (a * np.asarray(distances))
            elif self._gc_distance_bias == 'exponential':
                weights = np.exp(-a * np.asarray(distances))
            
            return random.choices(candidates, weights=weights)[0]
    
    
    def _gene_conversion(self, event_tstamp, branch):
        
        S_edge = branch.S_edge
        
        replaced_gene_branch = self._sample_replaced_gene_gc(event_tstamp,
                                                             branch)
        
        if replaced_gene_branch:
            gc_node = TreeNode(label=branch.label,
                               event='GC',
                               reconc=(S_edge.parent.label, S_edge.label),
                               tstamp=event_tstamp,
                               transferred=branch.transferred)
            branch.parent.add_child(gc_node)
            self._remove_branch(branch)
            
            self._new_branch(gc_node, S_edge, 0)
            self._new_branch(gc_node, S_edge, 0)
        
            # gene conversion leads to loss of a homologous gene
            self._loss(event_tstamp, replaced_gene_branch)
            # save replaced gene information in GC node
            gc_node.replaced_gene = replaced_gene_branch.label
        

# --------------------------------------------------------------------------
#                           PRUNE LOSS BRANCHES
# --------------------------------------------------------------------------

def prune_losses(tree):
    """Prune all loss branches.
    
    Returns a copy of the tree with all branches that lead to losses only
    removed and superfluous vertices suppressed. Additionally, if the root of
    the tree has only a single child, then this 'planted edge' is also removed.    
    
    Parameters
    ----------
    tree : Tree
        The tree to be pruned.
    
    Returns
    -------
    Tree
        A pruned version of the tree.
    """
    
    pruned_tree = delete_losses_and_contract(tree, inplace=False)
    
    remove_planted_root(pruned_tree, inplace=True)
    
    return pruned_tree


def observable_tree(tree):
    
    warn('This method is deprecated. Use prune_losses() instead.',
         DeprecationWarning, stacklevel=2)
    return prune_losses(tree)
