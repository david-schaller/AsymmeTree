# -*- coding: utf-8 -*-

"""
Assign asymmetric evolution rates to the branches of a gene tree.

Introduce evolution rate asymmetries and autocorrelation.
"""

import numpy as np

from asymmetree.treeevolve.GeneTree import GeneTreeSimulator
from asymmetree.tools.Sampling import Sampler


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                         USER INTERFACE FUNCTION
# --------------------------------------------------------------------------

def simulate_gene_trees(S, N=1,
                        dupl_rate=0.0,
                        loss_rate=0.0,
                        hgt_rate=0.0,
                        base_rate=1.0,
                        **kwargs):
    """Simulates dated gene trees with non-ultrametric edge lengths along a
    species tree.
    
    Keyword arguments:
        N -- number of gene trees to be simulated, default is 1, in which case
            a tree is returned, otherwise a list is returned
        dupl_rate -- (distribution for the) duplication rate,
            default is constant 0.0
        loss_rate -- (distribution for the) loss rate,
            default is constant 0.0
        hgt_rate -- (distribution for the) HGT rate,
            default is constant 0.0
        base_rate -- (distribution for the) evolution rate at the roots of
            the gene trees, default is constant 1.0
        kwargs -- see arguments of GeneTreeSimulator.simulate and assign_rates
        """
    
    gene_trees = []
    simulator = GeneTreeSimulator(S)
    
    dupl_rate_sampler = Sampler(dupl_rate)
    loss_rate_sampler = Sampler(loss_rate)
    hgt_rate_sampler = Sampler(hgt_rate)
    base_rate_sampler = Sampler(dupl_rate)
    
    # autocorrelation between genes of the same or related species
    autocorr_variance = kwargs.pop('autocorr_variance', 0.0)
    _, autocorr_factors = autocorrelation_factors(S, autocorr_variance)
    
    # main simulation and imbalancing
    for i in range(N):
        
        TGT = simulator.simulate(dupl_rate=dupl_rate_sampler.draw(),
                                 loss_rate=loss_rate_sampler.draw(),
                                 hgt_rate=hgt_rate_sampler.draw(),
                                 **kwargs)
        assign_rates(TGT, S,
                     base_rate=base_rate_sampler.draw(),
                     autocorr_factors=autocorr_factors,
                     **kwargs)
        gene_trees.append(TGT)
    
    if N == 1:
        return gene_trees[0]
    else:
        return gene_trees
    

    

# --------------------------------------------------------------------------
#                      IMBALANCING OF THE GENE TREE
# --------------------------------------------------------------------------


def assign_rates(T, S, base_rate=1.0,
                 autocorr_factors=None,
                 autocorr_variance=0.0,
                 rate_increase=('gamma', 0.5, 2.2),
                 CSN_weights=(1, 1, 1),
                 inplace=True,
                 **kwargs):
    """Assigns realistic evolution rates to a TRUE gene tree.
    
    The assigned rates are used to modify the length ('dist') of the edges of
    the (originally ultrametric) dated gene tree.
    
    Keyword arguments:
    base_rate -- mean of substitution rate for conserved genes
    autocorr_factors -- autocorrelation rate factors for the edges of S
    autocorr_variance -- autocorrelation variance factor for lognormal
        distribution, only relevant if 'autocorrelation_rates' are not supplied 
    rate_increase -- distribution of the (relative) rate increase (w.r.t. the
        base rate) for divergent genes, i.e. to a factor 1 + x, default is gamma
        distribution with shape 0.5 and scale 2.2
    CSN_weights -- weights for choice between conservation, subfunctionalization
        and neofunctionalization
    inplace -- if False, copy the tree before imbalancing
    """
    
    if not inplace:
        T = T.copy()
     
    # factors for subfunctionalization/neofunctionalization
    CSN_weights = np.asarray(CSN_weights) / sum(CSN_weights)
    sampler = Sampler(rate_increase, shift=1.0)
    _divergent_rates(T, S, sampler, CSN_weights)
    
    # autocorrelation
    if autocorr_factors:
        _apply_autocorrelation(T, autocorr_factors, inplace=True)
    elif autocorr_variance > 0.0:
        _, edge_rates = autocorrelation_factors(S, autocorr_variance)
        _apply_autocorrelation(T, edge_rates, inplace=True)
    
    # finally apply base rate
    for v in T.preorder():
        v.dist *= base_rate
    
    return T


# --------------------------------------------------------------------------
#                       EVOLUTION RATE ASYMMETRY
# --------------------------------------------------------------------------

def _adjust_distances(T, rates):
    
    for edge, rate_list in rates.items():
        time_points = np.asarray([tstamp for tstamp, _ in rate_list] + [edge[1].tstamp])
        rate_values = np.asarray([rate for _, rate in rate_list])
        edge[1].dist = np.dot(-np.diff(time_points), rate_values)
        

def _duplication_type(marked_as, CSN_weights):
    
    if marked_as == 'divergent':
        return 'divergent', 'divergent'
    else:
        r = np.random.choice(3, p=CSN_weights)
        if r == 0:                                  # conservation
            return 'conserved', 'conserved'
        elif r == 1:                                # subfunctionalization
            return 'divergent', 'divergent'
        else:                                       # neofunctionalization
            if np.random.uniform() < 0.5:
                return 'divergent', 'conserved'
            else:
                return 'conserved', 'divergent'


def _divergent_rates(T, S, sampler, CSN_weights):
    """
    Parameters:
        sampler -- sampler for rate increase for divergent genes
        CSN_weights -- weights for choice between conservation,
            subfunctionalization and neofunctionalization
    """
    
    T_nodes = T.sorted_nodes()
    rates = {edge: [] for edge in T.edges()}        # edge --> list of (tstamp, rate) tuples
    
    S_parents = {v.ID: v.parent.ID for v in S.preorder() if v.parent}
    gene_counter = {(e[0].ID, e[1].ID): [] for e in S.edges()}
    marked = {v: 'conserved' for v in T_nodes}      # marked as conserved or divergent
    
    for u in T_nodes:
        
        # ----------------- SPECIATION -----------------
        if u.label in ('S', ''):
            for v in u.children:
                marked[v] = marked[u]
                S_u = u.color
                S_v = v.color if not isinstance(v.color, (tuple, list)) else v.color[1]
                gene_counter[(S_u, S_v)].append(v)
                new_rate = sampler.draw() if marked[v] == 'divergent' else 1.0
                rates[(u,v)].append((u.tstamp, new_rate))
            
        # ---------------- DUPLICATION -----------------
        elif u.label == "D":
            marked[u.children[0]], marked[u.children[1]] = _duplication_type(marked[u],
                                                                             CSN_weights)
            gene_counter[u.color].remove(u)
            for v in u.children:
                gene_counter[u.color].append(v)
                new_rate = sampler.draw() if marked[v] == 'divergent' else 1.0
                rates[(u,v)].append((u.tstamp, new_rate))
        
        # ------------------- LOSS ---------------------
        elif u.is_loss():
            gene_counter[u.color].remove(u)
            if len(gene_counter[u.color]) == 1:
                v = gene_counter[u.color][0]
                if marked[v] == 'divergent':
                    marked[v] = 'conserved'
                    rates[(v.parent,v)].append((u.tstamp, 1.0))
        
        # ---------- HORIZONTAL GENE TRANSFER ----------
        elif u.label == "H":
            v1, v2 = u.children
            if v1.transferred:
                v1, v2 = v2, v1         # now v2 is the transferred copy
                
            # untransferred copy
            marked[v1] = marked[u]
            gene_counter[u.color].remove(u)
            gene_counter[u.color].append(v1)
            if u.parent:
                rates[(u,v1)].append((u.tstamp, rates[(u.parent,u)][-1][1]))
            else:
                new_rate = sampler.draw() if marked[v1] == 'divergent' else 1.0
                rates[(u,v1)].append((u.tstamp, new_rate))
            
            # transferred copy
            marked[v2] = 'divergent'
            if isinstance(v2.color, (tuple, list)):
                gene_counter[v2.color].append(v2)
            else:
                gene_counter[(S_parents[v2.color], v2.color)].append(v2)
            new_rate = sampler.draw() if marked[v2] == 'divergent' else 1.0
            rates[(u,v2)].append((u.tstamp, new_rate))
            
    _adjust_distances(T, rates)
    return T

# --------------------------------------------------------------------------
#                         AUTOCORRELATION
# --------------------------------------------------------------------------
    
def autocorrelation_factors(tree, variance):
    """Geometric Brownian motion process to assign rate factors to species tree.
    
    The parameter 'variance' is a hyperparameter for a log-normal distribution
    from which offspring rates are drawn. The overall variance of this
    distribution is 'variance' * divergence time.
    """
    
    node_rates = {}                 # maps node v --> rate of v
    edge_rates = {}                 # maps v of edge (u,v) --> rate of (u,v)
    
    for v in tree.preorder():
        if not v.parent:
            # assign factor 1.0 to root (= expected value for all other nodes
            # and edges)
            node_rates[v.ID] = 1.0
            edge_rates[v.ID] = 1.0
        else:
            var = variance * v.dist
            # ensure that exp. value is equal to parent's rate
            mu = np.log(node_rates[v.parent.ID]) - var/2
            
            node_rates[v.ID] = np.exp(np.random.normal(mu, np.sqrt(var)))
            
            # edge rate as arithmetic mean of u and v
            edge_rates[v.ID] = (node_rates[v.parent.ID] + node_rates[v.ID]) / 2
            
    return node_rates, edge_rates


def _apply_autocorrelation(T, edge_rates, inplace=True):
    
    if not inplace:
        T = T.copy()
    
    for v in T.preorder():
        if v.parent:
            edge_ID = v.color[1] if isinstance(v.color, (tuple, list)) else v.color
            v.dist *= edge_rates[edge_ID]
    
    return T
