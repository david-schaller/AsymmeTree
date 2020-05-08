# -*- coding: utf-8 -*-

"""
Tree Imbalancer.

Introduce evolution rate asymmetries and autocorrelation.

Methods in this module:
    - imbalance_tree
"""

import numpy as np


__author__ = "David Schaller"


# --------------------------------------------------------------------------
#                      IMBALANCING OF THE GENE TREE
#
# --------------------------------------------------------------------------


def imbalance_tree(T, S, baseline_rate=1.0,
                   autocorr_factors=None,
                   autocorr_variance=0.0,
                   gamma_param=(0.5, 1.0, 2.2),
                   CSN_weights=(1, 1, 1),
                   inplace=True):
    """Imbalances an (ultrametric) TRUE gene tree.
    
    Keyword arguments:
    baseline_rate -- mean of substitution rate for conserved genes
    autocorr_factors -- autocorrelation rate factors for the edges of S
    autocorr_variance -- autocorrelation variance factor for lognormal
        distribution, only relevant if 'autocorrelation_rates' are not supplied 
    gamma_param -- param. for gamma distribution (a, loc, scale)
    CSN_weights -- weights for choice between conservation, subfunctionalization
        and neofunctionalization
    inplace -- if False, copy the tree before imbalancing
    """
    
    if not inplace:
        T = T.copy()
     
    # factors for subfunctionalization/neofunctionalization
    CSN_weights = np.asarray(CSN_weights) / sum(CSN_weights)
    _divergent_rates(T, S, gamma_param, CSN_weights)
    
    # autocorrelation
    if autocorr_factors:
        _apply_autocorrelation(T, autocorr_factors, inplace=True)
    elif autocorr_variance > 0.0:
        _, edge_rates = autocorrelation_factors(S, autocorr_variance)
        _apply_autocorrelation(T, edge_rates, inplace=True)
    
    # finally apply baseline rate
    for v in T.preorder():
        v.dist *= baseline_rate
    
    return T


# --------------------------------------------------------------------------
#                       EVOLUTION RATE ASYMMETRY
#
# --------------------------------------------------------------------------

def _adjust_distances(T, rates):
    
    for edge, rate_list in rates.items():
        time_points = np.asarray([tstamp for tstamp, _ in rate_list] + [edge[1].tstamp])
        rate_values = np.asarray([rate for _, rate in rate_list])
        edge[1].dist = np.dot(-np.diff(time_points), rate_values)
        

def _get_rate(marked_as, gamma_param):
    
    if marked_as == "conserved":
        # --- conserved genes are set to rate 1 ---
        return 1.0
    
    elif marked_as == "divergent":
        # --- fitted gamma distribution ---
        return (np.random.gamma(gamma_param[0], scale=gamma_param[2]) +
                gamma_param[1])
        

def _duplication_type(marked_as, CSN_weights):
    
    if marked_as == "divergent":
        return "divergent", "divergent"
    else:
        r = np.random.choice(3, p=CSN_weights)
        if r == 0:                                  # conservation
            return "conserved", "conserved"
        elif r == 1:                                # subfunctionalization
            return "divergent", "divergent"
        else:                                       # neofunctionalization
            if np.random.uniform() < 0.5:
                return "divergent", "conserved"
            else:
                return "conserved", "divergent"


def _divergent_rates(T, S, gamma_param, CSN_weights):
    """
    Keyword arguments:
    gamma_param -- param. for gamma distribution (a, loc, scale)
    CSN_weights -- weights for choice between conservation,
        subfunctionalization and neofunctionalization
    """
    
    T_nodes = T.sorted_nodes()
    rates = {edge: [] for edge in T.edges()}        # edge --> list of (tstamp, rate) tuples
    
    S_parents = {v.ID: v.parent.ID for v in S.preorder() if v.parent}
    gene_counter = {(e[0].ID, e[1].ID): [] for e in S.edges()}
    marked = {v: "conserved" for v in T_nodes}      # marked as conserved or divergent
    
    for u in T_nodes:
        
        # ----------------- SPECIATION -----------------
        if u.label in ("S", ""):
            for v in u.children:
                marked[v] = marked[u]
                S_u = u.color
                S_v = v.color if not isinstance(v.color, (tuple, list)) else v.color[1]
                gene_counter[(S_u, S_v)].append(v)
                rates[(u,v)].append((u.tstamp, _get_rate(marked[v], gamma_param)))
#                rates[(u,v)].append((u.tstamp, rates[(u.parent,u)][-1][1] if u.parent else 1.0))
            
        # ---------------- DUPLICATION -----------------
        elif u.label == "D":
            marked[u.children[0]], marked[u.children[1]] = _duplication_type(marked[u],
                                                                             CSN_weights)
            gene_counter[u.color].remove(u)
            for v in u.children:
                gene_counter[u.color].append(v)
                rates[(u,v)].append((u.tstamp, _get_rate(marked[v], gamma_param)))
#                if marked[u] == "conserved":
#                    rates[(u,v)].append((u.tstamp, get_rate(marked[v]), gamma_param))
#                else:
#                    rates[(u,v)].append((u.tstamp, rates[(u.parent,u)][-1][1]))
        
        # ------------------- LOSS ---------------------
        elif u.is_loss():
            gene_counter[u.color].remove(u)
            if len(gene_counter[u.color]) == 1:
                v = gene_counter[u.color][0]
                if marked[v] == "divergent":
                    marked[v] = "conserved"
                    rates[(v.parent,v)].append((u.tstamp, _get_rate(marked[v], gamma_param)))
        
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
                rates[(u,v1)].append((u.tstamp, _get_rate(marked[v1], gamma_param)))
            
            # transferred copy
            marked[v2] = "divergent"
            if isinstance(v2.color, (tuple, list)):
                gene_counter[v2.color].append(v2)
            else:
                gene_counter[(S_parents[v2.color], v2.color)].append(v2)
            rates[(u,v2)].append((u.tstamp, _get_rate(marked[v2], gamma_param)))
            
    _adjust_distances(T, rates)
    return T

# --------------------------------------------------------------------------
#                         AUTOCORRELATION
#
# --------------------------------------------------------------------------
    
def autocorrelation_factors(tree, variance):
    """Geometric Brownian motion process to assign rate factors to species tree."""
    
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
