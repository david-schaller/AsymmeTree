# -*- coding: utf-8 -*-

"""
Tree Imbalancer.

Introduce evolution rate asymmetries and autocorrelation.

Methods in this module:
    - imbalance_tree
"""

import numpy as np

from asymmetree.tools.PhyloTree import PhyloTree


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


# --------------------------------------------------------------------------
#                      IMBALANCING OF THE GENE TREE
#
# --------------------------------------------------------------------------


def imbalance_tree(T, S, baseline_rate=1,
                   lognormal_v=0,
                   gamma_param=(0.5, 1.0, 2.2),
                   weights=(1/3, 1/3, 1/3), copy_tree=False):
    """Imbalances an (ultrametric) TRUE gene tree.
    
    Keyword arguments:
    baseline_rate   -- mean of substitution rate for conserved genes
    lognormal_v     -- variance factor for lognormal distribution
    gamma_param     -- param. for gamma distribution (a, loc, scale)
    weights         -- weights for choice between conservation,
                       subfunctionalization and neofunctionalization
    copy_tree       -- copy the tree before imbalancing
    """
    
    if copy_tree:
        T = T.copy()
    weights = np.asarray(weights) / sum(weights)
        
    _divergent_rates(T, S, gamma_param, weights)
    
    _autocorrelation(T, S, baseline_rate, lognormal_v)
    
    return T


# --------------------------------------------------------------------------
#                       EVOLUTION RATE ASYMMETRY
#
# --------------------------------------------------------------------------

def _adjust_distances(T, rates):
    for edge, rate_list in rates.items():
        time_points = np.asarray([tstamp for tstamp, _ in rate_list] + [edge[1].tstamp])
        rate_values = np.asarray([rate for _, rate in rate_list])
        edge[1].dist = np.dot(-np.diff(time_points),rate_values)


def _divergent_rates(T, S, gamma_param, weights):
    """
    Keyword arguments:
    gamma_param     -- param. for gamma distribution (a, loc, scale)
    weights         -- weights for choice between conservation,
                       subfunctionalization and neofunctionalization
    """
    
    def get_rate(marked_as):
        if marked_as == "conserved":
            # --- conserved genes are set to rate 1 ---
            return 1.0
        elif marked_as == "divergent":
            # --- fitted gamma distribution ---
            return (np.random.gamma(gamma_param[0], scale=gamma_param[2]) +
                    gamma_param[1])
    
    def duplication_type(marked_as):
        if marked_as == "divergent":
            return "divergent", "divergent"
        else:
            r = np.random.choice(3, p=weights)
            if r == 0:                                  # conservation
                return "conserved", "conserved"
            elif r == 1:                                # subfunctionalization
                return "divergent", "divergent"
            else:                                       # neofunctionalization
                if np.random.uniform() < 0.5:
                    return "divergent", "conserved"
                else:
                    return "conserved", "divergent"
    
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
                S_v = v.color if not isinstance(v.color, tuple) else v.color[1]
                gene_counter[(S_u, S_v)].append(v)
                rates[(u,v)].append((u.tstamp, get_rate(marked[v])))
#                rates[(u,v)].append((u.tstamp, rates[(u.parent,u)][-1][1] if u.parent else 1.0))
            
        # ---------------- DUPLICATION -----------------
        elif u.label == "D":
            marked[u.children[0]], marked[u.children[1]] = duplication_type(marked[u])
            gene_counter[u.color].remove(u)
            for v in u.children:
                gene_counter[u.color].append(v)
                rates[(u,v)].append((u.tstamp, get_rate(marked[v])))
#                if marked[u] == "conserved":
#                    rates[(u,v)].append((u.tstamp, get_rate(marked[v])))
#                else:
#                    rates[(u,v)].append((u.tstamp, rates[(u.parent,u)][-1][1]))
        
        # ------------------- LOSS ---------------------
        elif u.label == "*":
            gene_counter[u.color].remove(u)
            if len(gene_counter[u.color]) == 1:
                v = gene_counter[u.color][0]
                if marked[v] == "divergent":
                    marked[v] = "conserved"
                    rates[(v.parent,v)].append((u.tstamp, get_rate(marked[v])))
        
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
                rates[(u,v1)].append((u.tstamp, get_rate(marked[v1])))
            
            # transferred copy
            marked[v2] = "divergent"
            if isinstance(v2.color, tuple):
                gene_counter[v2.color].append(v2)
            else:
                gene_counter[(S_parents[v2.color], v2.color)].append(v2)
            rates[(u,v2)].append((u.tstamp, get_rate(marked[v2])))
            
    _adjust_distances(T, rates)
    return T

# --------------------------------------------------------------------------
#                         AUTOCORRELATION
#
# --------------------------------------------------------------------------

def _autocorrelation(T, S, baseline_rate, lognormal_v):
    """
    Keyword arguments:
    baseline_rate   -- mean of substitution rate for conserved genes
    lognormal_v     -- variance factor for lognormal distribution
    """
    rates = {}                                                          # maps node v --> rate
    S_parent = {}
    for v in S.preorder():
        if not v.parent:
            rates[v.ID] = baseline_rate
        else:
            var = lognormal_v * v.dist
            mu = np.log(rates[v.parent.ID]) - var/2                     # exp. value equal to parent's rate
            rates[v.ID] = np.exp(np.random.normal(mu, np.sqrt(var)))
            S_parent[v.ID] = v.parent.ID                                # save parent's ID for below
    
    for v in T.preorder():
        if v.parent:
            S_v_ID = v.color[1] if isinstance(v.color, tuple) else v.color
            effective_rate = (rates[S_parent[S_v_ID]] + rates[S_v_ID]) / 2
            v.dist = v.dist * effective_rate
