"""Simulation of species and gene trees.

The subpackage asymmetree.treeevolve contains modules for the simulation and
manipulation of species trees and gene trees. In terms of divergence time,
these trees define an ultrametric on the set of their (extant) leaves. Gene
trees, furthermore, can be manipulated with a realistic rate heterogeneity
among their branches resulting in general additive distances (but no longer
ultrametric).

A typical simulation consists of the following steps:

(1) dated species tree (models e.g. 'Yule', and '(episodic) birth-death
    process')

(2) dated gene tree(s) (birth-death process with speciations as additional
    branching events)

(3) assignment of asymmetric evolution rates to paralogous genes

(4) pruned gene tree(s) (removal of all branches that lead to losses only)
"""

from asymmetree.treeevolve.SpeciesTree import species_tree_n as species_tree_n
from asymmetree.treeevolve.SpeciesTree import species_tree_age as species_tree_age
from asymmetree.treeevolve.SpeciesTree import species_tree_n_age as species_tree_n_age
from asymmetree.treeevolve.SpeciesTree import nonbinary as nonbinary
from asymmetree.treeevolve.GeneTree import dated_gene_tree as dated_gene_tree
from asymmetree.treeevolve.GeneTree import GeneTreeSimulator as GeneTreeSimulator
from asymmetree.treeevolve.GeneTree import prune_losses as prune_losses
from asymmetree.treeevolve.RateHeterogeneity import rate_heterogeneity as rate_heterogeneity
from asymmetree.treeevolve.RateHeterogeneity import (
    autocorrelation_factors as autocorrelation_factors,
)
from asymmetree.treeevolve.RateHeterogeneity import gene_trees as gene_trees
from asymmetree.treeevolve.DistanceNoise import noisy_matrix as noisy_matrix
from asymmetree.treeevolve.DistanceNoise import convex_linear_comb as convex_linear_comb
from asymmetree.treeevolve.DistanceNoise import wrong_topology_matrix as wrong_topology_matrix
