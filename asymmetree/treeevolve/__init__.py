"""Simulation of species and gene trees.

The subpackage asymmetree.treeevolve contains modules for the simulation and manipulation of species
trees and gene trees. A typical simulation consists of the following steps:

(1) dated species tree (models e.g. 'Yule', and '(episodic) birth-death process')

(2) dated gene tree(s) (birth-death process with speciations as additional branching events)

(3) assignment of asymmetric evolution rates to paralogous genes

(4) pruned gene tree(s) (removal of all branches that lead to losses only)

In terms of divergence time, these trees define an ultrametric on the set of their (extant) leaves.
Gene trees, furthermore, can be manipulated with a realistic rate heterogeneity among their branches
(step 3) resulting in general additive distances (but no longer ultrametric).
"""

from __future__ import annotations

from asymmetree.treeevolve.species import species_tree_n as species_tree_n
from asymmetree.treeevolve.species import species_tree_age as species_tree_age
from asymmetree.treeevolve.species import species_tree_n_age as species_tree_n_age
from asymmetree.treeevolve.species import nonbinary as nonbinary
from asymmetree.treeevolve.genes import dated_gene_tree as dated_gene_tree
from asymmetree.treeevolve.genes import GeneTreeSimulator as GeneTreeSimulator
from asymmetree.treeevolve.genes import prune_losses as prune_losses
from asymmetree.treeevolve.rate_heterogeneity import rate_heterogeneity as rate_heterogeneity
from asymmetree.treeevolve.rate_heterogeneity import (
    autocorrelation_factors as autocorrelation_factors,
)
from asymmetree.treeevolve.rate_heterogeneity import gene_trees as gene_trees
from asymmetree.treeevolve.distance_noise import noisy_matrix as noisy_matrix
from asymmetree.treeevolve.distance_noise import convex_linear_comb as convex_linear_comb
from asymmetree.treeevolve.distance_noise import wrong_topology_matrix as wrong_topology_matrix
