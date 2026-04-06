# Phylogenetic tree simulation

The subpackage `asymmetree.treeevolve` contains modules for the simulation and manipulation of
 species trees and gene trees.
In terms of divergence time, these trees define an ultrametric on the set of their (extant) leaves.
Gene trees, furthermore, can be manipulated with a realistic rate heterogeneity among their branches
resulting in general additive distances (but no longer ultrametric).

A typical simulation consists of the following steps:

* dated species tree (models e.g. "Yule", and "(episodic) birth-death process"; conditioned on the
  number of leaves, the age of the tree, or both)
* dated gene tree(s) (gene duplications, losses, horizontal gene transfer (HGT), and gene
  conversion)
* assignment of asymmetric evolution rates to paralogous genes
* pruned gene tree(s) (removal of all branches that lead to losses only)

The resulting gene trees have edge lengths (`dist`) that correspond to the product of the divergence
time between the respective nodes and the evolutionary rates that were assigned to them.
Such a tree defines a distance matrix on its set of leaves (more precisely, an additive metric).
Noise can be added to this matrix by several methods.
Alternatively, sequences can be simulated along the tree, from which distances can be reestimated.

Example usage:

```python
import asymmetree.treeevolve as te
from asymmetree.utils.phylogenetic_trees import to_newick

# simulate and species tree with 10 leaves
S = te.species_tree_n_age(10, 1.0, contraction_probability=0.2)

# simulate a gene tree along the species tree S
T = te.dated_gene_tree(
    S,
    dupl_rate=0.5,
    loss_rate=0.3,
    hgt_rate=0.1,
    prohibit_extinction="per_species",
)

# simulate evolution rates for the  branches in the gene tree
# and update the distances in T accordingly
T = te.rate_heterogeneity(
    T,
    S,
    base_rate=1,
    autocorr_variance=0.2,
    rate_increase=("gamma", 0.5, 2.2),
)

# remove all gene loss branches
T = te.prune_losses(T)

# print the resulting tree in Newick format and save it to file
print(to_newick(T))
T.serialize("path/to/file.json")
```
