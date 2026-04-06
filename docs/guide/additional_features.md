# Additional features

## Analysis of the simulated scenarios

The subpackage `asymmetree.analysis` contains various functions to extract information from the
simulated scenarios.

### Best matches and orthology

Phylogenetic best matches of a gene x of species X are defined as those genes y of another species Y
that share the lowest common ancestor with x in the gene tree among all genes in that species.
In contrast, two genes are orthologs if their last common ancestor was a speciation event.
 Orthology and reciprocal best matches are closely related.
The module `asymmetree.analysis.best_matches` contains methods for extracting BMGs from trees, least
resolved trees of BMG, and orthology graphs from trees.

Example usage:

```python
from asymmetree.analysis import orthology_from_tree, bmg_from_tree

# tree is a pruned gene tree
orthology_graph = orthology_from_tree(tree)

# best match graph and reciprocal best match graph as a tuple
bmg, rbmg = bmg_from_tree(tree, supply_rbmg=True)
```

### Horizontal gene transfer (HGT)

The module `asymmetree.analysis.horizontal_gene_transfer` contains several functions for the
analysis of horizontal gene transfer events in the simulated scenarios.
In particular, the directed and undirected Fitch graph can be extracted, as well as the pairs of
genes that diverged later than the respective species in which they reside, i.e. the so-called
later-divergence-time (LDT) graph.
The latter situation is indicative for the presence of HGT events in the scenario.

An edge (u, v) in a gene tree is a "true" transfer edge if an HGT event happened on the path from u
to v. In the simulated trees, this is indicated by the attribute `transferred` of `TreeNode` v which
is set to `1`.
An edge (u, v) in the gene tree is an rs-transfer edge (named after the relaxed scenario
framework, Schaller et al. 2021) if the `reconc`s of u and v are incomparable in the corresponding
 species tree `S`.
True and rs-transfer edges may not be equivalent, e.g. when a transfer from branch A to B is
followed by a transfer from B to A and this gene lineage does not survive in branch B.

```python
from asymmetree.analysis import true_transfer_edges, rs_transfer_edges

# gene tree T and species tree S
transfer_edges1 = true_transfer_edges(T)
transfer_edges2 = rs_transfer_edges(T, S)
```

The directed Fitch graph of a tree T has as vertex set the leaves of T and a directed edge (x,y)
when the path from the last common ancestor of x and y to the leaf y contains a transfer edge
(Geiß et al. 2018).
The undirected Fitch graph of a tree T also has as vertex set the leaves of T and an undirected
edge xy when the path from x to y contains a transfer edge (Hellmuth et al. 2018).

```python
from asymmetree.analysis import fitch, undirected_fitch

# gene tree "T", and precomputed "transfer_edges"
dfg = fitch(T, transfer_edges)
ufg = undirected_fitch(T, transfer_edges)
```

As mentioned above, the situation in which two genes diverged later that their corresponding
species witnesses HGT events.
The graph that contains edges for any such gene pair has been termed the later-divergence-time
(LDT) graph.

```python
from asymmetree.analysis import ldt_graph

# gene tree "T", CORRESPONDING species tree "S"
ldt = ldt_graph(T, S, lca_T=None, lca_S=None)
```

In order to reduce runtime, precomputed instances of the class `LCA`
(see [tralda](https://github.com/david-schaller/tralda)) for `T` and `S` can be supplied using the
keyword parameters `lca_T` and `lca_S`, respectively.
Otherwise, such instances are initiated within the function.

## Distance Matrix and Noise

Distances derived from (real-life) gene or protein sequences are burdened with noise.
Such data can either be modeled by simulating sequences, or by disturbing the distances specified
by a given tree directly.
The latter alternative is described briefly in this section.

The additive (i.e., noiseless) distance from a **pruned** gene tree can be computed using the
function `distance_matrix(tree)` from `asymmetree.utils.phylogenetic_trees`.
It returns a tuple containing a list of leaves in the tree (corresponding to the row/column order)
and the distance matrix as a 2-dimensional `numpy` array.

```python
# tree is a pruned gene tree
leaves, D = tree.distance_matrix()
leaf_index_dict = {leaf: i for i, leaf in enumerate(leaves)}
```

In the next step, noise can be introduced into a distance matrix using the
`asymmetree.treeevolve.distance_noise` module.
Random noise can be simulated with the function `noisy_matrix(orig_matrix, sd)`.

The following keyword parameters are available:

| Parameter (with default value) | Description and type |
| --- | ----------- |
| `orig_matrix` | Matrix to be disturbed (2-dimensional `numpy` array). |
| `sd` | Standard deviation of a normal distribution with mean 1 from which noise factors are drawn (`float`). |
| `metric_repair="reject"` | Method to ensure that the resulting distance matrix is still a metric, available are the rejection of noise steps that violate the metric property (`"reject"`), the decrease-only metric repair (`"DOMR"`) and the general metric repair (`"general"`) algorithm. |

Alternatively, the function `convex_linear_comb(D1, D2)` can be used to simulate systematically
biased noise by computing a linear convex combination with a disturbance matrix.
The function thus takes two distance matrices (`numpy` arrays) not necessarily of the same size as
input and disturbs them with one another.
The contribution of the respective disturbance matrix is controlled by the keyword parameter
`alpha` (default is `0.5`).
If the keyword parameter `first_only` is `True`, only the first disturbed matrix is returned.
Otherwise, both are returned in a `tuple`.
