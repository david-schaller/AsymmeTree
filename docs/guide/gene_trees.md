# Gene tree simulation

## Birth-death process

Dated gene trees are simulated along a given species tree `S` using a birth-death process
(Kendall 1948, Gillespie 1976) with speciation events as additional branching events
(fixed time points given by the species tree).

At each time point, the total event rate is given by the sum of the event rates over all branches
that are currently active (not extinct).
Thus, the total event rate in general increases during the simulation if the loss rate does not
dominate the rates of the branching events.

To simulate gene tree, use the class `GeneTreeSimulator` or the function
`dated_gene_tree(S, **kwargs)` with a species tree of type `Tree`.

The following keyword parameters are available:

| Parameter (with default value) | Description and type |
| --- | ----------- |
| `dupl_rate=0.0` | Duplication rate (`float`) |
| `loss_rate=0.0` | Loss rate (`float`) |
| `hgt_rate=0.0` | Horizontal gene transfer rate (`float`) |
| `gc_rate=0.0` | Gene conversion rate (`float`) |
| `dupl_polytomy=0.0` | Allows non-binary duplication events by drawing a number from a Poisson distribution with rate parameter `dupl_polytomy` (copy number = 2 + drawn number) (`float`) |
| `prohibit_extinction="per_species"` | Avoid loss events for genes that are the last survivor in their species branch (`"per_species"`), the last survivor of the whole family (`"per_family"`); or no constraints  (`False`) |
| `replace_prob=0.0` | Probability that an HGT event is replacing (rather than additive), i.e., it replaces a homolog in the receiving species branch; default is `0.0` in which case all HGT events are additive (`float`) |
| `additive_transfer_distance_bias=False` | Specifies whether closer related species have a higher probability to be the recipient species in an additive HGT event. The default is False, in which case the recipient species is chosen at random among the co-existing species. The options `"inverse"` and `"exponential"` mean that a species branch is sampled weighted by 1/(a * t) or e^(-a * t), resp., where t is the elapsed time between the last common ancestor of the two species branches and the time of the event and a is a user-defined factor (see `transfer_distance_bias_strength` below) |
| `replacing_transfer_distance_bias=False` | Specifies whether closer related gene branches have a higher probability to be replaced in a replacing HGT event. The default is False, in which case the replaced gene is chosen at random among the co-existing gene branches. The options `"inverse"` and `"exponential"` mean that a species branch is sampled weighted by 1/(a * t) or e^(-a * t), resp., where t is the elapsed time between the last common ancestor of the two gene branches and the time of the event and a is a user-defined factor (see `transfer_distance_bias_strength` below) |
| `transfer_distance_bias=False` | Sets a common bias mode for additive and replacing HGT, see description of parameters `additive_transfer_distance_bias` and `replacing_transfer_distance_bias`. If the latter are no set to the default (False), then these optioned are prioritized. |
| `transfer_distance_bias_strength=1.0` | Intensity of the transfer distance bias (factor a) for additive and replacing HGT |
| `gc_distance_bias=False` | Specifies whether closer related gene branches have a higher probability to be replaced in a gene conversion event; the default is False, in which case the replaced gene is chosen at random among the paralogs in the respective species lineage. The options `"inverse"` and `"exponential"` mean that a paralog is sampled weighted by 1/(a * t) or e^(-(a * t)), resp., where t is the elapsed time between the last common ancestor of the two gene branches and the time of the event, see [1], and a is a user-defined factor |
| `gc_distance_bias_strength=1.0` | Intensity of the distance bias (factor a) for gene conversion |

For the constraints to avoid extinction, the loss rates in the respective branches are temporarily
set to zero.
Note that if replacing HGT is enabled (`replace_prob` >0.0), then the resulting tree may contain
loss leaves even if the loss rate is zero since replacement is modeled by a loss of a paralog.
Hence, constructing the pruned tree (see below) also becomes relevant in the latter setting.

AsymmeTree supports transfer distance bias for horizontal gene transfers similar as the tool
SaGePhy (Kundu and Bansal 2019).
If this bias is enabled, then

- in case of an **additive transfer** (no homolog in the recipient species is replaced), the time
  elapsed since the last common ancestor of the **species** branches and the HGT event is used for
  computing the weights, whereas
- in case of a **replacing transfer**, the time elapsed since the last common ancestor of the
  respective **gene** branches and the HGT event is used.

## Pruning gene trees

The function `prune_losses(tree)` returns the observable part of a gene tree, i.e., it copies the
tree, removes all branches that lead to loss events only and suppresses all inner nodes with only
one child.
It also removes the planted root.

Example usage:

```python
import asymmetree.treeevolve as te

# S is a species tree of type Tree
tree_simulator = te.GeneTreeSimulator(S)
tree = tree_simulator.simulate(dupl_rate=1.0, loss_rate=0.5, hgt_rate=0.1)

# or
tree = te.dated_gene_tree(
    S,
    dupl_rate=1.0,
    loss_rate=0.5,
    hgt_rate=0.1,
    replace_prob=0.5,
)

# prune all loss branches and the planted root
observable_gene_tree = te.prune_losses(tree)
```

## Rate Heterogeneity

The subpackage `treeevolve` contains functions to model realistic (asymmetric) evolution rates for
a given gene tree.
Moreover, correlation of the evolution rates between genes of the same (and closely related)
species is introduced (autocorrelation, Kishino 2001).
The function `rate_heterogeneity(T, S)` takes a gene tree `T` and the **corresponding** species
tree `S` as input, and manipulates the branch length of the gene tree.

The following keyword parameters are available:

| Parameter (with default value) | Description and type |
| --- | ----------- |
| `base_rate=1.0` | Starting value for the substitution rate (per time unit) and expected value for conserved genes (`float`) |
| `autocorr_factors=None` | A `dict` containing rate factors for the edges of the species tree |
| `autocorr_variance=0.0` | Variance factor for a lognormal distribution that controls autocorrelation between genes of the same (and closely related) species, the higher the lower the autocorrelation; only relevant if `autocorr_factors` is not directly supplied (`float`) |
| `rate_increase=("gamma", 0.5, 2.2)` | Distribution of the (relative) rate increase (w.r.t. the base rate) for divergent genes, i.e. to a factor 1 + x, the parameters the for default Gamma distribution are chosen to fit observed asymmetries between paralogs in yeast data (Byrne et al. 2007), (`tuple`) |
| `CSN_weights=(1, 1, 1)` | Weights for choice between conservation, subfunctionalization, and neofunctionalization after a duplication event (`tuple` of three `int`/`float`s) |
| `inplace=True` | Manipulate edge lengths (`dist`) of the gene tree in-place, otherwise copy the tree |

It is recommended to apply the rate assignment to the true gene tree that still contains loss
events.
Note that the rates are used to manipulate the `dist` attributes in the gene tree and not returned
explicitly.

The function `gene_trees(S)` combines the simulation of dated gene trees and the rate assignment
into one step.
It returns a `list` of gene trees that shared the same rate factors for the branches in the species
tree (autocorrelation factors) in the rate assignment procedure.
Moreover, a distribution for the base rate (assigned the planted edge of the gene tree) and for the
event rates can be specified with the parameters `base_rate`, `dupl_rate`, `loss_rate`, and
`hgt_rate`.
For available distributions and their syntax see the example below.

Example usage:

```python
import asymmetree.treeevolve as te

S = te.species_tree_n_age(10, 1.0)

# true gene tree (with losses)
tree = te.dated_gene_tree(S, dupl_rate=1.0, loss_rate=0.5, hgt_rate=0.1)
te.rate_heterogeneity(
    tree,  # gene tree
    S,  # corresponding species tree
    base_rate=1,  # base rate for conserved genes
    autocorr_variance=0.2,  # variance for autocorrelation
    rate_increase=("gamma", 0.5, 2.2),  # distribution for rate increase for divergent genes
)

# or in a single step
trees = gene_trees(
    S,  # species tree
    n=10,  # number of gene trees to simulate
    # use a distribution or a constant value for the rates, e.g.:
    dupl_rate=("gamma", 1.0, 0.8),  
    loss_rate=0.5,
    hgt_rate=0.1,
    base_rate=("constant", 1.0),
    autocorr_variance=0.2,  # variance for autocorrelation (no distribution!)
    rate_increase=("gamma", 0.5, 2.2),  # distribution for rate increase for divergent genes
)

# prune the simulated gene trees
pruned_trees = [te.prune_losses(tree) for tree in trees]
```
