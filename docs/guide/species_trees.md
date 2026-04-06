# Species Tree Simulation

AsymmeTree implements functions to simulate species trees conditioned on

* the number of (extant-species) leaves `n` (function `species_tree_n(n, **kwargs)`)
* the time span `age` covered by the tree (function `species_tree_age(age, **kwargs)`), or
* both `n` and `age` (function `species_tree_n_age(n, age, **kwargs)`)

under a pure-birth Yule model (Yule 1924), a constant-rate birth-death process (Kendall 1948,
Hagen & Stadler 2018), or an episodic birth-death process (Stadler 2011).
The latter is not available for conditioning on both `n` and `age`.
The model is specified by setting the `model` parameter of the respective function to `"yule"`
 (default), `"BDP"`, or `"EBDP"`.
Note that when conditioning to `n` alone, the time covered by the resulting trees varies between
the simulations. Conversely, when conditioning to `age`, the number of leaves varies.

All simulated trees are by default "planted", i.e., the root has a single child and the edge to
this child represents the ancestral lineage.
For any model, the root of the resulting tree has the maximal time stamp and all (extant) species
have time stamp 0.0.

The following keyword parameters are available for conditioning on `n`:

| Parameter (with default value) | Description |
| --- | ----------- |
| `model="yule"` | Model for the species tree simulation. |
| `innovation=False` | If True, use the innovation model (Keller-Schmidt & Klemm 2012) to sample a lineage for the next speciation event; only available for the Yule model; the default is False, in which case the lineage is chosen uniformly at random among the currently existing lineages. |
| `birthrate=None` | Birth rate, only relevant for models `"yule"` and `"BDP"`. The default is None, in which case the birth rate is set to 1.0 unless `episodes` are specified. |
| `deathrate=None` | Death rate, only relevant for model `"BDP"`. The default is None, in which case the death rate is set to 0.0 unless `episodes` are specified. |
| `episodes=None` | Episodes for episodic birth-death process, only relevant for `"EBDP"`, the episodes of the `"EBDP"` model must be supplied as a list of tuples/lists where each episode has the structure `(birthrate, deathrate, proportion_of_survivors, time_stamp)`, the first elements in this list correspond to the most recent ones, i.e., the first episode should have a time stamp of 0.0. |
| `planted=True` | Add a planted root that has the first true speciation node as its single neighbor, this way duplication (and loss) events can occur before the first speciation event in a subsequent gene tree simulation. |
| `remove_extinct=False` | Remove all branches leading to losses, only relevant for models with death events. |
| `contraction_probability=0.0` | Probability that an inner edge is contracted; the default is 0.0, in which case the tree is binary; only one of this parameter and `contraction_proportion` may be non-zero. |
| `contraction_proportion=0.0` | The proportion of inner edges to be contracted.; the default is 0.0, in which case the tree is binary; only one of this parameter and `contraction_probability` may be non-zero. |
| `contraction_bias=False` | Specifies whether shorter edges, i.e., with a smaller difference t of the time stamps, have a higher probability to be contracted; only relevant if `contraction_proportion > 0.0`; the default is False, in which case all edges have the same probability to be contracted, the options `"inverse"` and `"exponential"` mean that an edge is sampled weighted by 1/(a * t) or e^(-a * t), respectively, where a is a user-defined factor (`bias_strength`). |
| `bias_strength=1.0` | Intensity factor for preferring shorter edges to be contracted. |

The available parameters for the functions `species_tree_age(age, **kwargs)`) and
`species_tree_n_age(n, age, **kwargs)` are largely the same, but the available options differ
slightly.
See also the [API reference](../api/treeevolve.md) for details.

Example usage:

```python
import asymmetree.treeevolve as te

S1 = te.species_tree_n_age(10, 1.0, contraction_probability=0.2)
print(S1.to_newick())

S2 = te.species_tree_n(
    10,  # number of (extant) species leaves
    model="EBDP",  # episodic birth-death process
    episodes=[
        # birthrate, deathrate, proportion of survivors, time stamp
        (1.0, 0.3, 0.8, 0.0),
        (0.9, 0.4, 0.6, 0.3),
    ],
)
print(S2.to_newick())
```
