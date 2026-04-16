# Preliminaries

## Tree data structures

The class `Tree` in [tralda](https://github.com/david-schaller/tralda) implements a tree data
structure which is essential for most of the modules in the package.

The nodes in the simulated trees have attributes:

| Attribute | Description and type |
| --- | ----------- |
| `label` | Node label, usually `int` |
| `event` | The type of the event, in gene trees: "S" for speciation, "D" for duplication, "H" for horizontal gene transfer, "L" for loss, "GC" for gene conversion, `str` |
| `reconc` | Only gene trees; species in which the gene resides, i.e., `label` of some vertex in a species tree, `int` for extant genes, can be of type `tuple` (of two `int`s) for inner and loss vertices |
| `tstamp` | Time stamp of the event (`double`) |
| `dist` | Evolutionary distance from the parent vertex (`double`); if no evolution rates (see below) were simulated yet, then this value corresponds to the divergence time between the vertex and its parent |
| `transferred` | Only gene trees; indicates whether the edge from the parent is the transfer edge from an HGT event; `1` if this is the case and `0` otherwise |


All simulated trees can be serialized in either JSON format or the Python-specific serialization
format (using the library `pickle`).
By default, the serialization format is inferred from the file extension.
Alternatively, it can be specified as keyword argument, e.g. `mode="json"`.

Example usage:

```python
from tralda.datastructures import Tree

# T is of type "Tree"
T.serialize("path/to/tree.pickle")
# or
T.serialize("path/to/tree.json")

# reload serialized tree
T_from_file = Tree.load("path/to/tree.json", mode="json")
```

Moreover, the module `asymmetree.utils.phylogenetic_trees` provides the functions `to_newick()` and
`parse_newick()` for converting and parsing trees to and from Newick format, respectively.
In case of a gene tree, the reconciliation is represented in brackets, e.g.,
`(3<1>:0.534,2<2>:0.762)1<0>:0.273`.
To suppress this, use `to_newick(reconc=False)`.
Similarly, to suppress the distances, you can use `to_newick(distance=False)`.
The function `parse_newick()` can handle this customized format as well as the standard Newick
format.


## Distributions for sampling

The following distributions are available for sampling:

| Distribution | Syntax | Parameters |
| --- | --- | --- |
| constant | `x` or `("constant", x)` | `x` must be a number |
| uniform (continuous) | `("uniform", a, b)` | `a`<=`b` must be numbers |
| uniform (discrete) | `("discrete_uniform", a, b)` | `a`<=`b` must be `int`s |
| gamma | `("gamma", shape, scale)` | `shape` and `scale` must be `floats`>0 |
| gamma (mean) | `("gamma_mean", mean)` | `mean` must be a number>0, shape=1 and scale=`mean`/shape |
| exponential | `("exponential", rate)` | `rate` must be a `float`>=0 |
| Zipf | `("zipf", a)` | `x` must be a number |
| negative binomial | `("negative_binomial", r, q)` | `r`>=1 must be an `int`>=1, 0<`q`<1 must be a `float` |
