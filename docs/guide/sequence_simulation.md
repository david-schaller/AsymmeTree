# Sequence simulation

The subpackage `asymmetree.seqevolve` contains modules for the simulation of nucleotide or amino
acid sequences along a phylogenetic tree using time-continuous Markov chains, as usually applied
for this purpose (for textbooks, see e.g. Felsenstein 2004, Yang 2006, Yang 2014).

These models typically take a substitution-rate matrix and the equilibrium frequencies of the
states (i.e. the nucleotides or amino acids as input).
Moreover, insertions and deletions (indels) and heterogeneity among the sites can be simulated.

A typical simulation therefore is run with the following components:

* substitution model (mandatory; model for nucleotides e.g. "JC69", "K80", "GTR"; for amino acids
  "DAYHOFF", "BLOSUM62", "JTT", "WAG", "LG")
* indel model (based on the tool _Dawg_ by Cartwright 2005)
* heterogeneity model (constant / sitewise / number of classes; proportion of invariant sites)

Accordingly, the class `Evolver` in `asymmetree.seqevolve.evolver` requires the parameter
`subst_model` for its initialization, and has optional parameter `indel_model` and `het_model`.

Functions for outputting the true alignment are provided for several formats including phylip and
clustal.

Example usage:

```python
from tralda.datastructures import Tree
import asymmetree.seqevolve as se

T = Tree.load("path/to/file.json")

# specify models
subst_model = se.SubstModel("a", "JTT")
indel_model = se.IndelModel(0.01, 0.01, length_distr=("zipf", 1.821))

# initialize Evolver instance
evolver = se.Evolver(subst_model, indel_model=indel_model)

# simulate sequences along the tree
evolver.evolve_along_tree(T, start_length=150)

# print the node labels and the corresponding sequences
for node, sequence in evolver.sequences.items():
    print(node.label, subst_model.to_sequence(sequence))

# construct the true alignment (sequences with gaps) and write it to file
evolver.true_alignment(
    write_to="path/to/alignment.phylip",
    alignment_format="phylip",
)
```

## Substitution model

A substitution model usually comprises an exchangeability matrix S and a vector π containing the
equilibrium frequencies of the alphabet A of nucleobases, amino acids et cetera.
From this, the rate matrix Q can be computed as SΠ where Π=diag(π<sub>1</sub>,...,π<sub>|A|</sub>)
(Yang 2006).
The substitution probability matrix, in turn, is given by P=e<sup>Qt</sup> which is computed
efficiently by AsymmeTree using matrix diagonalization.

The following models are available:

Nucleotide models: `"n"`, amino acid models: `"a"`, codon models are not supported at the moment

| Model | Type | Reference | required parameters (`kwargs`) |
| --- | --- | --- | --- |
| `"JC69"` | `"n"`/`"a"` | Jukes & Cantor 1969 |  |
| `"K80"` | `"n"` | Kimura 1980 | `kappa` (transition/transversion rate ratio) |
| `"GTR"` | `"n"` | General time-reversable model (GTR), Tavare 1986 | `abcdef` (list of rates (a) C&#8596;T, (b) A&#8596;T, (c) G&#8596;T, (d) A&#8596;C, (e) C&#8596;G, (f) A&#8596;G); `f` (list of equilibrium frequencies A/C/G/T) |
| `"DAYHOFF"` | `"a"` | Dayhoff 1978 |  |
| `"BLOSUM62"` | `"a"` | BLOSUM62, Henikoff 1992 |  |
| `"JTT"` | `"a"` | Jones, Taylor & Thornton 1992 |  |
| `"WAG"` | `"a"` | Whelan & Goldman 2001 |  |
| `"LG"` | `"a"` | Le & Gascuel 2008 |  |
| `"CUSTOM"` | `"n"`/`"a"` |  | `filename` (path to a file with a model in PAML format) |

Note that a custom substitution model can be specified via `model_name="CUSTOM"`.
In this case, the path to the model in PAML format (Yang 1997) must be supplied.
Moreover, the model type (`"n"`/`"a"`) must fit this model.

Example usage:

```python
from asymmetree.seqevolve import SubstModel

# non-empirical model
subst_model_jc69 = SubstModel("n", "K80", kappa=2.0)

# empirical model
subst_model_jtt = SubstModel("a", "JTT")

# custom model
subst_model_custom = SubstModel(
    "a", "CUSTOM", filename="path/to/custom.paml"
)
```

## Indel model

Insertions and deletions are modeled as in _Dawg_ (Cartwright 2005).
An indel model requires sitewise rates for insertion `insertion_rate` and deletion `deletion_rate`.

The following keyword parameters are available:

| Parameter (with default value) | Description and type |
| --- | ----------- |
| `length_distr=("zipf", 1.821)` | Distribution of indel length, the default value is a Zipf distribution with empirically observed parameters (Chang 2004).  |
| `min_length=1` | Minimal value (`int`) at which the specified distribution is truncated, must be less than the expected value of the distribution, `None` means no limit.  |
| `max_length=None` | Maximal value (`int`) at which the specified distribution is truncated, must be greater than the expected value of the distribution, `None` means no limit.  |

For available length distributions and their syntax see below.
A zipf or negative binomial distribution are typically used for this purpose
(Cartwright 2005, Dalquen 2012).

Example usage:

```python
from asymmetree.seqevolve import IndelModel

# sitewise insertion rate 0.01 and deletion rate 0.008
indel_model = IndelModel(
    0.01, 0.008, length_distr=("zipf", 2.0)
)
```

## Heterogeneity model

Selective pressure usually varies among the sites of a sequence under evolution.
To model this, rate factors r for single sites or groups of sites are commonly drawn from a Gamma
distribution (+Γ) with mean 1 and parameter α (Cartwright 2005, Fletcher 2009, Dalquen 2012).
The rate matrix rQ is then used instead of Q.
Note that smaller values for α correspond to higher heterogeneity.

AsymmeTree supports two modes of the +Γ-model.
You can specify a number of classes to which the sites are assigned randomly and uniformly
distributed.
Sites of the same class share a common factor r.
The other possibility is a sitewise heterogeneity, i.e., every site has its own rate.
In both cases, the rate or class membership is inherited from the parent sites during the evolution
along a tree.
Note that the sitewise mode is expected to have a longer running time.

An other aspect of among site heterogeneity is the modeling of invariant sites (+I), i.e., sites
that never mutate at all as a result of very strong selective pressure.
The (expected) proportion p of invariant sites can be specified by the user,
and sites are assigned as `"invariant"` with probability p.
Note that p>0 affects the overall substitution rate.
In other words, the rates of the non-invariant sites are **not** adjusted to compensate the
decreased number of expected substitution over all sites.

The following keyword parameters are available:

| Parameter (with default value) | Description and type |
| --- | ----------- |
| `alpha` | Parameter α of the Gamma distribution (`float`). |
| `classes=5` | Number of classes; sites in the same class share the same rate factor (`int`). |
| `sitewise=False` | If `True`, factors are drawn sitewise; the number of classes is ignored in this case (`bool`). |
| `invariant=0.0` | (Expected) proportion of invariant sites (`float`). |

Note that the +I-model can be used without the +Γ-model by setting `classes=1` (the single class
will have a factor of 1) and `"invariant"` to some proportion greater than 0.

Example usage:

```python
from asymmetree.seqevolve import HetModel

# +G
het_model_g  = HetModel(1.5, classes=10)

# +G +I
het_model_gi = HetModel(2.0, sitewise=True, invariant=0.2)

# +I (value of alpha is irrelevant)
het_model_i  = HetModel(1.0, classes=1, invariant=0.25)
```

## The class Evolver

The class `Evolver` in the module `asymmetree.seqevolve.evolver` evolves a sequence according to
the specified models (see previous sections) along a phylogenetic tree.
In AsymmeTree, the `dist` attribute of the vertex v in an edge uv of the tree (u is closer to the
root) is always used as the **expected number of substitutions** along this edge.
Thus, PAM distances as e.g. used optionally in _ALF_ (Dalquen 2012) are not supported.
The `dist` attribute is also used as the duration of the Markov process in which insertions and
deletions are drawn.

The following parameters are available for initialization:

| Parameter (with default value) | Description and type |
| --- | ----------- |
| `subst_model` | Substitution model (`SubstModel`) |
| `indel_model=None` | Model for insertions and deletions (`IndelModel`). |
| `het_model=None` | Model for among site heterogeneity and invariant sites (`HetModel`). |
| `gillespie="auto"` | If `True`, the Gillespie algorithm is used for the substitution process instead of the computation of P=e<sup>Qt</sup> (`bool`); the default is "auto", in which case the exchange probability matrix is used except if rate heterogeneity is enabled. |

Once the `Evolver` is initialized, its function `evolve_along_tree(tree)` can be called to evolve a sequence along a tree.

The following parameters are available for this function:

| Parameter (with default value) | Description and type |
| --- | ----------- |
| `tree` | Phylogenetic tree with `dist` attributes set (`Tree`). |
| `start_length=200` | Length of the root sequence which is randomly drawn from the equilibrium frequencies in the specified substitution model (`int`). |
| `start_seq=None` | Root sequence; must be compatible with the specified substitution model (`model_type="n"`/`"a"`); if supplied, the `start_length` attribute is ignored (`str`). |

The sequences of a simulation run are returned by this function (and also available via the
attribute `sequences` as long as the function has not been called again), which is a `dict`
containing the nodes (inner and leaf nodes) as keys and instances of type `EvoSeq` as values.
The latter can be converted into `str` using `subst_model.to_sequence(evo_seq)`.
In addition, the `str` representations of the simulated sequences can be obtained and written to
file (in fasta format) using the functions `get_sequences()` and `write_sequences()`, resp., of
an `Evolver` instance.

The function `true_alignment()` can be used to compute (and optionally write into a file) the
"true" alignment of a simulation run.

The following parameters are available for this function:

| Parameter (with default value) | Description and type |
| --- | ----------- |
| `include_inner=True` | If `True`, include also inner nodes in the alignment; otherwise only sequences of leaf nodes are aligned (`bool`). |
| `write_to=None` | Path and filename for the output (`str`). |
| `alignment_format="phylip"` | Format of the alignment file; available are `"phylip"`, `"clustal"`, and `"pretty"` (`str`). |

Example usage:

```python
from asymmetree.seqevolve import Evolver, SubstModel, IndelModel, HetModel

# specify models (only SubstModel is mandatory)
subst_model = SubstModel("a", "WAG")
indel_model = IndelModel(0.01, 0.01, length_distr=("zipf", 1.821))
het_model  = HetModel(2.0, classes=7, invariant=0.1)

# evolve sequences along a phylogenetic tree
evolver = Evolver(
    subst_model, indel_model=indel_model, het_model=het_model
)
evolver.evolve_along_tree(tree, start_length=150)

# print resulting sequences
for node, evo_seq in evolver.sequences.items():
    print(node.label, subst_model.to_sequence(evo_seq))

# output true alignment
evolver.true_alignment(
    write_to="path/to/alignment.phylip",
    alignment_format="phylip",
)
```
