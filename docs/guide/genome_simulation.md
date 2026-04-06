# Genome simulation

The subpackage `asymmetree.genome` provides functions that combine the simulation of
phylogenetic trees and sequences.
In particular, the class `GenomeSimulator` combines multiple steps described in the previous
sections in order to conveniently simulate whole genomes/proteomes.
The (optional) output directory contains serialized trees, fasta files, and the true alignments.
The gene trees and the sequences are simulated in subsequent steps using the classes" functions
(i) `simulate_gene_trees(n, **kwargs)` and (ii) `simulate_sequences(subst_model, **kwargs)`.

The first step (i) takes the same keyword parameters as input as the function `gene_trees()` where
`n` is the number of gene families to be simulated.
Thus, rates for the three event types (`dupl_rate`, `loss_rate`, `hgt_rate`), autocorrelation
(`autocorr_variance`), the distribution of base rates (`base_rate_distr`) etc. can be specified.

The second step (ii) simulates the sequences along the pruned part (without loss branches) of the
simulated gene trees.

The function takes the following parameters as input:

| Parameter (with default value) | Description and type |
| --- | ----------- |
| `subst_model` | Substitution model (`SubstModel`) |
| `indel_model=None` | Model for insertions and deletions (`IndelModel`) |
| `het_model=None` | Model for among site heterogeneity and invariant sites (`HetModel`) |
| `root_genome=None` | `list` of sequences for the roots of the gene trees; must contain the same number of `str` sequences as trees were simulated in step (i) i.e. `n`; sequences must be compatible with the specified substitution model (`model_type="n"`/`"a"`) |
| `length_distr=("constant", 200)` | Distribution of the length of the root sequences if `root_genome` is not supplied; see below |
| `min_length=1` | Minimal value (`int`) at which the specified distribution is truncated, must be less than the expected value of the distribution, `None` means no limit  |
| `max_length=None` | Maximal value (`int`) at which the specified distribution is truncated, must be greater than the expected value of the distribution, `None` means no limit  |
| `write_fastas=True` | If `True` and an output directory was specified, write the `sequences` (one file **per species**) into the directory "fasta_files" in the output directory (`bool`) |
| `write_alignments=True` | If `True` and an output directory was specified, write the true alignments (one file **per gene tree**) into the directory "alignments" in the output directory (`bool`) |

After step (i), the `list`s of full and pruned gene trees are accessible via the attributes
 `true_gene_trees` and `pruned_gene_trees`, respectively.
Moreover, the full gene trees are serialized into the directory "true_gene_trees" if an output
directory was specified.
After step (ii), the `list`s of sequence `dict`s are accessible via the attribute `sequence_dicts`.

Example usage:

```python
from asymmetree.treeevolve import species_tree_n_age
from asymmetree.genome import GenomeSimulator
from asymmetree.seqevolve import SubstModel, IndelModel

# simulate the common species tree
S = species_tree_n_age(10, 1.0, model="yule")

# specify models for sequence evolution
subst_model = SubstModel("a", "JTT")
indel_model = IndelModel(0.01, 0.01, length_distr=("zipf", 1.821))

# initially GenomeSimulator instance
gs = GenomeSimulator(S, outdir="simulation_directory")

# simulate 50 gene trees along the species tree S (and write them to file)
gs.simulate_gene_trees(
    50,
    dupl_rate=1.0,
    loss_rate=0.5,
    base_rate=("gamma", 1.0, 1.0),
    prohibit_extinction="per_species",
)

# simulate sequences along the gene trees
gs.simulate_sequences(
    subst_model,
    indel_model=indel_model,
    het_model=None,
    length_distr=("constant", 200),
)

# results have been written to directories "true_gene_trees",
# "fasta_files" and "alignments" in "path/to/genome_directory"
```
