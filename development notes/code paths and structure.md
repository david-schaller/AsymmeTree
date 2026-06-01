This report is inside of the repository of Asymmetree, in a development branch. The main directory is `../`, there we can find the [main readme](../README.md) which includes links for [wiki](https://github.com/david-schaller/AsymmeTree/wiki/Manual) and [documentation](https://david-schaller.github.io/docs/asymmetree/). The following explanation is taken from the wiki.

Also see  [skill contributing.md](skill contributing.md)  and  [dependency_diagram.html](Example code/dependency_diagram.html) 

# AsymmeTree framework

AssymmeTree is a powerful tool for simularion of two-level evolutionary scenarios. The following text taken from the [paper of asymmetree](https://www.mdpi.com/2674-113X/1/3/13) expalins their implementation. For us, the most intersting part are steps (1) and (2).

Like many other tools for phylogenetic simulation [[25](https://www.mdpi.com/2674-113X/1/3/13#B25-software-01-00013),[26](https://www.mdpi.com/2674-113X/1/3/13#B26-software-01-00013),[40](https://www.mdpi.com/2674-113X/1/3/13#B40-software-01-00013),[48](https://www.mdpi.com/2674-113X/1/3/13#B48-software-01-00013)], `AsymmeTree` generates complex gene family histories in a multilevel process that consists of the following five steps. **(1)** A dated species tree S is simulated whose leaves correspond to extant species (thus having  time stamp 0) and possibly (depending on simulation models and  parameters) extinction events. **(2)** Along this species tree, one or  multiple gene trees T are simulated using a variant of a constant-rate birth-death process [[27](https://www.mdpi.com/2674-113X/1/3/13#B27-software-01-00013)] that considers speciations, gene duplications, and (additive) HGTs as  “birth events” and gene losses as “death events”. In addition, replacing HGT and gene conversion lead to a bifurcation in one gene lineage  while, at the same time, another lineage gets lost. **(3)** In a next step,  evolution rate heterogeneity between the branches is introduced in a  manner that accounts for both species effects as well as asymmetric  evolution of paralogs after gene duplication. **(4)** The pruned gene tree  is then obtained by removing all branches that lead to loss events only  and suppressing all vertices that only have a single child left. In a final step **(5)**, nucleotide or amino acid  sequences can be evolved along the gene tree using a continuous-time  Markov chain, which is the standard model for this purpose [[31](https://www.mdpi.com/2674-113X/1/3/13#B31-software-01-00013),[32](https://www.mdpi.com/2674-113X/1/3/13#B32-software-01-00013)]. `AsymmeTree` supports a variety of models for substitution, indels, and among-site  heterogeneity.

[...]

Gene trees are simulated along a species tree using a constant-rate birth-death process [[27](https://www.mdpi.com/2674-113X/1/3/13#B27-software-01-00013)] where rates for the four types of events duplication, loss, HGT, and  gene conversion are user-defined. To this end, the user must specify  rates for the four event types which serve as parameters for exponential distributions from which waiting times until the next events are drawn. By setting the respective rate to zero, an event type is disabled  completely. The simulation starts with a single gene in the root of the species tree and proceeds stepwise in time towards the leaves by drawing waiting  times until the next event and lineages in which these events take place. At each point in the simulation, the total rate is given by the  sum of the four event types over the currently existing lineages in the  gene tree under construction. Hence, the simulation of all gene lineages progresses synchronously. This avoids difficulties arising by the fact  that some processes, such as replacing HGT as introduced below,  introduce dependencies between the gene lineages. In particular, with  this method, it is not necessary that simulated branches have to be  invalidated as a consequence of an event in a lineage that is processed  later in the simulation as it, e.g., occurs in [[26](https://www.mdpi.com/2674-113X/1/3/13#B26-software-01-00013)].

In addition to the randomly-occurring events, speciations and species  extinctions (which are determined by the species tree and therefore  fixed) are included as branching and loss events, respectively, and  affect all currently existing lineages in the respective species  lineage. In case of a speciation, each offspring species receives one  copy of each original gene lineage. On the other hand, the extinction of a species leads to loss of all of its gene lineages. If a waiting time  for a duplication, loss, or HGT event is drawn such that the next event  in the species tree occurs earlier, then this waiting time is discarded, the time is updated to the next species tree event, and the latter is  executed.

## Current structure of AsymmeTree code

This report is inside of the repository of Asymmetree, in a development branch. The main directory is `../`, there we can find the [main readme](../README.md) which includes links for [wiki](https://github.com/david-schaller/AsymmeTree/wiki/Manual) and [documentation](https://david-schaller.github.io/docs/asymmetree/). The following explanation is taken from the wiki.

The module `asymmetree.genome.GenomeSimulation` provides functions that combine the simulation of phylogenetic trees and sequences. In particular, the class `GenomeSimulator` combines multiple steps in order to conveniently simulate whole genomes/proteomes. The (optional) output directory contains serialized trees, fasta files, and the true alignments. The gene trees and the sequences are simulated in subsequent steps using the classes' functions (i) `simulate_gene_trees(n, **kwargs)` and (ii) `simulate_sequences(subst_model, **kwargs)`.

The first step (i) takes the same keyword parameters as input as the function `gene_trees()` where `n` is the number of gene families to be simulated. Thus, rates for the three event types (`dupl_rate`, `loss_rate`, `hgt_rate`), autocorrelation (`autocorr_variance`), the distribution of base rates (`base_rate_distr`) etc. can be specified.

The second step (ii) simulates the sequences along the pruned part (without loss branches) of the simulated gene trees.

After step (i), the `list`s of full and pruned gene trees are accessible via the attributes `true_gene_trees` and `pruned_gene_trees`, respectively. Moreover, the full gene trees are serialized into the directory 'true_gene_trees' if an output directory was specified. After step (ii), the `list`s of sequence `dict`s are accessible via the attribute `sequence_dicts`.

Example usage:

```python
from asymmetree.treeevolve import species_tree_n_age
from asymmetree.genome import GenomeSimulator
from asymmetree.seqevolve import SubstModel, IndelModel

# simulate the common species tree
S = species_tree_n_age(10, 1.0, model='yule')

# specify models for sequence evolution
subst_model = SubstModel('a', 'JTT')
indel_model = IndelModel(0.01, 0.01, length_distr=('zipf', 1.821))

# initialy GenomeSimulator instance
gs = GenomeSimulator(S, outdir='simulation_directory')

# simulate 50 gene trees along the species tree S (and write them to file)
gs.simulate_gene_trees(50, dupl_rate=1.0, loss_rate=0.5,
                       base_rate=('gamma', 1.0, 1.0),
                       prohibit_extinction='per_species')

# simulate sequences along the gene trees
gs.simulate_sequences(subst_model,
                      indel_model=indel_model,
                      het_model=None,
                      length_distr=('constant', 200))

# results have been written to directories 'true_gene_trees',
# 'fasta_files' and 'alignments' in 'path/to/genome_directory'
```

