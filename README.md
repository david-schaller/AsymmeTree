# AsymmeTree
![Logo](manual/images/logo.png)

[![license: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![pypi version](https://img.shields.io/badge/pypi-v2.1.0-blue.svg)](https://pypi.org/project/asymmetree/)

AsymmeTree is an open-source Python library for the simulation and analysis of phylogenetic scenarios.
It includes a simulator for species and gene trees with heterogeneous evolution rates, nucleotide and amino acid sequences with or without indels, as well as whole genomes/proteomes.

Moreover, it includes tools for the extraction of information from the simulated scenarios such as orthology, best matches, and xenology, and a method to estimate rooted species trees from an ensemble of orthology/paralogy relations.

The library is primarily designed to explore and validate mathematical concepts, and to test inference methods for various steps on the way to more realistically-available data, i.e., dated gene trees, additive distances of gene sets, noisy distances and finally sequences.


## Installation

AsymmeTree requires Python 3.7 or higher.

#### Easy Installation with pip

The `asymmetree` package is available on PyPI:

    pip install asymmetree

For details about how to install Python packages see [here](https://packaging.python.org/tutorials/installing-packages/).
    
#### Dependencies

AsymmeTree has several dependencies (which are installed automatically when using `pip`):
* [NetworkX](https://networkx.github.io/)
* [Scipy and Numpy](http://www.scipy.org/install.html)
* [Matplotlib](https://matplotlib.org/)
* [tralda](https://github.com/david-schaller/tralda)

The package builds on top of the `tralda` library for tree algorithms and datastructures. In particular, all simulated trees can be serialized in either JSON format or the Python-specific serialization format (using the library `pickle`).

## Usage and Description

For a more detailed description of the usage and the implementation of the simulator please read the [manual](https://github.com/david-schaller/AsymmeTree/blob/master/manual/AsymmeTreeManual.pdf).

### Simulation of Phylogenetic Trees

The subpackage `treeevolve` contains modules for the simulation and manipulation of species trees and gene trees.

A typical simulation consists of the following steps:
* dated species tree (models e.g. 'innovation', 'Yule' and '(episodic) birth-death process')
* dated gene tree(s) (birth-death process with speciations as additional branching events)
* assignment of asymmetric evolution rates to paralogous genes
* observable gene tree(s) (removal of all branches that lead to losses only)

The resulting gene trees have edge lengths (`dist`) that correspond to the product of the divergence time between the respective nodes and the evolutionary rates that were assigned to them.
Such a tree defines a distance matrix on its set of leaves (more precisely, an additive metric).
Noise can be added to this matrix by several methods.
Alternatively, sequences can be simulated along the tree, from which distances can be reestimated.

<details>
<summary>Example usage: (Click to expand)</summary>
    
    import asymmetree.treeevolve as te
    
    # simulate and species tree with 10 leaves
    S = te.simulate_species_tree(10, planted=True, non_binary_prob=0.2)
    
    # simulate a gene tree along the species tree S
    T = simulate_dated_gene_tree(S, dupl_rate=D, loss_rate=L, hgt_rate=H,
                                   prohibit_extinction='per_species')
    
    # simulate evolution rates for the  branches in the gene tree
    # and update the distances in T accordingly
    T = te.assign_rates(T, S, base_rate=1,
                        autocorr_variance=0.2,
                        rate_increase=('gamma', 0.5, 2.2))
    
    # remove all gene loss branches
    T = te.observable_tree(TGT)
    
    # print the resulting tree in Newick format and save it to file
    print(to_newick(T))
    T.serialize('path/to/file.json')

</details>

### Simulation of Sequences

The subpackage `seqevolve` contains modules for the simulation of nucleotide or amino acid sequences along a phylogenetic tree.
The substitution of sites is modeled by continuous-time Markov chains.
These models typically take a substitution-rate matrix and the equilibrium frequencies of the states (i.e. the nucleotides or amino acids as input).
Moreover, insertions and deletions (indels) and heterogeneity among the sites can be simulated.

A typical simulation therefore is run with the following components (only the substitution model is mandatory):
* substitution model (model for nucleotidess e.g. 'JC69', 'K80', 'GTR'; for amino acids 'DAYHOFF', 'BLOSUM62', 'JTT', 'WAG', 'LG')
* indel model (based on the tool 'Dawg' by Cartwright 2005)
* heterogeneity model (constant / sitewise / number of classes; proportion of invariant sites)

Functions for outputting the true alignment are provided for several formats incl. phylip and clustal.

<details>
<summary>Example usage: (Click to expand)</summary>
    
    from tralda.datastructures import Tree
    import asymmetree.seqevolve as se
    
    T = Tree.load('path/to/file.json')
    
    # specify models
    subst_model = se.SubstModel('a', 'JTT')
    indel_model = se.IndelModel(0.01, 0.01, length_distr=('zipf', 1.821))
    
    # initialize Evolver instance
    evolver = se.Evolver(subst_model, indel_model=indel_model)
    
    # simulate sequences along the tree
    evolver.evolve_along_tree(T, start_length=150)
    
    # print the node labels and the corresponding sequences
    for node, sequence in evolver.sequences.items():
        print(node.label, subst_model.to_sequence(sequence))
    
    # construct the true alignment (sequences with gaps) and write it to file
    evolver.true_alignment(write_to='path/to/aligment.phylip',
                           alignment_format='phylip')
</details>

### Simulation of Genomes

The module `genome.GenomeSimulation` provides functions that combine the simulation of phylogenetic trees and sequences.
This way, whole genomes/proteomes can be simulated conveniently.
The (optional) output directory contains serialized trees, fasta files and the true alignments.

<details>
<summary>Example usage: (Click to expand)</summary>
    
    from asymmetree.treeevolve import simulate_species_tree
    from asymmetree.genome import GenomeSimulator
    from asymmetree.seqevolve import SubstModel, IndelModel
    
    # simulate the common species tree
    S = simulate_species_tree(10, model='innovation')
    
    # specify models for sequence evolution
    subst_model = SubstModel('a', 'JTT')
    indel_model = IndelModel(0.01, 0.01, length_distr=('zipf', 1.821))
    
    # initialy GenomeSimulator instance
    genome_sim = GenomeSimulator(S, outdir='simulation_directory')
    
    # simulate 50 gene trees along the species tree S (and write them to file)
    genome_sim.simulate_gene_trees(50, dupl_rate=1.0, loss_rate=0.5,
                                   base_rate=('gamma', 1.0, 1.0),
                                   prohibit_extinction='per_species')
    
    # simulate sequences along the gene trees (and write them to file)
    genome_sim.simulate_sequences(subst_model,
                                  indel_model=indel_model,
                                  het_model=None,
                                  length_distr=('constant', 200))
</details>

### Analysis of the Simulated Scenarios

The subpackage `analysis` contains various functions to extract information from the simulated scenarios.

Phylogenetic best matches of a gene x of species X are defined as those genes y of another species Y that share the lowest common ancestor with x in the gene tree among all genes in that species. In contrast, two genes are orthologs if their last common ancestor was a speciation event. Orthology and reciprocal best matches are closely related.
The module `analysis.BestMatches` contains methods for extracting BMGs from trees, least resolved trees of BMG, and orthology graphs from trees

The module `analysis.HGT` contains several functions for the analysis of horizontal gene transfer events in the simulated scenarios. In particular, the directed and undirected Fitch graph can be extracted, as well as the pairs of genes that diverged later than the respective species in which they reside, i.e. the so-called later-divergence-time (LDT) graph. The latter situation is indicative for the presence of HGT events in the scenario.

### Species Trees From Orthology

The subpackage `paraphylo` contains a method to compute rooted species tree from orthology/paralogy relations.
This is a reimplementation of [ParaPhylo](http://pacosy.informatik.uni-leipzig.de/208-0-ParaPhylo.html) which uses heuristics for the NP-hard optimization steps instead of exact ILP solutions.

## Citation and References

If you use AsymmeTree in your project or code from it, please consider citing:

* **Stadler, P. F., Geiß, M., Schaller, D., López Sánchez, A., González Laffitte, M., Valdivia, D., Hellmuth, M., and Hernández Rosales, M. (2020) From pairs of most similar sequences to phylogenetic best matches. Algorithms for Molecular Biology. doi: 10.1186/s13015-020-00165-2.**

For references to concepts and algorithms that were implemented please see the [manual](https://github.com/david-schaller/AsymmeTree/blob/master/manual/AsymmeTreeManual.pdf).

Please report any bugs and questions in the [Issues](https://github.com/david-schaller/AsymmeTree/issues) section.
Also, feel free to make suggestions for improvement and/or new functionalities.