![AsymmeTree Logo](manual/images/logo.png)
# AsymmeTree

[![license](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![pypi version](https://img.shields.io/badge/pypi-v0.0.1-blue.svg)](https://pypi.org/project/asymmetree/0.0.1/)

*Copyright © 2020 David Schaller*

AsymmeTree is Python library for the simulation and analysis of phylogenetic scenarios.
It includes a simulator for species and gene tree scenarios with asymmetric evolution rates, tools for the inference and analysis of best matches (resp. best hits) and orthology, as well as an algorithm to compute supertrees.

## Installation

#### Easy installation with pip

The `asymmetree` package is available on PyPI:

    pip install asymmetree

For details about how to install Python packages see [here](https://packaging.python.org/tutorials/installing-packages/).
    
#### Dependencies

AssymmeTree has several dependencies (which are installed automatically when using `pip`):
* [NetworkX](https://networkx.github.io/)
* [Scipy and Numpy](http://www.scipy.org/install.html)
* [Matplotlib](https://matplotlib.org/)

Furthermore, to use functions involving sequence simulation and alignment, the following packages must be installed (i.e., they are not installed automatically!):
* [Biopython](http://biopython.org/wiki/Download)
* [Pyvolve](https://github.com/sjspielman/pyvolve)

To use the tree reconstruction method for best match inference and the C++ implementation of the quartet method, resp., the following software must be installed
(I recommend that you compile these tools on your machine, place the binaries into a persistent location and add this location to your PATH environment variable):
* [RapidNJ](https://birc.au.dk/software/rapidnj/)
* [qinfer](https://github.com/david-schaller/qinfer)

## Overview

### Tree data structures

The two classes `Tree` (in `asymmetree.tools.Tree`) and `PhyloTree` (in `asymmetree.tools.PhyloTree`, inherits from `Tree`) implement tree data structures which are essential for most of the modules in the package.
The latter contains converters and parsers for the Newick format and a NetworkX graph format.

### Simulator for species and gene trees

The following steps are implemented in the Python package `asymmetree.simulator`:
* species tree simulation ('innovation model')
* gene tree simulation (Gillespie algorithm)
* gene tree imbalancing (asymmetric evolution rates of paralogous genes)
* computation of a (noisy) distance matrix from the gene tree

### Best Match Inference

Inference of the best match relation either directly from the gene tree or from a distance matrix (several methods).
* `asymmetree.best_matches`

### Supertree Computation

Implementation of the BuildST algorithm described by Deng & Fernández-Baca (2016) to compute a supertree from a given list of tree based on the leaf labels. The algorithm uses the dynamic graph data structure described by Holm, de Lichtenberg and Thorup in 2001 (HDT algorithm).
* `asymmetree.tools.BuildST`
* `asymmetree.tools.hdtgraph.DynamicGraph`

### Cograph editing and ParaPhylo

The subpackages `asymmetree.cograph` and `asymmetree.proteinortho` contain heuristics for cograph editing and a method to compute rooted species tree from orthology/paralogy relations.
The latter is a reimplementation of [ParaPhylo](http://pacosy.informatik.uni-leipzig.de/208-0-ParaPhylo.html) which uses heuristics for the NP-hard steps instead of exact ILP solutions.


## REFERENCES

* <small>**Stadler, P. F., Geiß, M., Schaller, D., López Sánchez, A., González Laffitte, M., Valdivia, D., Hellmuth, M., and Hernández Rosales, M. (2019) From Best Hits to Best Matches. Submitted to Algorithms for Molecular Biology.**</small>

* <small>Deng, Y. and Fernández-Baca, D. (2016) Fast Compatibility Testing for Rooted Phylogenetic Trees. 27th Annual Symposium on Combinatorial Pattern Matching (CPM 2016). doi: 10.4230/LIPIcs.CPM.2016.12.</small>

* <small>Geiß, M., Chávez, E., González Laffitte, M., López Sánchez, A., Stadler, B. M. R., Valdivia, D. I., Hellmuth, M., Hernández Rosales, M., and Stadler, P. F. (2019) Best match graphs. Journal of Mathematical Biology, 78(7):2015-2057. ISSN 0303-6812, 1432-1416. doi: 10.1007/s00285-019-01332-9.</small>

* <small>Hellmuth, M., Wieseke, N., Lechner, M., Lenhof, H.-P., Middendorf, M., and Stadler, P. F. (2015) Phylogenomics with paralogs. PNAS, 112(7):2058-2063. doi: 10.1073/pnas.1412770112.</small>

* <small>Holm, J., de Lichtenberg, K., and Thorup, M. (2001) Poly-logarithmic deterministic fully-dynamic algorithms for connectivity, minimum spanning tree, 2-edge, and biconnectivity. J. ACM, 48(4):723–760. doi: 10.1145/502090.502095.</small>

* <small>Lechner, M., Findeiß, S., Steiner, L., Marz, M., Stadler, P. F., and Prohaska, S. J. (2011). Proteinortho: detection of (co-) orthologs in large-scale analysis. BMC bioinformatics, 12(1), 124. ISSN 1471-2105. doi: 10.1186/1471-2105-12-124.</small>

* <small>Rauch Henzinger, M. and King, V. (1999) Randomized fully dynamic graph algorithms with polylogarithmic time per operation. J. ACM 46(4):502–536. doi: 10.1145/225058.225269.</small>