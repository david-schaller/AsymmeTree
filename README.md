# AsymmeTree

[![pypi version](https://img.shields.io/badge/pypi-v0.0.1-blue.svg)](https://pypi.org/project/asymmetree/0.0.1/)

AsymmeTree is Python3 library for the simulation and analysis of phylogenetic scenarios.
It includes a simulator for species and gene tree scenarios with asymmetric evolution rates and tools for the inference and analysis of best matches (resp. best hits) and orthology.

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

### Simulator for species and gene trees

The following steps are implemented in the Python package `asymmetree.simulator`:
* species tree simulation ('innovation model')
* gene tree simulation (Gillespie algorithm)
* gene tree imbalancing (asymmetric evolution rates of paralogous genes)
* computation of a (noisy) distance matrix from the gene tree

### Best Match Inference

Inference of the best match relation either directly from the gene tree or from a distance matrix (several methods).
* `asymmetree.best_matches`

# REFERENCES
* Stadler, P. F., Geiß, M., Schaller, D., López Sánchez, A., González Laffitte, M., Valdivia, D., Hellmuth, M., and Hernández Rosales, M. (2019) From Best Hits to Best Matches. Submitted to Algorithms for Molecular Biology.

* Geiß, M., Chávez, E., González Laffitte, M., López Sánchez, A., Stadler, B. M. R., Valdivia, D. I., Hellmuth, M., Hernández Rosales, M., and Stadler, P. F. (2019) Best match graphs. Journal of Mathematical Biology, 78(7):2015-2057. ISSN 0303-6812, 1432-1416. doi: 10.1007/s00285-019-01332-9.

* Hellmuth, M., Wieseke, N., Lechner, M., Lenhof, H.-P., Middendorf, M., and Stadler, P. F. (2015) Phylogenomics with paralogs. PNAS, 112(7):2058-2063. doi: 10.1073/pnas.1412770112.

* Lechner, M., Findeiß, S., Steiner, L., Marz, M., Stadler, P. F., and Prohaska, S. J. (2011). Proteinortho: detection of (co-) orthologs in large-scale analysis. BMC bioinformatics, 12(1), 124. ISSN 1471-2105. doi: 10.1186/1471-2105-12-124.
