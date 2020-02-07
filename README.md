# Gene Tree Simulator and Best Match Inference

## Overview

### Simulator for species and gene trees

Framework for the realistic simulation of species and gene tree scenarios with asymmetric evolution rates. The following steps are implemented in the Python package `asymmetree.simulator`:
* species tree simulation ('innovation model')
* gene tree simulation (Gillespie algorithm)
* gene tree imbalancing (asymmetric evolution rates of paralogous genes)
* computation of a (noisy) distance matrix from the gene tree

### Best Match Inference

Inference of the best match relation either directly from the gene tree or from a distance matrix (several methods).
* `asymmetree.best_matches`

## Dependencies

AssymmeTree has several dependencies:
* [NetworkX](https://networkx.github.io/)
* [Scipy and Numpy](http://www.scipy.org/install.html)
* [Matplotlib](https://matplotlib.org/)

Furthermore, to use functions involving sequence simulation and alignment, the following packages must be installed:
* [Biopython](http://biopython.org/wiki/Download)
* [Pyvolve](https://github.com/sjspielman/pyvolve)

To use the tree reconstruction method for best match inference and the C++ implementation of the quartet method, resp., the following software must be installed
(I recommend that you compile these tools on your machine, place the binaries into a persistent location and add this location to your PATH environment variable):
* [RapidNJ](https://birc.au.dk/software/rapidnj/)
* [qinfer](https://github.com/david-schaller/qinfer)
