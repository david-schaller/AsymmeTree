# Gene Tree Simulator and Best Match Inference

### Simulator for species and gene trees

Framework for the realistic simulation of species and gene tree scenarios with asymmetric evolution rates. The following steps are implemented in the Python package `asymmetree.simulator`:
* species tree simulation ('innovations model')
* gene tree simulation (Gillespie algorithm)
* gene tree imbalancing (asymmetric evolution rates of paralogous genes)
* computation of a (noisy) distance matrix from the gene tree

### Best Match Inference

Inference of the best match relation either directly from the gene tree or from a distance matrix (several methods).
* `asymmetree.best_match`_inference (Python implementation)
* `qinfer` (C++ implementation, inference from distance matrix only)
