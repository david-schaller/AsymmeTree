# AsymmeTree

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![pypi version](https://img.shields.io/pypi/v/asymmetree)](https://pypi.org/project/asymmetree/)
[![Python](https://img.shields.io/pypi/pyversions/asymmetree)](https://pypi.org/project/asymmetree/)
[![PyPI downloads](https://img.shields.io/pypi/dm/asymmetree)](https://pypi.org/project/asymmetree/)
[![DOI](https://img.shields.io/badge/DOI-10.3390%2Fsoftware1030013-blue)](https://doi.org/10.3390/software1030013)

AsymmeTree is an open-source Python library for the simulation and analysis of phylogenetic
scenarios.
It includes a simulator for species and gene trees with heterogeneous evolution rates, nucleotide
and amino acid sequences with or without indels, as well as whole genomes/proteomes.

Moreover, it includes a matplotlib-based visualization of the simulated trees as well as tools for
the extraction of information from the simulated scenarios such as orthology, best matches, and
xenology.

The library is primarily designed to explore and validate mathematical concepts, and to test
inference methods for various steps on the way to more realistically-available data, i.e., dated
gene trees, additive distances of gene sets, noisy distances and finally sequences.

## Installation and usage

The package requires Python 3.10 or higher.
The `asymmetree` package is available on [PyPI](https://pypi.org/project/asymmetree/):

```bash
pip install asymmetree
```

See also the page [Installation](https://david-schaller.github.io/AsymmeTree/installation/) in the
documentation.

The [AsymmeTree documentation](https://david-schaller.github.io/AsymmeTree/) contains a user manual
with example code as well as the
[API reference](https://david-schaller.github.io/AsymmeTree/api/treeevolve/).

## Citation

If you use AsymmeTree in your project or code from it, please consider citing:

> David Schaller, Marc Hellmuth, and Peter F. Stadler.
> AsymmeTree: A Flexible Python Package for the Simulation of Complex Gene Family Histories.
> Software 2022, 1(3), 276-298; doi: 10.3390/software1030013

Please report any bugs and questions in the
[Issues](https://github.com/david-schaller/AsymmeTree/issues) section.
Also, feel free to make suggestions for improvement and/or new functionalities.
