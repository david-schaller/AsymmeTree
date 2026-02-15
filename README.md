# AsymmeTree
![Logo](resources/images/logo.png)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![pypi version](https://img.shields.io/badge/pypi-v2.2.2-blue.svg)](https://pypi.org/project/asymmetree/)

AsymmeTree is an open-source Python library for the simulation and analysis of phylogenetic scenarios.
It includes a simulator for species and gene trees with heterogeneous evolution rates, nucleotide and amino acid sequences with or without indels, as well as whole genomes/proteomes.

Moreover, it includes a matplotlib-based visualization of the simulated trees as well as tools for the extraction of information from the simulated scenarios such as orthology, best matches, and xenology.

The library is primarily designed to explore and validate mathematical concepts, and to test inference methods for various steps on the way to more realistically-available data, i.e., dated gene trees, additive distances of gene sets, noisy distances and finally sequences.

## Installation

The package requires Python 3.10 or higher.
The `asymmetree` package is available on [PyPI](https://pypi.org/project/asymmetree/):

```bash
pip install asymmetree
```

Alternatively, you can clone the repo, go to the root folder of package and install it using the
command:

```bash
pip install .
```

AsymmeTree has several dependencies (which are installed automatically when using `pip`):

* [NetworkX](https://networkx.github.io/)
* [Scipy and Numpy](http://www.scipy.org/install.html)
* [Matplotlib](https://matplotlib.org/)
* [tralda](https://github.com/david-schaller/tralda)

## Usage and description

A user manual with example code can be found in the [Wiki](https://github.com/david-schaller/AsymmeTree/wiki/Manual).
AsymmeTree is divided into several subpackages and modules, see also the following [documentation](https://david-schaller.github.io/docs/asymmetree/) generated from the source code.

## Contributing

If you want to contribute to `asymmetree`, please use the package and project manager
[uv](https://docs.astral.sh/uv/).
See [this page](https://docs.astral.sh/uv/getting-started/installation/) for installation
instructions.

Moreover, please use [pre-commit](https://pre-commit.com) for automated code formatting and linting.

To install it and initialize it for your local `asymmetree` repository, follow these steps:

- Install `uv`
- Run the following command (after which you should be able to run `pre-commit` from anywhere)
    - `uv tool install pre-commit --with pre-commit-uv`
- Navigate to the root directory of your local `asymmetree` repository
    - `cd <MY_PATH_TO>/asymmetree`
- Install `pre-commit` as a git hook in the `asymmetree` repository
    - `pre-commit install`

## Citation

If you use AsymmeTree in your project or code from it, please consider citing:

> David Schaller, Marc Hellmuth, and Peter F. Stadler.
> AsymmeTree: A Flexible Python Package for the Simulation of Complex Gene Family Histories.
> Software 2022, 1(3), 276-298; doi: 10.3390/software1030013

Please report any bugs and questions in the
[Issues](https://github.com/david-schaller/AsymmeTree/issues) section.
Also, feel free to make suggestions for improvement and/or new functionalities.
