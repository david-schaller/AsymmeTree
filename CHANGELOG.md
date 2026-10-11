# Changelog

All notable changes to this project will be documented in this file.

## [2.4.0] - 2026-10-11

### 🚨 Breaking changes

- The tree visualization now uses the plotting functionality of
  [tralda](https://github.com/david-schaller/tralda) instead of an own matplotlib implementation.
  - The deprecated class `Visualizer` was removed; use `visualize()` instead.
  - The argument `save_as` of `visualize()` is deprecated in favor of `path`. The figure is no
    longer displayed by default; use `show=True` to display it.
  - `visualize()` now returns the figure (`(fig, ax)` for matplotlib, a plotly `Figure` for plotly).
- matplotlib is no longer a required dependency. Install it with
  `pip install asymmetree[matplotlib]` (or `asymmetree[viz]` for matplotlib and plotly).

### 🌟 Features

- `visualize()` now supports the plotly backend (`backend="plotly"`) for interactive figures.
- `visualize()` now supports vertical and circular layouts (`layout_mode`) and the edge length modes
  `"attr"`, `"uniform"`, `"even"`, and `"rank"` (`edge_length_mode`).
- Custom styles for nodes and edges can be passed to `visualize()` as a tralda `TreeStyle`.
- `visualize()` supports the arguments `species_info` (append the species to the leaf labels) and
  `distance_axis` (show the distance from the root; matplotlib backend with horizontal layout and
  `"attr"` edge lengths only).
- `visualize()`, `assign_colors()`, `assign_species_colors()`, and `assign_gene_colors()` can now be
  imported from `asymmetree.visualization`.

### ♻️ Refactorings

- Color assignment for the visualization no longer imports `matplotlib.pyplot`, and matplotlib is
  only imported when needed.
- Removed the `develop` branch; all development and pull requests now target `main`.

### 📚 Documentation

- Rewrote the visualization guide: it now describes the AsymmeTree default style, layouts, edge
  length modes, the distance axis, renderer options, custom styles, and the plotly backend, with
  example figures.
- Updated the example script `examples/tree_visualization.py` such that it generates the figures of
  the visualization guide with a fixed random seed.
- The contributing guide now refers to `main` instead of `develop` for branching and pull requests.
- Added a banner to the development documentation stating that it is built from `main` and may
  contain unreleased changes.

### 📦 Build

- Updated tralda to version 2.1.0.
- Added the optional dependencies `matplotlib`, `plotly`, and `viz`.
- Bumped indirect dependency `urllib3` to 2.8.0.

### 🔖 Release

- The stable documentation is now deployed only when a release is published, so it always matches
  the latest version on PyPI.
- The development documentation (`/dev/`) is now deployed on every push to `main`.
- Pull requests now run a strict documentation build (without deployment) to catch broken
  documentation early.

## [2.3.2] - 2026-09-12

### 📦 Build

- Bumped indirect dependencies: `pymdown-extensions` to 11.0.2, `pillow` to 12.3.0.

## [2.3.1] - 2026-05-28

### 📚 Documentation

- Added a section on contributing to the documentation, including guidelines for code style,
  testing, and pull requests.

### 📦 Build

- Added the `classifiers` field to `pyproject.toml` for better visibility on PyPI.
- Bumped indirect dependencies: `urllib3` to 2.7.0, `pymdown-extensions` to 10.21.3, `idna` to 3.15.

## [2.3.0] - 2026-04-16

### 🚨 Breaking changes

- Renamed all modules using snake_case instead of camelCase.
  - For example, `asymmetree.treeevolve.SpeciesTree` was renamed to `asymmetree.treeevolve.species`.
  - However, most of the API functions are not affected as long as re-imports on subpackage
    level are used (e.g., `from asymmetree.treeevolve import species_tree_n`).

### ♻️ Refactorings

- The functions in module `asymmetree.treeevolve.species` now use explicit arguments instead of
  `**kwargs` for better readability and usability.

### 🐛 Bug fixes

- Edges for edge contraction to generate multifurcations are now correctly sampled without
  replacement.
- Fixed undefined variable `mean` and bounds checks falsely treating 0/0.0 as unset in
  `sampling.py`.

### 📚 Documentation

- Moved the documentation from the Wiki to a folder `docs` in the repository.
- Updated the documentation to reflect the changes in the API and the new features.
- Set up [MkDocs](https://www.mkdocs.org/) for generating the documentation from the `docs`
  folder.

### 🎨 Style

- Changed to Google style for docstrings.
- Introduced typing hints for all functions and methods.
- Maximal line length is now 100 characters.

### 📦 Build

- Introduced [uv](https://docs.astral.sh/uv/) as package and project manager.
- Updated tralda to version 2.0.0.

### 🔖 Release

- Added a release workflow using GitHub Actions that automatically builds and publishes the
  package to PyPI on new releases.

## [2.2.2] - 2025-10-26

### 🐛 Bug fixes

- Updated tralda to version 1.1.1, which fixes a bug occurring with Python 3.13.

## [2.2.1] - 2023-04-26

### 🐛 Bug fixes

- Fixed a bug in the module for rate heterogeneity: The `base_rate_sampler` falsely used
  `dupl_rate` instead of `base_rate`.

## [2.2.0] - 2022-05-30

### 🌟 Features

- The sampling of species trees can now be conditioned on the number of leaves and the age of the
  tree.
- Gene conversion is now supported as an event shaping the gene tree (in addition to duplication,
  loss, and HGT).
- Additive and replacing HGT with transfer distance bias.
- Matplotlib-based visualization of trees now also accepts species trees as input.

### ♻️ Refactorings

- The innovation model is now an option that modifies the sampling of species lineages for the
  next bifurcation rather than its own model.
- A number of functions have been re-named to be shorter and/or more meaningful, e.g.,
  `observable_tree()` is now `prune_losses()`.
- The reconciliation attribute of the gene tree nodes has been re-named from `color` to `reconc`.

### ⚡️ Performance

- The Gillespie mode (formerly jump chain) for sequence simulation now runs faster.

### 🐛 Bug fixes

- The distance attribute of the trees is now updated after contraction of edges to generate
  multifurcations.
- Evolution rate heterogeneity now correctly handles multifurcating duplication nodes.

### 📚 Documentation

- Added docstrings to all functions and methods.

## [2.1.0] - 2021-09-28

### ♻️ Refactorings

- Major refactoring of the subpackage and module structure.

## [2.0.1] - 2021-07-01

### ♻️ Refactorings

- The tree data structures are now in the package
  [tralda](https://github.com/david-schaller/tralda).

## [1.0.2] - 2021-01-12

## [1.0.1] - 2021-01-11

## [1.0.0] - 2020-12-17

## [0.1.0] - 2020-06-10

### 🌟 Features

- Added the simulation of nucleotide and amino acid sequences
- Added the generation of genomes/proteomes.

## [0.0.5] - 2020-03-12

## [0.0.1] - 2020-02-18
