# Changelog

All notable changes to this project will be documented in this file.

## [2.3.0] - unreleased

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
