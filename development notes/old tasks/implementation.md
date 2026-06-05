# Implementation

The code should follow all AsymmeTree standards, in particular, have a look to  [skill contributing.md](skill contributing.md).

## Initialize `hologenome` module

- [x] Create a module `asymmetree.hologenome` inside of the directory `development notes/Example code/`
- [x] It should have the same structure as `asymmetree.genome`
- [x] It should be initialized with the code in `development notes/Archive/simulate_HS_level_legacy.py`. Until now, it simulates only the 'holobiont'.
- [x] Create the auxiliary tree $T_A$​ (One per each host-symbiont pair of trees).
- [x] Create an example script `development notes/Example code/example_simulations`. The result should be:
  - [x] One pandas series with Host trees trees (in general nhx format)
  - [x] One pandas dataframe with with symbiont trees  (in general nhx format)
- [x] Create a new copy of the whole package inside `development notes/Example code/` and add the new code.

### AI implementation notes

- The prototype now lives in `development notes/Example code/hologenome` with the same package
  pattern as `asymmetree.genome`: a thin `__init__.py` and a single simulation module.
- `HologenomeSimulator` currently covers the holobiont layer only and keeps the prototype logic
  close to `simulate_HS_level_legacy.py`, but removes the dependency on `revolutionhtl` by exporting
  trees with a local NHX serializer.
- The auxiliary tree uses collision-safe `H*`, `S*`, and `A*` labels and stores the host-side map
  in `reconc`, which gives us a clean hand-off for the later gene-level extension.
- `T_A` is built from the full symbiont tree rather than the pruned one so that later
  gene-level work still has access to extinct symbiont branches.
- The example script returns a host-tree `Series` and a symbiont-scenario `DataFrame`; the
  dataframe also includes the serialized auxiliary tree for each simulated host-symbiont pair.



## Add `holoevolve` module

- [x] Create a copy of `asymmetree.treeevolve` inside of the directory `development notes/Example code/asymmetree/`. Name this new module `holoevolve`
- [x] Add code to the `hologenome` module: simulate a single gene tree inside of the auxiliary tree $T_A$​, use the functions from `holoevolve` instead of `treeevolve`.
- [x] Make sure this new code is effectively called at the example script `development notes/Example code/example_simulations`. The result should be again a pandas dataframe.
- [ ] We need to obtain the *inverse* maps of $\mu$ and $\mu'$
  - This should be done when simulation the trees, i.e. when simulation symbiont tree we construct $\mu'^{(-1)} : T_S \rightarrow T_H$, and when simulating gene trees we construct $\mu^{(-1)} : T_G \rightarrow T_S$.
  - These *inverse maps* 
- [ ] We need a function that given a gene branch $g$ in $L^0$ during gene tree simulation​, returns a dictionary of interactions, where keys and values are defined as in the table at  [Math and Algorithms.md](Math and Algorithms.md)  lines 116-123.
  - [ ] To fulfill this task, we need to first create an auxiliary inverse map for $\mu$  and $\mu'$​.
    - Note that the maps $\mu$ and $\mu'$ are encoded as attributes of nodes in trees.
    - Make sure this is not already implemented.
  - [ ] Follow algorithm at line XXX...


### AI implementation notes

- The holobiont stage still uses `treeevolve`: host trees come from `treeevolve.species_tree_n()`
  and `HologenomeSimulator.simulate_symbiont_trees()` keeps using
  `treeevolve.dated_gene_tree()` / `treeevolve.prune_losses()`.
- `HologenomeSimulator.simulate_gene_trees()` simulates one gene tree per stored auxiliary tree and
  keeps both the unpruned and pruned versions alongside the host-symbiont results, using
  `holoevolve.dated_gene_tree()` / `holoevolve.prune_losses()`.
- `example_simulations.py` now runs the auxiliary-tree gene simulation for every scenario and adds
  `T_gene_unpruned` and `T_gene_pruned` to the returned scenario dataframe.
- The example currently reuses each scenario's DTL, replacement, and transfer-bias settings for
  the auxiliary-tree gene simulation so that every row remains self-contained.

## Improve auxiliary tree

- [ ] Update creation of auxiliary tree (at the script `development notes/Example code/example_simulations.py`, and possibly also at dependences), following the updated algorithm at  [Math and Algorithms.md](Math and Algorithms.md), lines 97-106.

### AI implementation notes



## Add homolog-homolog interactions

Right now asymmetree assumes that all the branches of the gene tree have exactly the same set of rates, which are the user-provided. We will change this.

- [x] Edit the module `holoevolve` in such a way that now each gene-branch in `self.branches` has their own set of DTLG rates. At the beginning of the simulation the rates are those provided by the user.
- [x] Update `_get_branch_and_type` to choose branch and event considering that not all of the branches have the same rate.
- [x] If there is a transfer event to a auxiliary-branch (i.e. a branch in the auxiliary tree), then all those gene-branches (directly or indirectly) residing on that branch should have an updated loss rates defined by rules (A0) to (A2), the new rate is $l'=\delta_l l$ for $\delta_l$​ being an user-provided value, by default it's value is two. If there are consecutive transfers to the same branch, the rate can increase multiple times
  - [ ] Correction: the loss rate should increase after a symbiont-transfer event instead of gene transfer.
    - [x] Comment the increase of rate after gene transfer
    - [ ] In the `_run` function, at line 490, the rates should be updated for each branch by checking con conditions (A0) to (A2), or after detecting that a gene were inherited to a new branch of the auxiliary tree who was caused by transfer event, this have to be checked at line 469.

- [x] If there is a gene loss for a branch with rate bigger than the user provided, other gene-branches following conditions (A0) to (A2) should decrease the rate as $l'=l/\delta_l$​. If there are consecutive losses the rate should be updated each time. The minimum value for loss rate is the user provided.

### AI implementation notes

- Active branches in `holoevolve` now carry their own DTLG rates, and descendant branches inherit the
  current rates of the lineage they split from.
- `_get_branch_and_type()` now samples branches proportionally to their current total rate and then
  samples the event type from that branch-specific rate vector.
- Homolog interactions are enabled only on auxiliary-tree branches annotated as `host` or
  `symbiont`; ordinary `holoevolve` runs on unannotated trees keep the old behavior.
- The previous gene-HGT rate-increase hook is now commented out because loss-rate increases should
  be triggered by symbiont-level transfer inheritance in `_run()` instead.
- A loss still divides surviving partners by `delta_l` without going below the user-provided base
  loss rate.
- The current A2 handling is symmetric: host and symbiont branches are treated equally for the
  multiplier update. The optional `S-keeps` / `H-keeps` asymmetry is still open for a later step.

## Prioritize intra-transfers

- [ ] 

## Separate auxiliary tree

- [x] After the simulation of genes inside of the holobiont, separate the auxiliary tree into host and symbiont again (maybe is better to just recall the original trees),
- [x] additionally separate the gene tree into two subtrees, one whose root is mapped to the symbiont side , and the other whose root is mapped to the host side of the auxiliary tree.
- [x] Keep the auxiliary tree and corresponding gene tree

### AI implementation notes

- `split_auxiliary_tree()` returns copied host and symbiont components from the already annotated
  auxiliary tree, so the original auxiliary tree remains available unchanged.
- `split_gene_tree_by_auxiliary_level()` projects a gene tree onto host-side and symbiont-side
  reconciliations using the auxiliary-tree `level` annotations.
- `HologenomeSimulator.simulate_gene_trees()` keeps the full auxiliary tree and full gene tree, and
  additionally stores host/symbiont auxiliary components plus unpruned and pruned side-specific
  gene trees.
- The example dataframe now includes the full `T_auxiliary` / `T_gene_*` columns as well as
  `T_auxiliary_host`, `T_auxiliary_symbiont`, and host/symbiont gene-subtree columns.

## Save trees in three different formats

- [x] In the example simulations script, save the trees in three different formats, like in the script `simulate_HS_level_legacy.py`. Use revolutionhtl if necesary.
- [x] Also save species trees as a csv table

### AI implementation notes

- Added a local `to_simple_newick()` serializer for compact `event|label->reconc` Newick output,
  with an optional dated variant that emits `dist` as branch lengths.
- `run_example_simulations(..., output_dir=...)` now writes NHX, simple, and simple-dated CSVs
  for host trees, scenario trees, and a separate species-tree table. The function return value
  remains the original NHX host series and scenario dataframe.
