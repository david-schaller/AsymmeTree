

---

- [x] I Have to adjust the parameters of $\alpha$ and $\beta$ in order to generate evolutionary histories with only 'host keps' and only 'symbiont keeps' simulations.
- [ ] ~~Change parameters for the number of simulations:~~
  - how many species trees? hoy many replications per species tree? how many simbionts? 
  - define these numebrs based on Sanchita data
- [x] We have to add a human-readable docstring to 'example simulations' and inliene comments
- [ ] We have to create a documentation of the generated dataset
- [ ] Should we refactor `example_simulations.py`?

---

- [x] Add a new task to correct `example_simulations.py`

  Separate rates bellow

  > \- ``symbiont_dtl_rates`` provides candidate ``(duplication, transfer, loss)``
  >
  >   triples. Each triple is used for the symbiont simulation and then reused for
  >
  >   the gene-tree simulation inside the resulting auxiliary tree.

- [x] `example_simulations.py` standalone

  - [x] to provide input parameters, let's provide an input file with parameters.

- [ ] ~~Also output a file describing the output (number of trees)~~

- [x] If output directory is non-specified, generate a new one with date and unique number.

- [x] How is handled the map $\mu$ when the auxiliary tree is separated into host and symbiont again?

  - [x] Add to the teoretical markdown a note about

    > the host-side and symbiont-side projections of those gene trees.

    

Describe in the theoretical markdown how projection of maps are handled in the code. Since you write in the documentation 'the host-side and symbiont-side projections of those gene trees.'

---

1. `dated_gene_tree_auxiliary` vs `create_auxiliary_tree`

---

Take context from `development notes/General context.md`, no other action required. Then go to `development notes/implementation.md`  and implement the task at lines 

