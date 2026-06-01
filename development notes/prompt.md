0. ~~Define new set of cases for interactors~~
   0.1 Add functions in code: `get_interactors`

1. Create the inverse maps mu_inv and mup_inv
2. Create function to obtain interactors via inverse map
3. ~~Generate an interactive dependency graph~~
   1. `dated_gene_tree_auxiliary` vs `create_auxiliary_tree`

---

Take context from `development notes/General context.md`, no other action required. Then go to `development notes/implementation.md`  and implement the task at lines 


---

- [ ] ~~Update creation of auxiliary tree (at the script `development notes/Example code/example_simulations.py`, and possibly also at dependences), following the updated algorithm at lines 137-144.~~
- [ ] ~~We need a function that given a gene branch $g$ in $L^0$, returns a dictionary of interactors, where keys and values are defined as in the table at lines 154-161.~~
  - [ ] ~~To fullfill this task, we need to first create an auxiliary inverse map for $\mu$  and $\mu'$.~~
  - [ ] ~~Follow algorithm at line XXX~~
- [ ] Ignore the 'level' attribute of branches, get rid of it.
