

Each first level title ('#') is an issue.

Important to check the context:  [General context.md](General context.md) 

I paste the error at each title. AI agent can add notes in-line. Also add subsection 'AI notes'.

AI can write to suggest changes and can code to implement them.

Human curator has to approve to add '[solved]' to error title.

# 'DLListIterator' object is not iterable [solved]

I have two scripts:

-  [example_simulations.py](Example code/example_simulations.py) 
-  [simulate_HS_level.py](Archive/simulate_HS_level_legacy.py) 

- [x] Which one is the latest?
- [x] Explain the cause of the error and suggest how to fix it.

## AI notes

- [x] Which one is the latest?
`example_simulations.py` is the newer script in this repository. Its latest commit is `1db244d` on `2026-04-27` (`format corrections`), while `simulate_HS_level.py` was last changed in `f7a3342` on `2026-04-20` (`add mu map for host tree`). The larger diff between both files also suggests that `example_simulations.py` is the successor or refactor.

- [x] Explain the cause of the error and suggest how to fix it.
The error is caused by an obsolete `tralda` version in the failing environment. In `tralda==1.1.0`, `node.children` resolves to a `DLListIterator` that does not fully satisfy the Python iterable protocol under Python 3.13, so expressions such as `for child in node.children` and `",".join(... for child in node.children)` raise `TypeError: 'DLListIterator' object is not iterable`.

Suggested fix:

1. Upgrade `tralda` to a current release where the iterator bug is already fixed.
2. Keep the `asymmetree` code unchanged, because the original loops over `node.children` work again with the updated dependency.

Implementation note:
The temporary local shim can be removed after upgrading `tralda`. In the reported environment, upgrading from `tralda==1.1.0` to `tralda==2.0.3` is the correct fix.

Verification note:
I could not run the simulations in the current workspace because `tralda` is not installed here, but the upstream history confirms that this iterator bug was fixed in `tralda`, and the user reported upgrading successfully.

Both scripts raise the same error:

## Example errors

### `simulate_HS_level.py`

```bash
$ python simulate_HS_level.py 
Running simulations:   0%|                                                                           | 0/270 [00:00<?, ?it/s]
Traceback (most recent call last):
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/simulate_HS_level.py", line 108, in <module>
    Tg_full, Tg_prunned= simulate_HG_scenario(D_h_trees[hostTree_id], D, T, L, R_probability, T_bias)
                         ~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/simulate_HS_level.py", line 64, in simulate_HG_scenario
    gene_tree = dated_gene_tree(speciens_tree,
                                dupl_rate= D,
    ...<4 lines>...
                                transfer_distance_bias=T_bias,
                                )
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/asymmetree/treeevolve/genes.py", line 658, in dated_gene_tree
    gene_tree_simulator = GeneTreeSimulator(S)
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/asymmetree/treeevolve/genes.py", line 75, in __init__
    self._analyze_species_tree()
    ~~~~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/asymmetree/treeevolve/genes.py", line 215, in _analyze_species_tree
    self.S_subtree_survivors[u] = sum(self.S_subtree_survivors[v] for v in u.children)
                                  ~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/asymmetree/treeevolve/genes.py", line 215, in <genexpr>
    self.S_subtree_survivors[u] = sum(self.S_subtree_survivors[v] for v in u.children)
                                                                           ^^^^^^^^^^
TypeError: 'DLListIterator' object is not iterable

```



### `example_simulations.py`

```bash
$ python example_simulations.py 
Traceback (most recent call last):
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/example_simulations.py", line 330, in <module>
    host_series, symbiont_dataframe = run_example_simulations(output_dir= 'simulations/')
                                      ~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/example_simulations.py", line 70, in run_example_simulations
    host_records[host_tree_id] = to_nhx(host_tree)
                                 ~~~~~~^^^^^^^^^^^
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/asymmetree/hologenome/hologenome_simulation.py", line 201, in to_nhx
    return _serialize(tree.root) + ";"
           ~~~~~~~~~~^^^^^^^^^^^
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/asymmetree/hologenome/hologenome_simulation.py", line 188, in _serialize
    children = ",".join(_serialize(child) for child in node.children)
  File "/scr/k80san/jantonio/AsymmeTree/development notes/Example code/asymmetree/hologenome/hologenome_simulation.py", line 188, in <genexpr>
    children = ",".join(_serialize(child) for child in node.children)
                                                       ^^^^^^^^^^^^^
TypeError: 'DLListIterator' object is not iterable

```



