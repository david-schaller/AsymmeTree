

Each first level title ('#') is an issue.

Important to check the context:  [General context.md](General context.md) 

I paste the error at each title. AI agent can add notes in-line. Also add subsection 'AI notes'.

AI can write to suggest changes and can code to implement them.

Human curator has to approve to add '[solved]' to error title.

# 'DLListIterator' object is not iterable

I have two scripts:

-  [example_simulations.py](Example code/example_simulations.py) 
-  [simulate_HS_level.py](Example code/simulate_HS_level.py) 

- [ ] Which one is the latest?

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



