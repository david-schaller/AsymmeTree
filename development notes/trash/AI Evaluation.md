The proposal is compatible with the current architecture of AsymmeTree, and the prototype in `development notes/Example code/simulate_HS_level.py` already shows that a first host-symbiont layer can be built by reusing the existing tree simulation code, the `reconc` attribute for reconciliation maps, and the dated-event information stored in `event` and `tstamp`. This makes the proposal more than a purely conceptual extension: the main simulation ingredients already exist.

At the same time, the prototype also shows the main challenge. A research script can assemble host and symbiont scenarios and post-process labels afterwards, but a library-quality implementation will need explicit collision-safe identifiers, helper functions for coexistence queries, and a clean API in a new `asymmetree.holobiont` module. The strongest part of the proposal is still the reuse of the current birth-death machinery, especially by modifying loss sampling and HGT recipient selection. The main technical risk is no longer basic feasibility, but precise model specification: cases (A0)-(A3) should be translated into explicit weighting rules and tested carefully, otherwise the behavior of the simulator will be difficult to validate.

# Plan for modification of the code



1. Use the prototype in `development notes/Example code/simulate_HS_level.py` as the starting point for the host-symbiont level, since it already reuses `species_tree_n`, `dated_gene_tree`, `prune_losses`, and the `reconc` metadata.
2. Move that logic into a new public module `asymmetree.holobiont`, following the same high-level structure as `asymmetree.genome`.
3. Replace script-level relabeling with package-level helper functions for collision-safe node identifiers and reconciliation formatting.
4. Add helper functions that detect coexisting host, symbiont, and gene lineages from the dated trees using `reconc`, `event`, and `tstamp`.
5. Add a gene-tree simulation step on the auxiliary tree $T_A$, reusing the current logic of `asymmetree.treeevolve.genes.GeneTreeSimulator` as much as possible.
6. Add configurable interaction parameters for cases (A0), (A1), (A2), and (A3), including the host-keeps or symbiont-keeps asymmetry.
7. Modify event sampling so loss rates can be increased for interacting coexisting genes, and HGT recipient choice can be decreased across different hosts.
8. Add documentation and tests for the new module, reconciliation maps, coexistence queries, collision-safe labels, and the new rate modifiers.

# AsymmeTree vs SaGePhy

SaGePhy is attractive because it already supports nested simulations beyond the usual species-gene setting. In particular, it can simulate gene trees inside species trees and domain or subgene trees inside one or more gene trees, including additive and replacing transfers and distance-biased transfers. This means that, if the goal were mainly to simulate another generic nested reconciliation process, SaGePhy already provides ideas and machinery that are close in spirit to the present project.

However, the match is not exact. SaGePhy is designed around species-gene-domain evolution, while the present proposal is about host-symbiont-gene evolution with host-dependent coexistence effects on losses and transfers. Those biological rules are not just another copy of domain evolution, so using SaGePhy directly would still require adapting its model assumptions. In addition, our current work is already inside the AsymmeTree repository, where tree simulation, reconciliation information, pruning, rate heterogeneity, and sequence simulation are available in one Python package.

Advantages of modifying AsymmeTree:

1. It preserves the current workflow and codebase, including `treeevolve`, `genome`, pruning, and sequence simulation.
2. The existing node attributes `reconc`, `event`, and `tstamp` already provide much of the information needed for a three-level extension.
3. The prototype in `development notes/Example code` shows that a first host-symbiont layer can already be assembled with the current tools.
4. It is easier to introduce a custom host-symbiont-gene interpretation than to retrofit that meaning into a domain-evolution framework.

Pitfalls of modifying AsymmeTree:

1. The current simulator is built around a single guiding tree for gene evolution, so a clean three-level implementation will require non-trivial refactoring.
2. Coexistence-aware rate updates are new behavior and will need careful specification and testing.
3. More development effort is needed before the new functionality reaches the same level of polish as the existing genome workflow.

Advantages of using SaGePhy:

1. It already supports nested phylogenetic simulation beyond two levels.
2. It already includes additive and replacing transfers, distance-biased transfer models, and birth away from the root.
3. Its structure may provide useful design inspiration for how to organize a multi-level simulator.

Pitfalls of using SaGePhy:

1. Its native model is species-gene-domain rather than host-symbiont-gene, so the biological interpretation is different from the target problem.
2. It would introduce an external toolchain based on Java plus supplementary Python scripts, instead of extending the current Python package directly.
3. It is less naturally integrated with AsymmeTree's current sequence simulation and package structure.
4. Moving to SaGePhy would likely reduce immediate reuse of the code already prototyped in this repository.

Overall, modifying AsymmeTree appears to be the better path if the main objective is to add a native host-symbiont-gene simulator to this package. SaGePhy is still valuable as a conceptual reference, especially for nested simulation design and transfer modeling, but it is not a drop-in replacement for the proposed feature.

# Rates

In SaGePhy, the duplication rate is the event-rate parameter for duplications during guest/gene-tree simulation inside the guide species tree. In the manual, GuestTreeGen takes dup rate, loss rate, and trans rate, and says that a higher duplication rate gives more frequent duplications. The SaGePhy paper also treats these as rates per unit branch length / time in the birth-death simulation sense.

So for your practical question: yes, it is essentially the same kind of quantity as dupl_rate in AsymmeTree.

In AsymmeTree, dupl_rate is also an event-rate parameter in the dated gene-tree birth-death process. It contributes to the total event rate together with loss, HGT, and gene conversion, and higher dupl_rate means duplications are sampled more often during simulation.

The important caveat is this:

- Conceptually equivalent: yes
- Guaranteed numerically identical across tools: no

They are not necessarily interchangeable one-to-one because the simulators differ in surrounding mechanics, for example:

- AsymmeTree also includes gc_rate for gene conversion.
- AsymmeTree handles some extra behaviors like prohibit_extinction and replace_prob.
- SaGePhy has its own options such as gene birth away from the root and its own transfer/replacement handling.

So if you see [0.133, 0.266, 0.266] in your script, the first number is meaningfully the same kind of “duplication intensity” as AsymmeTree’s dupl_rate, but matching outputs between the two tools would still require aligning the whole model, not just that one parameter.

Sources:

- SaGePhy manual PDF: https://compbio.engr.uconn.edu/wp-content/uploads/sites/2447/2018/10/sagephy-manual.pdf
- SaGePhy paper: https://academic.oup.com/bioinformatics/article/35/18/3496/5305633
- AsymmeTree gene-tree docs: docs/guide/gene_trees.md
- AsymmeTree implementation: asymmetree/treeevolve/genes.py

If you want, I can also map SaGePhy’s dup/loss/trans parameters to the closest AsymmeTree settings for your current script.
