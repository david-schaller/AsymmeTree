# Implementation

The code should follow all AsymmeTree standards, in particular, have a look to  [skill contributing.md](skill contributing.md).

It is important to be aware of the  [general context.md](General context.md).

AI agent can add consice in-line comments and reply to questions in this file.

AI agent should generate a file `implemented.md` as a reply to the following list of tasks.

## Review theory and code

- [x] Have a look to the file  [Math and Algorithms.md](Math and Algorithms.md).
  - [x] Is this well detailed and consistent?
  - [x] Is there room for theoretical imprvement?
- [x] Have a look to the code.
  - [x] How complete is with respect to the theory?
  - [x] What are the missing steps?
  - [x] Is there something difficult to translate from theory to code?

### AI implementation notes

Reviewed in [implemented.md](implemented.md). Theory is detailed enough, but edge notation, the
time-restricted inverse map, and interaction-rate weights need tightening before implementation.
The code now has the main scaffolding, corrected auxiliary tree, `gamma_prime`, and six-class
interactor collection; C0 transfer weighting and weighted loss updates are still missing.

## Correct auxiliary tree in code

Te current implementation of the auxiliary tree is incorrect

- [x] Spot the mismatch and write it down.
- [x] Write the plan to upgrade the current implementation.
- [x] Correct the auxiliary tree so it follows lines 101-110 of `Math and Algorithms.md`.

### AI implementation notes

Updated `create_auxiliary_tree()` so `T_A` has shape `0_A -> rho_A -> {rho_H, 0_S}`. The host
component now starts at the canonical host root, `rho_A` and `0_A` use separate positive offsets
above `0_S`, and the symbiont planted root is marked as the transfer event.

## Collect interactors during gene tree simulation

- [x] Plan the changes in the code to:
  - [x] Construct the inverse maps
  - [x] Collect interactors
- [x] Construct `gamma_prime`, the inverse host-to-symbiont edge map, from the unpruned symbiont tree or from the annotated auxiliary tree.
- [x] Replace the current boolean `_branches_interact()` helper with an interactor collector returning
  the six classes `A0`, `A1`, `A2`, `B0`, `B1`, and `B2`.

### AI implementation notes

`GeneTreeSimulator` now builds `gamma_prime` from the annotated auxiliary tree and keeps
`species2genes` as the active `gamma` map. `_branches_interact()` was replaced by
`_collect_interactors()`, which returns `A0`, `A1`, `A2`, `B0`, `B1`, and `B2`; the old partner API
now unions those classes for existing rate hooks.

## Update rates based on interactions

- [x] Spot the places in the code where rates are updated or should be updated
- [x] plan the changes.

### AI implementation notes

Rate hooks are in `_new_branch()`, `_get_branch_and_type()`, `_run()`, `_speciation()`, `_loss()`,
and `_sample_recipient()`. The safer plan is to recompute loss-rate multipliers from current
interactor classes after branch-changing events, and apply the C0 penalty during recipient sampling.
