# Implementation

The code should follow all AsymmeTree standards, in particular, have a look to  [skill contributing.md](skill contributing.md).

It is important to be aware of the  [general context.md](General context.md).

AI agent can add consice in-line comments and reply to questions in this file.

AI agent should generate a file `implemented.md` as a reply to the following list of tasks.
In this file, AI should only change checkbox states unless the user explicitly asks for another kind of edit.
All AI notes and explanations should be written in `implemented.md`.

## Review theory and code

- [x] Have a look to the file  [Math and Algorithms.md](Math and Algorithms.md).
  - [x] Is this well detailed and consistent?
  - [x] Is there room for theoretical imprvement?
- [x] Have a look to the code.
  - [x] How complete is with respect to the theory?
  - [x] What are the missing steps?
  - [x] Is there something difficult to translate from theory to code?

- [ ] Provide a deeper explanation of the following. Add it to 'implemented.md'.

    > A conservative first default is alpha = 0.25 * r_l and beta = 0.10 * r_l; a neutral baseline is alpha = beta = 0.15 * r_l.

    It looks good, nevertheless I need a deeper explanation of the meaning of multiplying alpha/beta with the 'relative gene content' $\frac{|\gamma(s)|}{M/N}$, and the implications for the probability of this event happening, i.e. the event-sampling process using a wheighted random choice. What probability distribution genes follow, given the time-sampling process? based on this, explain why to use alpha = 0.25 * r_l and beta = 0.10 * r_l.

    > A cap such as
    > `crowding = min(|gamma(.)| * N / M, c_max)` with `c_max` around `2` or `3` is worth keeping in
    > reserve.

    why to choose 2 and 3 from a theoretical point of view under the same arguments above?


## Correct auxiliary tree in code

Te current implementation of the auxiliary tree is incorrect

- [x] Spot the mismatch and write it down.
- [x] Write the plan to upgrade the current implementation.
- [x] Correct the auxiliary tree so it follows lines 101-110 of `Math and Algorithms.md`.


## Collect interactors during gene tree simulation

- [x] Plan the changes in the code to:
  - [x] Construct the inverse maps
  - [x] Collect interactors
- [x] Construct `gamma_prime`, the inverse host-to-symbiont edge map, from the unpruned symbiont tree or from the annotated auxiliary tree.
- [x] Replace the current boolean `_branches_interact()` helper with an interactor collector returning
  the six classes `A0`, `A1`, `A2`, `B0`, `B1`, and `B2`.


## Update rates based on interactions

- [x] Spot the places in the code where rates are updated or should be updated
- [x] plan the changes.
- [x] Check updated theoretical description to update rates, in  [Math and Algorithms.md](Math and Algorithms.md) section *Update rates based on gene-gene interactions*.
  - [x] Consider that the rates are a parameter for the process of sampling events given a distribution.
  - [x] Analyze how good is the new idea for loss rate update, identify pitfalls and advantages. Should we proceed with this?
  - [x] Tell me which would be a good value for `alpha` and `beta`.

- [x] Create a backup of current progress.
- [x] Plan refactor of code, we do not need to collect interactors. We still need the inverse maps.
- [x] Proceed with the implementation plan in the file  [implemented.md](implemented.md), section 'Refactor plan without interactor collection' (lines 261-274)
- [x] Remove the redundant `_gamma_star_at()` call inside `_host_system_summary()` and use a direct host-edge activity check instead.
- [x] Remove the redundant coexistence-plus-host-system double filtering in `_sample_recipient()` and derive `valid_species` directly from `host_system_edges` for host-system donors.
- [x] Update rate refreshing so it follows the new rationale `\Gamma \in E(T_A)` for affected auxiliary-tree branches and their host systems.
  - [x] Improve the theoretical description of how event times, event types, and the affected branches are drawn from the current rates together with the next fixed species-tree event.
  - [x] Keep the new rationale for `\Gamma`, i.e. treat it as the set of affected auxiliary-tree branches that determines which host systems should have their rates refreshed after an event.
  - [x] Ensure theoretical consistency across `The big picture`, `Where to update rates`, and the rate-update sections, and add a concrete implementation plan for changing the code from global refresh to `\Gamma`-restricted refresh.

## Improve `example_simulations.py` script

- [x] Separate rates

  In the docstring we have the parameter

  > \- ``symbiont_dtl_rates`` provides candidate ``(duplication, transfer, loss)``
  >
  > triples. Each triple is used for the symbiont simulation and then reused for
  >
  > the gene-tree simulation inside the resulting auxiliary tree.

  There should be two separated parameters; one for symbiont dtl rates, and one for gene dtl rates.

  Fix this in the code and the docs.

- [x] The script `example_simulations.py` should run standalone

  Since the input parameters are complex, let's use an input file to specify the simulations to perform.

  The input file should be easy to read/write by humans. I recommend something like this:

  > \> host_n_species
  >
  > 5, 10, 30
  >
  > \> symbiont_dtl_rates
  >
  > 0.133, 0.266, 0.266
  >
  > 0.3, 0.6, 0.6
  >
  > 0.6, 1.2, 1.2
  >
  > \> replace_probabilities
  >
  > 0.0, 0.5, 1.0
  >
  > \> n_simulations
  >
  > 5

  The only parameters with default values are:

  - alpha_values
  - beta_values
  - transfer_distance_biases
  - seed
  - output_dir (the default value for this is the date of the simulation and an unique integer to avoid colition)

- [x] Allow values for `alpha` and `beta` to be passed through `example_simulations.py`.

  These should be configurable from the standalone input file and forwarded to the gene-tree simulation inside the auxiliary tree.
