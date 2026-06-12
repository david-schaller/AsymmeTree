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
- [ ] Proceed with the implementation plan in the file  [implemented.md](implemented.md), section 'Refactor plan without interactor collection' (lines 261-274)

