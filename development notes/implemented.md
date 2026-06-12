# Implemented

This file replies to the task list in `implementation.md` and records the current implementation
state.

## Review Theory And Code

### Theory review

The theory in `Math and Algorithms.md` is detailed enough to guide the next implementation step.
The main objects are present: the host-symbiont scenario, the auxiliary tree `T_A`, the growing
branch map `kappa`, inverse maps `gamma` and `gamma_prime`, the host-inclusive map `gamma^*`, and
the current system-level rules for loss-rate refresh and transfer restriction.

There are a few translation points to keep in mind before translating it directly into code:

- The theory now distinguishes the symbiont-only inverse map `gamma'` from the host-inclusive
  auxiliary map `gamma^*`. The next code pass should mirror that distinction explicitly instead of
  overloading one helper for both notions.
- The notation for edges should be normalized for the implementation. The code represents an edge
  by its lower node in most places, while the math alternates between edge pairs and lower-node
  shorthand.
- The interactor classes are now secondary: they still describe overlap structure, but the theory no
  longer uses them as the primary mechanism for loss-rate updates. The next code pass should keep
  them, if at all, as optional diagnostics rather than as the core rate path.
- The loss-rate rule should be translated into code-facing host-system summaries: per-host counts,
  `M`, `N`, activation only when `N > 1`, and an optional crowding cap kept in reserve.
- The transfer restriction for class C0 is now a hard host-system constraint in the theory, so the
  next code pass should filter recipient branches by host-system membership before any optional
  distance-bias weighting is applied.

### Code review

The prototype is now partly complete with respect to the theory:

- `asymmetree.hologenome.hologenome_simulation.HologenomeSimulator` simulates symbiont trees inside
  host trees, creates one auxiliary tree per pair, simulates gene trees inside the auxiliary trees,
  and stores full/pruned as well as host/symbiont component trees.
- `create_auxiliary_tree()` now follows the intended auxiliary-tree root shape
  `0_A -> rho_A -> {rho_H, 0_S}` and marks the symbiont planted root as a transfer event.
- `asymmetree.holoevolve.genes.GeneTreeSimulator` already has active gene branches with
  branch-specific DTLG rates, and `_get_branch_and_type()` samples from the branch-specific totals.
- `species2genes` in `holoevolve.genes` is already the active inverse of `kappa` for the current
  growing branches, although it is not named as `gamma`.
- `GeneTreeSimulator` now builds `gamma_prime` from annotated auxiliary-tree edges and collects
  active interactors in the six classes `A0`, `A1`, `A2`, `B0`, `B1`, and `B2`.

The remaining missing steps are:

- Add host-system summary helpers and refresh loss rates from current `gamma` / `gamma_prime`
  counts instead of from per-class interactor weights.
- Update transfer-recipient sampling with the C0 host-system restriction.
- Add focused tests for auxiliary-tree shape, inverse maps, host-system rate refresh, and transfer
  restriction. Interactor-class tests can remain as secondary structural checks.

The hardest translation point is that reconciliation edges are stored as node attributes and often
encoded by lower-node labels or label tuples. The implementation should introduce small helper
functions that normalize "edge-like" values before comparing them.

### AI Implementation Notes

The math is now coherent enough to implement directly. Edge normalization is already in place for
`gamma_prime` and interactor collection; the next code pass should reuse those helpers to build
host-system summaries, derive a host-inclusive view matching `gamma_t^*`, and enforce the C0
recipient filter.

## Correct Auxiliary Tree In Code

### Mismatch

The previous auxiliary tree was built in
`development notes/Example code/asymmetree/hologenome/hologenome_simulation.py`, lines 17-80,
and did not match the construction in `Math and Algorithms.md`, lines 101-110.

Observed mismatches:

- The code copies the entire host tree at line 48 and attaches `host_copy.root` under `AR` at line
  64. The theory requires the host subtree rooted at `rho_H`, not the host planted root `0_H`.
- The code attaches both copied planted roots under `AR`. The theory requires `AR -> rho_H` and
  `AR -> 0_S`.
- The code sets `AR.tstamp` to `max(host_copy.root.tstamp, symbiont_copy.root.tstamp)` at line 54.
  The theory requires `tau_A(rho_A) = tau_A(0_S) + epsilon_0`.
- The code sets `A0.tstamp = merged_tstamp + root_offset` at line 71. The theory requires a second
  offset, `tau_A(0_A) = tau_A(0_S) + epsilon_0 + epsilon_1`.
- The code maps `AR` and `A0` to themselves, but it does not enforce `mu_A(0_S) = rho_A rho_H`.
- The code does not convert `0_S` into a transfer event with a transferred outgoing edge.

### Implemented upgrade

- `create_auxiliary_tree()` now accepts `root_offset` for `0_S -> rho_A` and
  `planted_root_offset` for `rho_A -> 0_A`; if the second offset is omitted, it reuses
  `root_offset`.
- The host component is copied from the canonical host root by `_extract_host_component_root()`,
  so the old host planted root is no longer attached below `rho_A`.
- The symbiont tree remains rooted at its planted root `0_S`.
- The symbiont planted root is set to event `H`, reconciled to `("AR", rho_H)`, and its child edge
  is marked as transferred.
- Distances are recomputed with `distance_from_timing()` after the new root structure is assembled.
- `holoevolve.genes._speciation()` now preserves auxiliary transfer and loss event labels when
  creating the corresponding gene-tree node.

### AI Implementation Notes

The auxiliary-tree fix is complete. The root shape and transfer annotation were verified with a
small host/symbiont/auxiliary/gene simulation.

## Collect Interactors During Gene Tree Simulation

### Construct inverse maps

Use two inverse maps during `holoevolve.genes.GeneTreeSimulator._run()`:

- `gamma`: active auxiliary edge to active gene branches. This already exists as
  `self.species2genes`, initialized at lines 448-449 and updated in `_new_branch()` /
  `_remove_branch()`.
- `gamma_prime`: host auxiliary edge to symbiont auxiliary edges. This is now built in
  `GeneTreeSimulator.__init__()` by `_build_auxiliary_inverse_maps()`.

Implemented `gamma_prime` details:

- `host_edges_by_key` maps normalized host lower-node labels to host edge nodes.
- `symbiont_edges_by_key` maps normalized symbiont lower-node labels to symbiont edge nodes.
- `gamma_prime` maps host edge nodes to the symbiont edge nodes reconciled into them.
- `_gamma_prime_at(host_edge, tstamp)` provides the current time-restricted symbiont-only inverse map.

### Collect interactors

`_branches_interact()` has been replaced by `_collect_interactors()`, which returns the six
interaction classes:

```python
def _collect_interactors(self, branch: _Branch) -> dict[str, set[_Branch]]:
    interactors = {key: set() for key in ("A0", "A1", "A2", "B0", "B1", "B2")}
    s = branch.S_edge
    h = self._host_edge_for_species_edge(s)

    if self._branch_level(branch) == "host":
        interactors["A0"] = set(self.species2genes[s]) - {branch}
        for symbiont_edge in self.gamma_prime.get(s, set()):
            interactors["B0"].update(self.species2genes[symbiont_edge])
    elif self._branch_level(branch) == "symbiont":
        interactors["A1"] = set(self.species2genes[s]) - {branch}
        interactors["B1"] = set(self.species2genes[h]) - {branch}
        for symbiont_edge in self.gamma_prime.get(h, set()) - {s}:
            interactors["B2"].update(self.species2genes[symbiont_edge])
    else:
        interactors["A2"] = set(self.species2genes[s]) - {branch}

    for partners in interactors.values():
        partners.discard(branch)
    return interactors
```

The existing `_interaction_partners()` API now returns the union of these sets so the current
loss-rate hooks continue to work.

### AI Implementation Notes

The interactor collection step is complete. These classes can stay as debugging or structural aids,
but the planned rate refactor should no longer depend on per-class weights.

## Update Rates Based On Interactions

### Current rate update locations

Rates are currently initialized, sampled, or changed in these places:

- `_Branch.total_rate`, lines 54-57, computes branch-specific total event rate.
- `simulate()`, lines 162-218, stores the base DTLG rates and validates `delta_l`.
- `_new_branch()`, lines 247-285, initializes or inherits branch-specific rates.
- `_get_branch_and_type()`, lines 323-347, samples events from current branch-specific rates.
- `_run()`, lines 454-483, samples event times from the sum of current branch totals.
- `_increase_loss_rates_after_transfer()`, lines 414-418, exists but is not currently called.
- `_decrease_loss_rates_after_loss()`, lines 420-426, divides surviving partner loss rates.
- `_speciation()`, lines 515-519, contains a commented-out transfer-inheritance hook.
- `_loss()`, lines 587-599, collects partners before removal and decreases their loss rates.
- `_sample_recipient()`, lines 623-691, chooses transfer recipients and is the right place for the
  C0 host-system restriction.

### Reading the revised theory

The revised section in `Math and Algorithms.md` is an improvement over the earlier class-weight idea
for loss updates. Instead of trying to weight `A0` through `B2` separately, it uses occupancy
summaries for one host system at a time:

- For a symbiont branch `s` inside host branch `h`, update each active gene on `s` with
  `r_l^* = r_l + alpha * |gamma(s)| * N / M`.
- For the host branch `h`, update each active host-resident gene with
  `r_l^* = r_l + beta * |gamma(h)| * N / M`.

The theory now uses the corrected host-side formula `|gamma(h)|`, so the host loss rate depends on
the host copy count rather than on one symbiont copy count.

### Rates are hazard parameters

This point is important for the implementation: these rates are event intensities, not direct
per-event probabilities. Changing `loss_rate` affects the simulation in two ways at once:

- it changes the branch-local chance that the next event on that branch is a loss rather than D/H/GC;
- it changes the total process intensity because `_run()` samples the next waiting time from the
  sum of all active branch rates.

So the new update rule should be implemented as a recomputation of current branch hazards from the
current active state, not as a one-off additive probability correction.

### Advantages

- It uses only the inverse maps that the simulator already needs: `gamma = species2genes` and
  `gamma_prime`.
- It removes six user-facing interaction-weight parameters and replaces them with just `alpha` and
  `beta`.
- It is easier to recompute from scratch after state changes than the earlier pairwise
  multiply/divide idea.
- The normalization by `M / N` naturally compares each species to the current average gene load in
  the host system.

### Pitfalls and safeguards

- If the formula is applied when `N = 1`, a lone host lineage would still receive an interaction
  penalty even though no holobiont interaction exists. The update should therefore activate only
  for true holobiont systems with at least one live symbiont edge.
- The extra loss hazard in one species grows like `n_s^2 * N / M` at the system level, so large
  duplication bursts can make losses dominate too quickly. A cap such as
  `crowding = min(|gamma(.)| * N / M, c_max)` with `c_max` around `2` or `3` is worth keeping in
  reserve.
- Gene conversion, replacing HGT, species loss, and symbiont loss can change many active counts at
  once. Trying to undo old multipliers locally will be fragile; refreshing from scratch by host
  system is safer.
- Independent auxiliary branches should stay at the base loss rate because they are outside any
  host-symbiont system.

### Should we proceed?

Yes, with those safeguards. The revised formula is closer to the current code structure than the
old interactor-weight approach, and it gives a cleaner first implementation target.

### Suggested values for alpha and beta

It is better to choose `alpha` and `beta` as fractions of the base loss rate `r_l` rather than as
fixed absolute constants. That keeps the effect comparable across the current example parameter
sets, where `loss_rate` ranges from `0.266` to `1.2`.

Recommended first pass:

- Host-keeps regime: `alpha = 0.25 * r_l`, `beta = 0.10 * r_l`.
- Neutral baseline: `alpha = beta = 0.15 * r_l`.
- Stronger exploratory regime: `alpha = 0.50 * r_l`, `beta = 0.20 * r_l`.

The first option is the safest default here because it makes symbiont genes somewhat more
expendable without overwhelming the existing D/T/H dynamics. Around average occupancy
(`|gamma(.)| * N / M ~= 1`), it increases symbiont losses by about `25%` of `r_l` and host losses
by about `10%` of `r_l`.

### Backup

A pre-refactor snapshot was created in
`development notes/Archive/backup_2026-06-07_pre_rate_refactor/`. It contains the current theory
notes, task files, and the main hologenome code paths that will be affected by the next
implementation pass.

### Refactor plan without interactor collection

1. Keep `species2genes` as the active inverse map `gamma` and keep `gamma_prime` /
   `_gamma_prime_at()` as the host-to-symbiont inverse maps.
2. Add a helper that summarizes one host system: live host copy count, live symbiont copy counts,
   `M`, and `N`.
3. Add `_refresh_loss_rates_for_host_system(host_edge)` and `_refresh_all_loss_rates()` that reset effective loss rates to the inherited or base value and then apply the new system-level increment.
4. Refresh rates after every event that changes multiplicities or placements: duplication, loss,
   HGT, gene conversion, speciation, and species or symbiont loss. HGT may require refreshing both
   the donor and recipient host systems.
5. Remove `_collect_interactors()`, `_interaction_partners()`,
   `_increase_loss_rates_after_transfer()`, and `_decrease_loss_rates_after_loss()` from the rate path. They can remain temporarily only if they are still useful for debugging or for checking the
   C0 host-transfer condition.
6. Add a host-inclusive helper mirroring `\gamma_t^*` and use it in `_sample_recipient()` so that recipient branches are filtered to the donor host system before any optional distance-bias logic is applied.

### AI Implementation Notes

The updated theory now points to a cleaner refactor: interaction classes were useful for exploring
the structure of overlaps, but they are no longer the right primitive for loss-rate updates. The
next coding pass should center the rate logic on host-system summaries derived from `gamma` and
`gamma_prime`, and should enforce the C0 host-system restriction directly in `_sample_recipient()`.
