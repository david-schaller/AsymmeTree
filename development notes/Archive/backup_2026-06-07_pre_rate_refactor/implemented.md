# Implemented

This file replies to the task list in `implementation.md` and records the current implementation
state.

## Review Theory And Code

### Theory review

The theory in `Math and Algorithms.md` is detailed enough to guide the next implementation step.
The main objects are present: the host-symbiont scenario, the auxiliary tree `T_A`, the growing
branch map `kappa`, inverse maps `gamma` and `gamma_prime`, the six interactor classes, and the
rate changes based on those classes.

There are a few consistency issues to fix before translating it directly into code:

- The time-restricted inverse map is written as `gamma'_t(ab) = {uv in gamma(b) ...}`. This should
  refer to `gamma_prime(ab)` or to an explicitly defined lower-node key for the host edge.
- The notation for edges should be normalized for the implementation. The code represents an edge
  by its lower node in most places, while the math alternates between edge pairs and lower-node
  shorthand.
- In the interactor algorithm, the symbiont case assigns `A1 <- gamma(s)` but the table defines all
  classes as excluding the focal branch `g`. The implementation should remove `g` from every class.
- The condition `s <= 0_S` should be restated as "the lower node of `s` is in the symbiont
  component" or `level == "symbiont"`. As written, it is easy to misread for rooted-tree order.
- The rate section says the six classes should be weighted differently, but the exact combination
  rule is not specified. The code needs an explicit rule such as multiplicative weights per partner,
  additive increments, or a capped transform.
- The `S-keeps` / `H-keeps` asymmetry is described qualitatively. It should be converted into
  parameters that state which B-class multipliers apply to host-resident and symbiont-resident
  genes.
- The transfer-probability decrease for class C0 needs a numeric multiplier and a rule for combining
  it with the existing distance-bias weights.

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

- Update loss rates from the collected interactor classes using explicit per-class weights.
- Update transfer-recipient sampling with the C0 host-system penalty.
- Add focused tests for auxiliary-tree shape, inverse maps, interactor classes, and rate updates.

The hardest translation point is that reconciliation edges are stored as node attributes and often
encoded by lower-node labels or label tuples. The implementation should introduce small helper
functions that normalize "edge-like" values before comparing them.

### AI Implementation Notes

The math is usable. Edge normalization is now implemented for `gamma_prime` and interactor
collection; the next code pass should reuse those helpers for weighted rate updates and transfer
recipient penalties.

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
- `_gamma_prime_at(host_edge, tstamp)` provides the time-restricted inverse map.

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

The interactor collection step is complete. Differential weights and the `S-keeps` / `H-keeps`
asymmetry are still open in the later rate-update task.

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
- `_sample_recipient()`, lines 623-691, chooses transfer recipients and is the right place for C0
  transfer weighting.

### Plan for loss rates

1. Add a base loss-rate field to `_Branch`, or otherwise store the branch's non-interaction loss
   rate separately from the effective loss rate.
2. Replace incremental multiply/divide logic with a recomputation helper:
   `_refresh_interaction_loss_rates(branches: Iterable[_Branch] | None = None)`.
3. For each active branch, call `_collect_interactors(branch)` and compute an effective multiplier
   from the six class sizes and user-provided class weights.
4. Add user-facing parameters for the six weights, for example:
   `loss_interaction_multipliers={"A0": 2.0, "A1": 2.0, "A2": 1.0, "B0": 2.0, "B1": 2.0, "B2": 2.0}`.
5. Add an asymmetry mode for B0/B1, for example `host_symbiont_loss_bias in (None, "S-keeps",
   "H-keeps")`, and apply it inside the multiplier calculation.
6. Refresh rates after any event that changes active branch placement or multiplicity:
   speciation, duplication, loss, HGT, gene conversion, and species/symbiont loss.
7. In `_speciation()`, specifically refresh after a new branch is inherited onto a transferred
   symbiont edge (`getattr(S_w, "transferred", 0)`), because this is the symbiont-transfer trigger
   called out in the old implementation notes.
8. Keep the minimum effective loss rate at the branch's base loss rate.

The recomputation approach is safer than pairwise divide-on-loss because it avoids stale multipliers
after duplications, speciations, replacing HGT, gene conversion, or mass removal at species-loss
nodes.

### Plan for transfer rates

1. Keep `_coexisting_species_edges()` as the time-consistency filter.
2. In `_sample_recipient()`, compute the donor host `z0 = host_edge(branch.S_edge)`.
3. Compute `same_host_edges = gamma_prime_at(z0, event_tstamp)` plus the host edge itself.
4. When additive HGT samples a recipient edge, multiply existing distance-bias weights by a C0
   penalty if the candidate is outside `same_host_edges`.
5. When replacing HGT samples recipient gene branches, apply the same C0 penalty through each
   candidate gene branch's `S_edge`.
6. Add a parameter such as `inter_host_transfer_penalty`, defaulting to `1.0` to preserve current
   behavior. Values below `1.0` decrease transfers between different host systems.

### AI Implementation Notes

Rate updates should be driven by current active interactor sets rather than only by the event that
created them. This matches the theory more closely and makes future differential weights much less
fragile.
