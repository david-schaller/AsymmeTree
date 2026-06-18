# Implemented

This file replies to the task list in `implementation.md` and records the current implementation
state.

## Review Theory And Code

### Theory review

The theory in `Math and Algorithms.md` is now consistent enough to map directly into code. The main
objects are all present in both places: the host-symbiont scenario, the auxiliary tree `T_A`, the
growing branch map `kappa`, the symbiont-only inverse map `gamma'`, the host-inclusive map
`gamma^*`, and the system-level rules for loss-rate refresh and transfer restriction.

The most important theory-to-code alignments are now in place:

- The code now mirrors the distinction between symbiont-only `gamma'` and host-inclusive `gamma^*`
  via `_gamma_prime_at()` and `_gamma_star_at()`.
- The loss-rate rule is implemented through host-system summaries with `M`, `N`, optional
  `crowding_cap`, and activation only when `N > 1`.
- The C0 transfer rule is implemented as a hard same-host-system filter before any optional
  distance-bias weighting is applied.
- The interactor classes remain useful for describing overlap structure, but they are no longer the
  primary mechanism for loss-rate updates.

### Code review

The prototype now covers the main theory path:

- `asymmetree.hologenome.hologenome_simulation.HologenomeSimulator` simulates symbiont trees inside
  host trees, creates one auxiliary tree per pair, simulates gene trees inside the auxiliary trees,
  and stores full/pruned as well as host/symbiont component trees.
- `create_auxiliary_tree()` follows the intended auxiliary-tree root shape
  `0_A -> rho_A -> {rho_H, 0_S}` and marks the symbiont planted root as a transfer event.
- `species2genes` in `holoevolve.genes` remains the active inverse of `kappa` for the current
  growing branches.
- `GeneTreeSimulator` builds `gamma_prime` from annotated auxiliary-tree edges, derives a
  host-inclusive live view through `_gamma_star_at()`, and keeps the six interactor classes as
  secondary/debugging structure.
- Each active gene branch now stores both `base_loss_rate` and effective `loss_rate`, and the
  simulator refreshes effective loss hazards from current host-system summaries after initialization,
  speciation, duplication, loss, HGT, and gene conversion.
- `_sample_recipient()` now filters coexisting recipient branches to the donor host system before
  applying additive or replacing transfer distance bias.

The remaining missing steps are:

- Add focused tests for auxiliary-tree shape, inverse maps, host-system rate refresh, and transfer
  restriction.
- Decide later whether the legacy interactor-based helpers should stay as diagnostics or be removed
  entirely from production code.

The hardest translation point remains that reconciliation edges are stored as node attributes and
often encoded by lower-node labels or label tuples. The normalization helpers introduced around
`gamma_prime` are still the key tool for keeping those comparisons consistent.

### AI Implementation Notes

The host-system rate refactor is now implemented in `holoevolve.genes`. A local `py_compile`
syntax check passed, but I could not run a runtime simulator smoke test here because `tralda` is
not installed in this workspace.

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

Rates are now initialized, sampled, or changed in these places:

- `_Branch.total_rate` still computes the branch-specific total event rate, while
  `_Branch.base_loss_rate` stores the inherited/base hazard separately from the current effective
  `loss_rate`.
- `simulate()` stores the base DTLG rates, validates `delta_l`, and now also validates/stores
  `alpha`, `beta`, and the optional `crowding_cap`.
- `_new_branch()` inherits `base_loss_rate` from the template branch, resets the effective
  `loss_rate` to that base value, and lets the global refresh recompute any host-system increment.
- `_host_system_summary()`, `_crowding_factor()`, `_refresh_loss_rates_for_host_system()`, and
  `_refresh_all_loss_rates()` implement the active host-system refresh path.
- `_run()` now refreshes loss hazards after initialization, each speciation event, and each sampled
  duplication, loss, HGT, or gene-conversion event.
- `_sample_recipient()` now uses `_gamma_star_at()` to enforce the C0 same-host-system restriction
  before optional additive or replacing transfer distance bias is applied.
- The legacy pairwise helpers `_collect_interactors()`, `_interaction_partners()`,
  `_increase_loss_rates_after_transfer()`, and `_decrease_loss_rates_after_loss()` remain in the
  file, but they are no longer part of the active rate-update path.

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

### Implemented refactor without interactor collection

1. `species2genes` remains the active inverse map `gamma`, and `gamma_prime` /
   `_gamma_prime_at()` remain the host-to-symbiont inverse maps.
2. `_host_system_summary()` now computes live host and symbiont copy counts together with `M` and
   `N`, using the current simulation time to decide which host-system edges are live.
3. `_refresh_loss_rates_for_host_system(host_edge)` and `_refresh_all_loss_rates()` reset every
   active branch to its `base_loss_rate` and then apply the host-system increment derived from the
   current state.
4. `_run()` now refreshes effective loss rates after initialization, every speciation, and every
   sampled duplication, loss, HGT, or gene-conversion event.
5. The legacy interactor helpers remain available for diagnostics, but the active rate path no
   longer depends on them.
6. `_gamma_star_at()` mirrors the host-inclusive live view needed for `\gamma_t^*`, and
   `_sample_recipient()` uses it to keep HGT recipients inside the donor host system before any
   optional distance-bias logic is applied.

### AI Implementation Notes

This implementation also adds user-facing `alpha`, `beta`, and `crowding_cap` parameters to
`holoevolve.dated_gene_tree()` / `GeneTreeSimulator.simulate()`, with defaults
`alpha = 0.25 * loss_rate` and `beta = 0.10 * loss_rate`. A later micro-cleanup removed the
redundant `_gamma_star_at()` call inside `_host_system_summary()`: host-edge activity is now
checked directly from the host branch interval, while `_gamma_prime_at()` is still called once to
get the live symbiont edges. Another small cleanup simplified `_sample_recipient()` so that, for
host and symbiont donors, recipient candidates are derived directly from `host_system_edges`
instead of first building a global coexisting set and then intersecting it with the same-host
constraint. Verification so far is limited to a local syntax check; I could not run an end-to-end
simulator smoke test here because `tralda` is not installed in this workspace.

### Gamma-restricted rate refresh

The theory notes now describe event sampling as the comparison between the next stochastic
birth-death waiting time and the next fixed event in the auxiliary tree. Loss-rate updates are
therefore phrased as hazard updates: changing an effective loss rate changes both the weighted
choice of the next stochastic event type and the exponential waiting-time parameter.

The active implementation now treats `Gamma` as a set of affected auxiliary-tree branches. After
initialization the simulator still performs one all-system refresh, but later events return only
the branches whose gene content changed:

- speciation returns the child auxiliary branches that became active,
- duplication and loss return the branch where the sampled gene branch lived,
- HGT returns donor and recipient branches, and
- gene conversion returns the empty set because the per-species copy count is unchanged.

The main loop maps these affected branches to their host systems, resets active genes in those
systems to their base loss rates, and reapplies the host-system crowding formula only there. This
replaces the previous global refresh after every event. Continuation branches now inherit the
template branch's current effective loss rate, which keeps gene-conversion branches consistent
when `Gamma` is empty.

### Example simulation script updates

`example_simulations.py` now has separate `symbiont_dtl_rates` and `gene_dtl_rates` inputs. Output
metadata records the two triples separately as `symbiont_*_rate` and `gene_*_rate` columns, and
scenario identifiers include both triples.

The script can now be run as a standalone command with a human-editable config file using
`> section` headers. Required sections are `host_n_species`, `symbiont_dtl_rates`,
`gene_dtl_rates`, `replace_probabilities`, and `n_simulations`; optional sections are
`transfer_distance_biases`, `seed`, and `output_dir`. A sample file was added at
`development notes/Example code/example_simulations_input.txt`. The script defers heavy simulation
imports until after CLI parsing, so `python3 example_simulations.py --help` works even in a partial
environment. An end-to-end simulation was not run here because `pandas` is not installed in this
workspace.

### Alpha and beta in example simulations

`example_simulations.py` now accepts `alpha_values` and `beta_values` as optional simulation-grid
parameters. Config files can provide floats or `default`/`None`; the default value keeps the
underlying `holoevolve` simulator behavior, where `alpha` and `beta` are derived from the gene
loss rate. Provided values are forwarded only to `simulate_gene_trees()`, because they control
host-system loss-rate refresh inside the auxiliary tree rather than symbiont-tree simulation.

### Task: TQDM progress bar for example simulations

Addressed. `example_simulations.py` now pre-computes the full Cartesian-product scenario count and
passes it as the `total` argument to `tqdm`. The progress bar is advanced once per completed
holobiont scenario. The `tqdm` runtime dependency was added to `pyproject.toml`; `uv.lock` was not
refreshed because `uv` is not installed in this workspace.

### Task: Projection maps in the theory notes

Addressed. `Math and Algorithms.md` now includes a section describing how the code obtains
host-side and symbiont-side projections. The description covers auxiliary-tree `level`
annotations, normalization of gene-node `reconc` values to lower auxiliary-node labels, retention
of nodes by host/symbiont level, suppression of connector nodes, synthetic roots, empty projected
components, and the fact that retained `reconc` values stay in the prefixed auxiliary-tree label
system.
