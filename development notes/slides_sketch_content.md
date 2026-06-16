# Sketch content description

This file is a first-pass interpretation of the hand-drawn material in
[draft_figures.pdf](draft_figures.pdf).

The draft is handwritten and lightly colored, so the section numbering below
follows reading order rather than any uncertain handwritten label:

1. page 1, upper band
2. page 1, middle band
3. page 1, lower band
4. page 2, upper band
5. page 2, middle band
6. page 2, lower band
7. page 3, full-page composite

The mathematical identifications below are working inferences based on the
sketch and the formalism in [Math and Algorithms.md](Math and Algorithms.md).
They should be confirmed again during render/check against the PDF.

## Shared visual conventions

- Main geometry is drawn in dark graphite or black pencil.
- The first page uses blue-cyan accents.
- The second page uses orange-ochre accents.
- The third page uses magenta accents.
- All sections are informal pencil sketches, but they are already structured as
  scientific diagrams: trees, correspondences, highlighted branches, and
  comparison layouts.

## Section 1

### Mathematical objects shown

- Likely a species tree $T_S$ and a gene tree $T_G$ in the embedded format.
- The relevant formal object is a relaxed evolutionary scenario
  $(T_S,T_G,t,\sigma,\mu)$, with the visual nesting standing in for the
  reconciliation map $\mu: V(T_G) \to V(T_S) \cup E(T_S)$.

### Content

- Location: page 1, upper band.
- Likely role: first and simplest example of the "tree embedded in tree"
  representation.
- Content: the larger graphite scaffold should be read as the ambient species
  tree $T_S$, while the smaller blue-cyan structure represents the embedded gene
  tree $T_G$. The nesting suggests the reconciliation map
  $\mu: V(T_G) \to V(T_S) \cup E(T_S)$, with only a few emphasized branches or
  placements highlighted in color so the basic two-level scenario remains easy
  to parse.
- Pencil/style cues: graphite carries the full outline, while blue-cyan is used
  only for the branches or regions of $T_G$ that should stand out.

## Section 2

### Mathematical objects shown

- Again most likely a species tree $T_S$ together with a gene tree $T_G$.
- Compared with Section 1, this sketch seems more explicit about internal event
  structure, so the natural formal reading is a relaxed or dated evolutionary
  scenario $(T_S,T_G,t,\sigma,\mu,\tau_S,\tau_G)$, where $t$ marks gene-level
  events and $\tau_S, \tau_G$ encode temporal placement.

### Content

- Location: page 1, middle band.
- Likely role: denser embedded scenario, still in the nested-tree format.
- Content: this looks like a richer realization of the same pair $(T_S,T_G)$,
  now with more internal event points that can be read as vertices of $T_G$
  labeled by $t$, and with a more explicit temporal feel consistent with
  $\tau_S$ and $\tau_G$. The graphite drawing still supplies most of $T_S$ and
  the full scaffold, while blue-cyan marks selected branches, local clusters, or
  paths of $T_G$ that are important in the dated scenario.
- Pencil/style cues: graphite gives the dense base structure, and blue-cyan is
  reserved for emphasized routes inside the embedded gene history.

## Section 3

### Mathematical objects shown

- Most likely still a two-level pair $(T_S,T_G)$, but arranged as a focused or
  comparative view rather than a full overview.
- The best mathematical reading is a restricted or emphasized subscenario, for
  example subtrees such as $T_S|_{L'}$ and $T_G|_{L'}$ together with the
  inherited maps $\sigma$ and $\mu$.

### Content

- Location: page 1, lower band.
- Likely role: final blue-page variant of the embedded representation.
- Content: the compartmentalized layout suggests a focused subscenario in which
  a restricted host/species part $T_S|_{L'}$ is compared with or aligned to a
  restricted gene part $T_G|_{L'}$, still governed by the inherited maps
  $\sigma$ and $\mu$. The graphite blocks define the ambient comparison space,
  while the blue-cyan accents concentrate on the substructure that has been
  selected for emphasis.
- Pencil/style cues: the same graphite plus blue-cyan pairing as Sections 1 and
  2 is kept, but the color is more localized around the singled-out subscenario.

## Section 4

### Mathematical objects shown

- This is most naturally read as two rooted phylogenetic trees shown separately
  and related by a reconciliation-type map.
- In the notation of the project, the strongest current guess is a host tree
  $T_H$ and a symbiont tree $T_S$ linked by $\mu'$, i.e. the backbone of a
  holobiont scenario $(T_H,T_S,t',\sigma',\mu',\tau_H,\tau_S)$.

### Content

- Location: page 2, upper band.
- Likely role: first example of the "two trees related by a map"
  representation.
- Content: the separated graphite trees are best read as $T_H$ and $T_S$ shown
  side by side, with the relation between them standing for the map $\mu'$ in a
  holobiont scenario. Orange-ochre highlights seem to identify the key
  branches, nodes, or correspondences where the placement of $T_S$ inside $T_H$
  should be visually tracked.
- Pencil/style cues: graphite outlines both trees, while orange-ochre is used
  only on the main mapped elements.

## Section 5

### Mathematical objects shown

- Likely a second, more explicit view of the same kind of mapped-tree data as in
  Section 4.
- In formal terms this is still best read as a holobiont scenario
  $(T_H,T_S,t',\sigma',\mu',\tau_H,\tau_S)$, now with more visible emphasis on
  where branches or nodes of $T_S$ are placed under the map $\mu'$ and how the
  leaves are related by $\sigma'$.

### Content

- Location: page 2, middle band.
- Likely role: larger and clearer mapped-tree example.
- Content: this appears to expand the mapped pair $(T_H,T_S)$ into a clearer
  reconciliation view in which the placements given by $\mu'$ and the leaf-level
  relation $\sigma'$ are easier to read across the page. The two graphite
  branching objects remain distinct, while the stronger orange-ochre marks pick
  out the branches, links, or assignments that the viewer should follow first.
- Pencil/style cues: orange-ochre is more prominent than in Section 4, but it
  still functions as selective emphasis over a graphite base drawing.

## Section 6

### Mathematical objects shown

- This is the strongest candidate among the page-2 panels for introducing the
  auxiliary-tree side of the construction.
- The likely objects are the auxiliary tree $T_A$ derived from the host and
  symbiont trees, together with the inherited host map $\mu'$ and possibly the
  auxiliary placement map $\mu_A$. If the dense annotations correspond to the
  simulation workflow, this panel may also be hinting at $\Pi$, $\kappa$,
  $\gamma$, $\gamma'$, and $\gamma^*$.

### Content

- Location: page 2, lower band.
- Likely role: most detailed orange-page example in the mapped representation.
- Content: the dense layout makes this the best candidate for a view centered on
  the auxiliary tree $T_A$, with nearby annotations hinting at the inherited map
  $\mu'$, the auxiliary placement map $\mu_A$, or even the simulation-side
  objects $\Pi$, $\kappa$, $\gamma$, $\gamma'$, and $\gamma^*$. Most branches
  and local structure are drawn in graphite, while orange-ochre calls attention
  to the active branches, correspondences, or directional cues that organize the
  construction.
- Pencil/style cues: graphite carries most of the detail, while orange-ochre is
  reserved for the branches, nodes, or arrows that should drive the viewer's
  attention.

## Section 7

### Mathematical objects shown

- This composite page is best read as the full three-level construction.
- The central formal object is the hologenome scenario
  $(T_H,T_S,T_G,t,\sigma,\mu,t',\sigma',\mu',\tau_H,\tau_S,\tau_G)$, likely
  accompanied by the auxiliary tree $T_A$ used for the actual gene simulation.
- Because the page is composite, it is also the best candidate for displaying
  the simulation-side maps and sets together: the growing branches $\Pi$, the
  growing map $\kappa$, the inverse maps $\gamma$ and $\gamma'$, and the
  host-system map $\gamma^*$, possibly with their time-restricted forms.

### Content

- Location: page 3, full-page composite.
- Likely role: summary figure for the new three-level extension rather than one
  small isolated panel.
- Content: the full-page composition is most naturally read as one overview of
  the hologenome scenario
  $(T_H,T_S,T_G,t,\sigma,\mu,t',\sigma',\mu',\tau_H,\tau_S,\tau_G)$ together
  with the auxiliary construction around $T_A$. Its internal zones likely
  separate the tree layers from the simulation-side objects such as $\Pi$,
  $\kappa$, $\gamma$, $\gamma'$, and $\gamma^*$, while the magenta accents mark
  the structures that are new or conceptually central in the three-level
  extension.
- Pencil/style cues: graphite carries the full branching geometry and most
  annotations, while magenta is reserved for the most important composite
  structures.
- SVG implication: when this figure is redrawn, it should remain one composite
  explanatory figure with internal subareas, not be flattened into unrelated
  mini-panels.

## Notes for later SVG work

- Preserve the page-wise accent-color families: blue-cyan, orange-ochre,
  magenta.
- Keep the reading order and the "nested versus mapped" contrast explicit.
- Keep the mathematical-object subsection at the top of each section, since it
  gives the design target before visual polish.
- During render/check, verify the exact boundaries and any handwritten numbering
  directly against the PDF before final polishing.
