# Sketch hand-style description

This file describes the visual style that the SVG reconstructions should follow
when translating [draft_figures.pdf](draft_figures.pdf) into presentation-ready
figures.

The target is not a literal scan-like imitation of the sketch, and not a fully
generic polished vector diagram either. The target style is a controlled
scientific redraw: formal enough for a research talk, but still recognizably
close to the author's natural hand style.

## Core design goal

- Final figures must be clear enough for a scientific presentation.
- Final figures must still feel like they originate from a hand-drawn draft.
- The redraw should therefore keep mild human irregularity while removing
  accidental messiness.

## What to preserve from the hand style

- Slight curvature in lines and branches instead of perfectly rigid geometry.
- Mild asymmetry in branching and spacing, so the figures do not feel generated
  from a strict template.
- Uneven visual rhythm in local details, especially where branches split or
  correspondences are marked.
- Selective use of accent color as if added by hand after the graphite
  structure.
- A restrained sketch feeling: expressive, but not cartoonish.

## What to formalize for slides

- Branches and connectors must remain crisp and readable at presentation scale.
- Relative spacing between trees, labels, and maps must be cleaner than in the
  draft.
- Event marks, arrows, and correspondences must be unambiguous.
- Color emphasis must support interpretation, not decoration.
- Visual clutter from the raw sketch should be reduced during SVG design.

## Stroke language

- Use dark graphite or near-black as the main structural stroke color.
- Keep stroke widths mostly consistent, with only slight variation where a
  natural hand-drawn feel is useful.
- Avoid perfectly mechanical straightness; slight bends or soft deviations are
  desirable when they do not reduce clarity.
- Avoid noisy scribble textures. The hand style should appear deliberate, not
  shaky.

## Color language

- Preserve the draft's page-wise accent families:
  blue-cyan, orange-ochre, and magenta.
- Accent colors should be applied selectively to the branches, nodes, mapped
  relations, or substructures that deserve attention.
- The graphite structure should remain dominant, with color acting as emphasis.
- Avoid saturated neon colors or heavy fills that would overpower the formal
  content.

## Geometry and composition

- Tree layouts should remain mathematically legible first.
- Embedded figures should clearly show one structure living inside or along the
  other.
- Mapped figures should clearly separate the trees and make the correspondence
  legible.
- Composite figures may retain internal zoning, but each zone should feel part
  of a single coherent plate.
- Minor irregularity is welcome, but not if it obscures the mathematical
  objects.

## Labels and annotations

- Text and symbols should be cleaner than the handwritten draft.
- Mathematical notation should be typeset clearly and consistently.
- Label placement may keep a slightly hand-placed feel, but should avoid
  collisions and ambiguity.
- Annotations should look integrated with the figure rather than pasted on top.

## Explicit balance to maintain

- Too rough is wrong:
  if the SVG looks like a notebook scan, it is not formal enough.
- Too clean is also wrong:
  if the SVG looks like a generic software-generated phylogeny, it has lost the
  author's hand style.
- The correct balance is a polished hand-drawn scientific diagram.

## Practical SVG guidance

- Start from a clean vector scaffold with correct mathematical structure.
- Reintroduce hand-style character through small curvature, spacing choices, and
  restrained accent-color placement.
- Prefer consistency across the figure set, so all seven figures feel like part
  of one talk.
- When uncertain, choose clarity first, then re-add only subtle hand-made
  character.
