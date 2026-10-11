"""Generate the figures shown in docs/guide/visualization.md.

The images are written to examples/results/visualization/, which is not tracked by git. Copy this
directory to the corresponding location in the documentation (e.g., docs/images/visualization).
"""

import random
from pathlib import Path

import numpy as np
from tralda.visualization import NodeStyle, TreeStyle

import asymmetree.treeevolve as te
from asymmetree.visualization import visualize, assign_colors


SEED = 0

OUTPUT_DIR = Path(__file__).parent / "results" / "visualization"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

random.seed(SEED)
np.random.seed(SEED)


# --------------------------------------------------------------------------
#                              SIMULATION
# --------------------------------------------------------------------------

S = te.species_tree_n_age(6, 1.0)

dated_gene_tree = te.dated_gene_tree(S, dupl_rate=0.3, loss_rate=0.3, hgt_rate=0.3)
full_gene_tree = te.rate_heterogeneity(
    dated_gene_tree,
    S,
    base_rate=1.0,
    autocorr_variance=0.1,
    rate_increase=("gamma", 0.5, 2.2),
    inplace=False,
)
pruned_gene_tree = te.prune_losses(full_gene_tree)

species_colors, gene_colors = assign_colors(S, dated_gene_tree)


# --------------------------------------------------------------------------
#               HORIZONTAL LAYOUT WITH DISTANCE AXIS (DEFAULT)
# --------------------------------------------------------------------------

visualize(
    S,
    color_dict=species_colors,
    path=OUTPUT_DIR / "species_tree.png",
    distance_axis=True,
)
visualize(
    dated_gene_tree,
    color_dict=gene_colors,
    path=OUTPUT_DIR / "gene_tree_dated.png",
    distance_axis=True,
)
visualize(
    full_gene_tree,
    color_dict=gene_colors,
    path=OUTPUT_DIR / "gene_tree_full.png",
    distance_axis=True,
)
visualize(
    pruned_gene_tree,
    color_dict=gene_colors,
    path=OUTPUT_DIR / "gene_tree_pruned.png",
    distance_axis=True,
)


# --------------------------------------------------------------------------
#                       OTHER LAYOUTS AND EDGE LENGTHS
# --------------------------------------------------------------------------

# circular layout, all leaves end at the same distance from the root
visualize(
    pruned_gene_tree,
    color_dict=gene_colors,
    path=OUTPUT_DIR / "gene_tree_circular_even.png",
    layout_mode="circular",
    edge_length_mode="even",
)

# vertical layout, nodes are placed according to their topological depth
visualize(
    pruned_gene_tree,
    color_dict=gene_colors,
    path=OUTPUT_DIR / "gene_tree_vertical_rank.png",
    layout_mode="vertical",
    edge_length_mode="rank",
)

# horizontal layout, all edges have length 1
visualize(
    pruned_gene_tree,
    color_dict=gene_colors,
    path=OUTPUT_DIR / "gene_tree_uniform.png",
    edge_length_mode="uniform",
)


# --------------------------------------------------------------------------
#                           LABELS AND STYLE
# --------------------------------------------------------------------------

# the species of each gene is appended to the leaf label
visualize(
    pruned_gene_tree,
    color_dict=gene_colors,
    path=OUTPUT_DIR / "gene_tree_species_info.png",
    species_info=True,
    distance_axis=True,
)

# custom style: other symbols and thicker edges (leaf colors have to be part of the style)
custom_style = TreeStyle.from_maps(
    symbol_attr="event",
    symbol_map={"D": "diamond", "H": "triangle-right", "GC": "hexagon", "L": "x-mark"},
    node_color_attr="label",
    node_color_map=gene_colors,
    root_symbol="circle-inner-ring",
    leaf_symbol="circle",
    internal_symbol="circle",
    default=NodeStyle(edge_lw=2.0),
)
visualize(
    pruned_gene_tree,
    tree_style=custom_style,
    path=OUTPUT_DIR / "gene_tree_custom_style.png",
    distance_axis=True,
    scale_symbols=1.3,
    fontsize="large",
)
