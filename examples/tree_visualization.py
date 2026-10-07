from pathlib import Path

import asymmetree.treeevolve as te
from asymmetree.visualization import visualize, assign_colors


OUTPUT_DIR = Path(__file__).parent


# --------------------------------------------------------------------------
#                            SPECIES TREE
# --------------------------------------------------------------------------

S = te.species_tree_n_age(
    6,
    1.0,
    model="BDP",
    innovation=True,
    birth_rate=1.0,
    death_rate=0.5,
    contraction_probability=0.2,
)


# --------------------------------------------------------------------------
#                             GENE TREE
# --------------------------------------------------------------------------

T = te.dated_gene_tree(
    S,
    dupl_rate=0.7,
    loss_rate=0.7,
    hgt_rate=0.7,
    gc_rate=0.7,
    dupl_polytomy=0.5,
    replace_prob=0.5,
    transfer_distance_bias="inverse",
)


# --------------------------------------------------------------------------
#                         RATE HETEROGENEITY
# --------------------------------------------------------------------------

te.rate_heterogeneity(T, S, inplace=True)


# --------------------------------------------------------------------------
#                           VISUALIZATION
# --------------------------------------------------------------------------

species_colors, gene_colors = assign_colors(S, T)

visualize(S, color_dict=species_colors, path=OUTPUT_DIR / "testfile_speciestree.pdf")
visualize(T, color_dict=gene_colors, path=OUTPUT_DIR / "testfile_genetree.pdf")
visualize(
    T,
    color_dict=gene_colors,
    path=OUTPUT_DIR / "testfile_genetree_circular.pdf",
    layout_mode="circular",
    edge_length_mode="even",
)
