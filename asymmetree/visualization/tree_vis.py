"""Species and gene tree visualization."""

from __future__ import annotations

import os
from functools import partial
from typing import Any, Literal, TYPE_CHECKING

import numpy as np
from tralda.datastructures import Tree
from tralda.datastructures import TreeNode
from tralda.visualization import LayoutMode
from tralda.visualization import NodeStyle
from tralda.visualization import plot_tree
from tralda.visualization import TreeStyle

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from plotly.graph_objects import Figure as PlotlyFigure


def assign_species_colors(tree: Tree) -> dict:
    """Assign a unique color to each (extant) species in a species tree.

    Args:
        tree: The species tree with uniquely labeled (non-loss) leaves.

    Returns:
        A dictionary containing the labels of the extant species as keys and the assigned colors as
            values.
    """
    return _assign_cmap({v.label for v in tree.leaves() if v.event != "L"})


def assign_gene_colors(tree: Tree, species_colors: dict | None = None) -> dict:
    """Assign a color to each (extant) gene in a gene tree according to its species.

    Args:
        tree: The gene tree whose (non-loss) leaves have the 'reconc' attribute, i.e, the
            information to which species they belong.
        species_colors: A dictionary containing the values of the 'color' attribute appearing in the
            tree as keys and assigned colors as values.

    Returns:
        A dictionary containing the labels of the extant genes as keys and the assigned colors as
            values.
    """
    if species_colors is None:
        species_colors = _assign_cmap({v.reconc for v in tree.leaves() if v.event != "L"})

    return {v.label: species_colors[v.reconc] for v in tree.leaves() if v.event != "L"}


def assign_colors(species_tree: Tree, gene_tree: Tree) -> tuple[dict, dict]:
    """Assign a color to each (extant) species and gene.

    Args:
        species_tree: The species tree with uniquely labeled (non-loss) leaves.
        gene_tree: The gene tree whose (non-loss) leaves have the 'reconc' attribute, i.e, the
            information to which species they belong, that furthermore appear as labels of the
            (non-loss) leaves in the species tree.

    Returns:
        A `dict` containing the labels of the extant species as keys and the assigned colors as
            values.
        A `dict` containing the labels of the extant genes as keys and the assigned colors as
            values.
    """
    species_colors = assign_species_colors(species_tree)

    return species_colors, assign_gene_colors(gene_tree, species_colors=species_colors)


def _assign_cmap(colors: set) -> dict:
    """Assign a color to each label in the given set of labels.

    Args:
        colors: A set of labels for which colors should be assigned.

    Returns:
        A dictionary containing the labels as keys and the assigned colors as values.
    """
    color_dict = {}

    try:
        from matplotlib import colormaps
    except ImportError as exc:
        raise ImportError(
            "Color assignment requires matplotlib. Install it with: pip install asymmetree[matplotlib]"
        ) from exc

    if len(colors) <= 10:
        cmap = colormaps["tab10"](np.arange(len(colors), dtype=int))
    else:
        cmap = colormaps["jet"](np.linspace(0, 1.0, len(colors)))

    for label, color in zip(colors, cmap):
        color_dict[label] = color

    return color_dict


def visualize(
    tree: Tree,
    *,
    backend: Literal["matplotlib", "plotly"] = "matplotlib",
    path: str | os.PathLike | None = None,
    show: bool = False,
    edge_length_mode: Literal["attr", "uniform", "even", "rank"] = "attr",
    edge_length_attr: str = "dist",
    layout_mode: Literal["horizontal", "vertical", "circular"] = "horizontal",
    node_rank_mode: Literal["mean", "first", "last", "node"] = "mean",
    tree_style: TreeStyle | None = None,
    color_dict: dict | None = None,
    species_info: bool = False,
    scale_symbols: float = 1.0,
    fontsize: str | float | int = "medium",
    save_as: str | os.PathLike | None = None,
    **renderer_kwargs: Any,
) -> tuple[Figure, Axes] | PlotlyFigure:
    """Render a tree.

    This function calls the ``plot_tree`` function in tralda for rendering the tree according to the
    specified backend and options.
    If no ``tree_style`` is provided, the default AsymmeTree style is used:
    - Leaf colors are determined by the ``color_dict`` if provided.
    - Node symbols are determined by the event type:
        - Leaf nodes are represented by a circle with a dot inside.
        - Speciation events are represented by a circle.
        - Duplication events are represented by a square.
        - Loss events are represented by a cap.
        - Transfer events are represented by a triangle.
        - Gene conversion events are represented by a star.

    Args:
        tree: The tree to render.
        backend: ``"matplotlib"`` (default) or ``"plotly"``.
        path: If given, save the resulting figure to this path. For the matplotlib backend the
            format is inferred from the file extension by
            :meth:`~matplotlib.figure.Figure.savefig` (e.g. ``.png``, ``.pdf``, ``.svg``). For the
            plotly backend, a ``.html`` extension writes a self-contained interactive file via
            :meth:`~plotly.graph_objects.Figure.write_html`; any other extension writes a static
            image via :meth:`~plotly.graph_objects.Figure.write_image`.
        show: If ``True``, display the figure (``fig.show()``) after rendering.
        edge_length_mode: How edge lengths are determined. See tralda documentation for details.
        edge_length_attr: Node attribute read as edge length in ``"attr"`` mode.
        layout_mode: Tree orientation. See tralda documentation for details.
        node_rank_mode: How the rank of internal nodes is derived from their children. See tralda
            documentation for details.
        tree_style: Styling configuration. When ``None`` the AsymmeTree default style is used.
        color_dict: For the default style function, a dictionary containing as values the assigned
            color for each label of the (non-loss) leaves. The default is None, in which case all
            leaves are colored white.
        species_info: If True, append the reconciliation information (the 'reconc' attribute) to
            the label of each (non-loss) leaf for which it is available.
        scale_symbols: Adjust the size of the event type symbols.
        fontsize: Adjust the size of the text.
        save_as: Deprecated alias for ``path``. Use ``path`` instead.
        **renderer_kwargs: Additional keyword arguments forwarded to the renderer constructor
            (e.g. ``figsize`` / ``ax`` for matplotlib, ``width`` / ``height`` / ``hover_attrs`` for
            plotly).

    Returns:
        ``(fig, ax)`` for the matplotlib backend, or a plotly ``Figure`` for the plotly backend.
    """
    if save_as is not None:
        if path is not None:
            raise ValueError("Cannot specify both 'path' and 'save_as'. Use 'path' only.")
        path = save_as

    if species_info:
        tree = _with_species_info_labels(tree)

    if tree_style is None:
        tree_style = TreeStyle(style_fn=partial(_default_style_fn, color_map=color_dict))
    elif color_dict is not None:
        raise ValueError("Cannot specify color_dict when providing a custom tree_style.")

    # the default symbol size in tralda is 9, so we scale it by the scale_symbols factor
    base_node_style = NodeStyle(symbol_size=9.0 * scale_symbols, label_fontsize=fontsize)
    tree_style.default = tree_style.default.overlay(base_node_style)

    return plot_tree(
        tree,
        backend=backend,
        path=path,
        show=show,
        edge_length_mode=edge_length_mode,
        edge_length_attr=edge_length_attr,
        layout_mode=layout_mode,
        node_rank_mode=node_rank_mode,
        tree_style=tree_style,
        **renderer_kwargs,
    )


def _with_species_info_labels(tree: Tree) -> Tree:
    """Return a copy of the tree with reconciliation info appended to extant leaf labels.

    The original label is kept on the copied nodes as the 'color_key' attribute, so that color
    lookups (keyed by the original label) in the default style function still work.

    Args:
        tree: The tree to copy and relabel.

    Returns:
        A copy of the tree with the modified leaf labels.
    """
    new_tree, orig_to_new = tree.copy(mapping=True)

    for orig, new in orig_to_new.items():
        if orig.is_leaf() and getattr(orig, "event", None) != "L" and hasattr(orig, "reconc"):
            new.color_key = orig.label
            new.label = f"{orig.label} <{orig.reconc}>"

    return new_tree


def _default_style_fn(
    node: TreeNode,
    _mode: LayoutMode,
    color_map: dict | None = None,
) -> NodeStyle:
    """The default style function for tree nodes.

    This function determines the visual style of a tree node based on its attributes and the
    provided color map. Speciation events are represented by circles, duplications by squares,
    horizontal gene transfers by triangles, gene conversions by stars, and losses by caps.

    Args:
        node: The tree node to style.
        _mode: The layout mode (not used in this function).
        color_map: A dictionary mapping node labels to colors.

    Returns:
        A NodeStyle object with the appropriate styling for the node.
    """
    color_map = color_map or {}
    ns = NodeStyle()

    event = getattr(node, "event", None)

    if node.is_leaf() and event != "L":  # extant leaf
        ns.symbol = "circle-dot"
        ns.symbol_color = color_map.get(getattr(node, "color_key", node.label), "white")
    elif not node.parent:  # root
        ns.symbol = "circle-inner-ring"
        ns.symbol_color = "white"
        ns.edge_color = "black"
    elif event == "S":  # internal speciation
        ns.symbol = "circle"
        ns.symbol_color = "black"
        ns.symbol_edge_color = "none"
    elif event == "D":  # duplication
        ns.symbol = "square"
        ns.symbol_color = "white"
        ns.symbol_edge_color = "black"
    elif event == "H":  # horizontal gene transfer
        ns.symbol = "triangle-up"
        ns.symbol_color = "white"
        ns.symbol_edge_color = "black"
    elif event == "GC":  # gene conversion
        ns.symbol = "star"
        ns.symbol_color = "black"
        ns.symbol_edge_color = "none"
    elif event == "L":  # loss
        ns.symbol = "cap"

    # Highlight transferred edges (gene tree only)
    transferred = getattr(node, "transferred", 0)
    if transferred:
        ns.edge_color = "crimson"
        ns.edge_ls = "--"

    return ns
