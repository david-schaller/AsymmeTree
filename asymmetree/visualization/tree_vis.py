"""Species and gene tree visualization."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from tralda.datastructures import Tree


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
        Two dictionaries, one for the species tree and one for the gene tree, containing as keys the
        labels of the extant species/genes and as values the assigned colors.
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

    if len(colors) <= 10:
        cmap = plt.get_cmap("tab10")(np.arange(len(colors), dtype=int))
    else:
        cmap = plt.get_cmap("jet")(np.linspace(0, 1.0, len(colors)))

    for label, color in zip(colors, cmap):
        color_dict[label] = color

    return color_dict


def visualize(
    tree: Tree,
    color_dict: dict | None = None,
    save_as: str | bool = False,
    scale_symbols: float = 1.0,
    fontsize: str | float | int = "medium",
) -> None:
    """Visualize a tree.

    Args:
        tree: A tree whose nodes have the 'dist' attribute.
        color_dict: A dictionary containing as values the assigned color for each label of the
            (non-loss) leaves. The default is None, in which case all leaves are colored white.
        save_as: Save the image to a file. The default is False, in which case the image is not
            saved.
        scale_symbols: Adjust the size of the event type symbols.
        fontsize: Adjust the size of the text.
    """
    vis = Visualizer(tree, color_dict=color_dict, scale_symbols=scale_symbols, fontsize=fontsize)
    vis.draw(save_as=save_as)


class Visualizer:
    """Species and gene tree visualization."""

    def __init__(
        self,
        tree: Tree,
        color_dict: dict | None = None,
        scale_symbols: float = 1.0,
        fontsize: str | float | int = "medium",
        species_info: bool = False,
    ) -> None:
        """Construct a visualizer for a tree.

        Args:
            tree: A tree whose nodes have the 'dist' attribute.
            color_dict: A dictionary containing as values the assigned color for each label of the
                (non-loss) leaves. The default is None, in which case all leaves are colored white.
            scale_symbols: Adjust the size of the event type symbols.
            fontsize: Adjust the size of the text.
            species_info: If True, write the reconciliation information behind the leaf label if
                available.
        """
        self.tree = tree

        if color_dict:
            self.colors = color_dict
        else:
            self.colors = {v.label: "white" for v in tree.preorder()}

        self.species_info = species_info

        self.symbolsize = 0.03 * scale_symbols
        self.fontsize = fontsize
        self.symbollw = 0.04
        self.leafs_per_vertical_unit = 15
        self.symbol_zorder = 3

        self.distance_dict = {}
        self.leaf_counter = sum(1 for _ in tree.leaves())
        self.node_positions = {}

    def draw(self, distance_mode: str = "evolutionary", save_as: str | Path | None = None) -> None:
        """Draw the tree and optionally save it to file.

        Args:
            distance_mode: Information that is used for the length of the species edges. Currently
                only 'evolutionary' is supported, in which case the 'dist' attribute is used.
            save_as: Save the image to a file. The default is None, in which case the image is not
                saved.
        """
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect("equal")
        self.ax.invert_yaxis()

        self.initial_traversal(distance_mode)
        # self.assign_colors()
        self.assign_positions()
        self.draw_edges()
        self.draw_nodes()

        self.ax.axvline(x=0, linewidth=1, color="grey", linestyle="--")
        self.ax.set_yticks([])
        self.ax.spines["top"].set_visible(False)
        self.ax.spines["left"].set_visible(False)
        self.ax.spines["right"].set_visible(False)
        self.ax.tick_params(axis="x", labelsize=self.fontsize)

        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        self.fig.set_size_inches(5 * abs(xmax - xmin), 5 * abs(ymax - ymin) + 0.4)
        plt.tight_layout()
        if save_as:
            plt.savefig(save_as, dpi=300)
        else:
            plt.show()

    def initial_traversal(self, distance_mode: str) -> None:
        """Traversal to compute the distance of each node from the root and the maximum distance.

        Args:
            distance_mode: Information that is used for the length of the species edges. Currently
                only 'evolutionary' is supported, in which case the 'dist' attribute is used.

        Raises:
            ValueError: If the given distance mode is not supported.
        """
        xmax = 0.0

        for v in self.tree.preorder():
            if distance_mode == "evolutionary":
                if not v.parent:
                    self.distance_dict[v] = 0.0
                else:
                    self.distance_dict[v] = self.distance_dict[v.parent] + v.dist
                    if self.distance_dict[v] > xmax:
                        xmax = self.distance_dict[v]
            elif distance_mode == "divergence time":
                raise ValueError("divergence time not yet implemented")
            else:
                raise ValueError(f"distance mode '{distance_mode}' not supported")

        self.ax.set_xlim(-0.1, xmax + 0.5)

    def assign_positions(self):
        """Assign x and y coordinates to each node."""
        ymax = (self.leaf_counter - 1) / self.leafs_per_vertical_unit
        self.ax.set_ylim(ymax + 0.1, -self.symbolsize * 0.6)

        yposition = 0
        for v in self.tree.postorder():
            if not v.children:
                self.node_positions[v] = (self.distance_dict[v], yposition)
                yposition += 1 / self.leafs_per_vertical_unit
            else:
                ymean = (
                    self.node_positions[v.children[0]][1] + self.node_positions[v.children[-1]][1]
                ) / 2
                self.node_positions[v] = (self.distance_dict[v], ymean)

    def draw_edges(self):
        """Draw edges between nodes."""
        for v in self.tree.preorder():
            if v.parent:
                if hasattr(v, "transferred") and v.transferred:
                    color = "red"
                else:
                    color = "black"
                self.ax.plot(
                    [self.node_positions[v.parent][0], self.node_positions[v][0]],
                    [self.node_positions[v][1], self.node_positions[v][1]],
                    color=color,
                    linestyle="-",
                    linewidth=1,
                )
            if v.children:
                self.ax.plot(
                    [self.node_positions[v][0], self.node_positions[v][0]],
                    [self.node_positions[v.children[0]][1], self.node_positions[v.children[-1]][1]],
                    color="black",
                    linestyle="-",
                    linewidth=1,
                )

    def draw_nodes(self):
        """Draw nodes according to their event type."""
        for v in self.tree.preorder():
            if not v.parent and not v.event:
                self.draw_root(*self.node_positions[v])
            elif v.children:
                if v.event == "D":
                    self.draw_dupl(*self.node_positions[v])
                elif v.event == "S":
                    self.draw_spec(*self.node_positions[v])
                elif v.event == "H":
                    self.draw_hgt(*self.node_positions[v])
                elif v.event == "GC":
                    self.draw_gc(*self.node_positions[v])
            else:
                if v.event == "L":
                    self.draw_loss(*self.node_positions[v])
                else:
                    x, y = self.node_positions[v]
                    self.draw_leaf(x, y, color=self.colors[v.label])
                    if not self.species_info:
                        text = str(v.label)
                    elif hasattr(v, "reconc"):
                        text = "{} <{}>".format(v.label, v.reconc)
                    self.write_label(x + self.symbolsize + 0.02, y, text)

    def draw_leaf(self, x: float, y: float, color: str = "white", leftalign: bool = False) -> None:
        """Draw a leaf node.

        Args:
            x: The x-coordinate of the leaf node.
            y: The y-coordinate of the leaf node.
            color: The color to fill the leaf node.
            leftalign: If True, align the symbol to the left of the given coordinates. The default
                is False, in which case the symbol is centered on the given coordinates.
        """
        if leftalign:
            x += self.symbolsize / 2
        fill = mpatches.Circle(
            (x, y), self.symbolsize / 2, color=color, fill=True, zorder=self.symbol_zorder
        )
        self.ax.add_patch(fill)
        outer = mpatches.Circle(
            (x, y),
            self.symbolsize / 2,
            color="black",
            fill=False,
            lw=self.symbolsize / self.symbollw,
            zorder=self.symbol_zorder,
        )
        self.ax.add_patch(outer)
        little = mpatches.Circle(
            (x, y), self.symbolsize / 8, color="black", fill=True, zorder=self.symbol_zorder
        )
        self.ax.add_patch(little)

    def draw_loss(self, x, y):
        self.ax.plot(
            [x, x],
            [y - self.symbolsize / 2, y + self.symbolsize / 2],
            color="black",
            linestyle="-",
            linewidth=1,
        )

    def draw_root(self, x: float, y: float, rightalign: bool = False) -> None:
        """Draw a root node.

        Args:
            x: The x-coordinate of the root node.
            y: The y-coordinate of the root node.
            rightalign: If True, align the symbol to the right of the given coordinates. The default
                is False, in which case the symbol is centered on the given coordinates.
        """
        if rightalign:
            x -= self.symbolsize / 2
        fill = mpatches.Circle(
            (x, y), self.symbolsize / 2, color="white", fill=True, zorder=self.symbol_zorder
        )
        self.ax.add_patch(fill)
        outer = mpatches.Circle(
            (x, y),
            self.symbolsize / 2,
            color="black",
            fill=False,
            lw=self.symbolsize / self.symbollw,
            zorder=self.symbol_zorder,
        )
        self.ax.add_patch(outer)
        little = mpatches.Circle(
            (x, y),
            self.symbolsize / 5,
            color="black",
            fill=False,
            lw=self.symbolsize / self.symbollw,
            zorder=self.symbol_zorder,
        )
        self.ax.add_patch(little)

    def draw_spec(self, x: float, y: float) -> None:
        """Draw a speciation node.

        Args:
            x: The x-coordinate of the speciation node.
            y: The y-coordinate of the speciation node.
        """
        fill = mpatches.Circle(
            (x, y), self.symbolsize / 2, color="black", fill=True, zorder=self.symbol_zorder
        )
        self.ax.add_patch(fill)

    def draw_dupl(self, x: float, y: float) -> None:
        """Draw a duplication node.

        Args:
            x: The x-coordinate of the duplication node.
            y: The y-coordinate of the duplication node.
        """
        square = mpatches.Rectangle(
            (x - self.symbolsize / 2, y - self.symbolsize / 2),
            width=self.symbolsize,
            height=self.symbolsize,
            color="white",
            fill=True,
            zorder=self.symbol_zorder,
        )
        self.ax.add_patch(square)
        border = mpatches.Rectangle(
            (x - self.symbolsize / 2, y - self.symbolsize / 2),
            width=self.symbolsize,
            height=self.symbolsize,
            color="black",
            fill=False,
            lw=self.symbolsize / self.symbollw,
            zorder=self.symbol_zorder,
        )
        self.ax.add_patch(border)

    def draw_hgt(self, x: float, y: float) -> None:
        """Draw a horizontal gene transfer (HGT) node.

        Args:
            x: The x-coordinate of the HGT node.
            y: The y-coordinate of the HGT node.
        """
        coord = np.asarray(
            [
                [x - self.symbolsize / 2, y],
                [x + self.symbolsize / 2, y - self.symbolsize / 2],
                [x + self.symbolsize / 2, y + self.symbolsize / 2],
            ]
        )

        inner = mpatches.Polygon(
            coord, closed=True, color="white", fill=True, zorder=self.symbol_zorder
        )
        self.ax.add_patch(inner)
        outer = mpatches.Polygon(
            coord,
            closed=True,
            color="black",
            fill=False,
            lw=self.symbolsize / self.symbollw,
            zorder=self.symbol_zorder,
        )
        self.ax.add_patch(outer)

    def draw_gc(self, x: float, y: float) -> None:
        """Draw a gene conversion (GC) node.

        Args:
            x: The x-coordinate of the GC node.
            y: The y-coordinate of the GC node.
        """
        coord = []
        R = self.symbolsize / 2
        for i in range(5):
            coord.append(
                [
                    x + 0.5 * R * np.cos(np.radians(0 + i * 72)),
                    y - 0.5 * R * np.sin(np.radians(0 + i * 72)),
                ]
            )
            coord.append(
                [x + R * np.cos(np.radians(36 + i * 72)), y - R * np.sin(np.radians(36 + i * 72))]
            )

        inner = mpatches.Polygon(
            coord, closed=True, color="black", fill=True, zorder=self.symbol_zorder
        )
        self.ax.add_patch(inner)

    def write_label(self, x: float, y: float, text: str) -> None:
        """Write a label at the specified coordinates.

        Args:
            x: The x-coordinate of the label.
            y: The y-coordinate of the label.
            text: The text of the label.
        """
        self.ax.text(
            x,
            y,
            text,
            fontsize=self.fontsize,
            horizontalalignment="left",
            verticalalignment="center",
        )
