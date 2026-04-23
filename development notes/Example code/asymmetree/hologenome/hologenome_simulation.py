"""Module for the simulation of hologenomes."""

from __future__ import annotations

from pathlib import Path

from tralda.datastructures import Tree
from tralda.datastructures import TreeNode

from asymmetree.holoevolve import dated_gene_tree as dated_gene_tree_auxiliary
from asymmetree.holoevolve import prune_losses as prune_auxiliary_losses
from asymmetree.treeevolve import dated_gene_tree as dated_gene_tree_holobiont
from asymmetree.treeevolve import prune_losses as prune_holobiont_losses
from asymmetree.utils.phylogenetic_trees import distance_from_timing


def create_auxiliary_tree(
    host_tree: Tree,
    symbiont_tree: Tree,
    root_offset: float = 1e-6,
) -> Tree:
    """Create the auxiliary tree for a host-symbiont scenario.

    The auxiliary tree is constructed from copies of the host and symbiont trees. The copied host
    labels are prefixed with ``H`` and the copied symbiont labels with ``S`` to keep the merged
    tree collision-safe. Their planted roots are joined below a new auxiliary root ``AR`` and a new
    planted root ``A0`` is added above it. For later three-level simulations, the attribute
    ``reconc`` is used as the host map ``mu_A``.

    Args:
        host_tree: The host tree.
        symbiont_tree: The symbiont tree reconciled with the host tree.
        root_offset: Positive offset used for the new planted root of the auxiliary tree.

    Returns:
        The auxiliary tree for the given host-symbiont pair.

    Raises:
        TypeError: If one of the input trees is not a non-empty tree.
        ValueError: If ``root_offset`` is not positive.
    """
    _validate_tree(host_tree, "host tree")
    _validate_tree(symbiont_tree, "symbiont tree")

    if root_offset <= 0.0:
        raise ValueError("root_offset must be > 0")

    host_copy = host_tree.copy()
    symbiont_copy = symbiont_tree.copy()

    host_label_map = _annotate_host_tree(host_copy)
    _annotate_symbiont_tree(symbiont_copy, host_label_map)

    merged_tstamp = max(host_copy.root.tstamp, symbiont_copy.root.tstamp)

    merged_root = TreeNode(
        label="AR",
        event="S",
        reconc="AR",
        tstamp=merged_tstamp,
        level="auxiliary",
        source_label="rho_A",
    )
    merged_root.add_child(host_copy.root)
    merged_root.add_child(symbiont_copy.root)

    planted_root = TreeNode(
        label="A0",
        event=None,
        reconc="A0",
        tstamp=merged_tstamp + root_offset,
        level="auxiliary",
        source_label="0_A",
    )
    planted_root.add_child(merged_root)

    auxiliary_tree = Tree(planted_root)
    distance_from_timing(auxiliary_tree)

    return auxiliary_tree


def to_nhx(
    tree: Tree,
    annotation_attributes: tuple[str, ...] = ("event", "reconc", "tstamp", "transferred", "level"),
    include_inner_labels: bool = True,
) -> str:
    """Serialize a tree in general NHX format.

    Args:
        tree: The tree to serialize.
        annotation_attributes: Node attributes that shall be exported in the NHX annotation block.
        include_inner_labels: If True, inner-node labels are included in the output.

    Returns:
        The serialized tree in NHX format.
    """
    _validate_tree(tree, "tree")

    def _format_value(value: object) -> str:
        if isinstance(value, (tuple, list)):
            return ",".join(str(item) for item in value)
        return str(value).replace(" ", "_")

    def _annotation(node: TreeNode) -> str:
        entries = []
        for attr in annotation_attributes:
            value = getattr(node, attr, None)
            if value in (None, ""):
                continue
            entries.append(f"{attr.upper()}={_format_value(value)}")
        return "[&&NHX:{}]".format(":".join(entries)) if entries else ""

    def _serialize(node: TreeNode) -> str:
        if node.children:
            children = ",".join(_serialize(child) for child in node.children)
            label = str(getattr(node, "label", "")) if include_inner_labels else ""
            token = f"({children}){label}"
        else:
            token = str(getattr(node, "label", ""))

        if hasattr(node, "dist"):
            token += f":{node.dist}"

        token += _annotation(node)

        return token

    return _serialize(tree.root) + ";"


class HologenomeSimulator:
    """Class for the simulation of host-symbiont scenarios."""

    def __init__(self, host_tree: Tree, outdir: str | Path | None = None) -> None:
        """Constructor for the HologenomeSimulator class.

        Args:
            host_tree: The host tree along which symbiont trees are simulated.
            outdir: The path to a directory into which the results of the simulation are written.

        Raises:
            TypeError: If ``host_tree`` is not a non-empty tree.
        """
        _validate_tree(host_tree, "host tree")

        self.host_tree = host_tree
        self.outdir = Path(outdir) if outdir is not None else None

        if self.outdir:
            self._check_outdir()
            self.host_tree.serialize(self.outdir / "host_tree.json", mode="json")

        self.true_symbiont_trees: list[Tree] = []
        self.pruned_symbiont_trees: list[Tree] = []
        self.auxiliary_trees: list[Tree] = []
        self.true_gene_trees: list[Tree] = []
        self.pruned_gene_trees: list[Tree] = []

    def simulate_symbiont_trees(
        self,
        n: int,
        **kwargs: object,
    ) -> tuple[list[Tree], list[Tree], list[Tree]]:
        """Simulate symbiont trees and create one auxiliary tree per host-symbiont pair.

        The simulation currently covers the holobiont layer only. The arguments in ``kwargs`` are
        forwarded to :func:`treeevolve.dated_gene_tree`.

        Args:
            n: Number of symbiont trees to simulate.
            **kwargs: Parameters forwarded to :func:`treeevolve.dated_gene_tree`.

        Returns:
            The simulated unpruned symbiont trees, their pruned versions, and the auxiliary trees.
            The auxiliary trees are constructed from the unpruned symbiont trees so that later
            extensions can still access extinct symbiont lineages.

        Raises:
            ValueError: If ``n`` is not an integer greater than or equal to zero.
        """
        if not isinstance(n, int) or n < 0:
            raise ValueError("n must be an int >= 0")

        self.true_symbiont_trees = [
            dated_gene_tree_holobiont(self.host_tree, **kwargs) for _ in range(n)
        ]
        self.pruned_symbiont_trees = [
            prune_holobiont_losses(tree) for tree in self.true_symbiont_trees
        ]
        self.auxiliary_trees = [
            create_auxiliary_tree(self.host_tree, tree) for tree in self.true_symbiont_trees
        ]
        self.true_gene_trees = []
        self.pruned_gene_trees = []

        if self.outdir:
            for i, tree in enumerate(self.true_symbiont_trees):
                tree.serialize(self.outdir / "true_symbiont_trees" / f"symbiont_tree{i}.json")

            for i, tree in enumerate(self.pruned_symbiont_trees):
                tree.serialize(
                    self.outdir / "pruned_symbiont_trees" / f"symbiont_tree{i}.json"
                )

            for i, tree in enumerate(self.auxiliary_trees):
                tree.serialize(self.outdir / "auxiliary_trees" / f"auxiliary_tree{i}.json")

        return self.true_symbiont_trees, self.pruned_symbiont_trees, self.auxiliary_trees

    def simulate_gene_trees(self, **kwargs: object) -> tuple[list[Tree], list[Tree]]:
        """Simulate one gene tree inside each previously generated auxiliary tree.

        The arguments in ``kwargs`` are forwarded to :func:`holoevolve.dated_gene_tree`.

        Args:
            **kwargs: Parameters forwarded to :func:`holoevolve.dated_gene_tree`.

        Returns:
            The simulated unpruned gene trees and their pruned versions.

        Raises:
            RuntimeError: If no auxiliary trees are available yet.
        """
        if not self.auxiliary_trees:
            raise RuntimeError(
                "simulate_symbiont_trees() must be called before simulate_gene_trees()"
            )

        self.true_gene_trees = [
            dated_gene_tree_auxiliary(tree, **kwargs) for tree in self.auxiliary_trees
        ]
        self.pruned_gene_trees = [prune_auxiliary_losses(tree) for tree in self.true_gene_trees]

        if self.outdir:
            for i, tree in enumerate(self.true_gene_trees):
                tree.serialize(self.outdir / "true_gene_trees" / f"gene_tree{i}.json")

            for i, tree in enumerate(self.pruned_gene_trees):
                tree.serialize(self.outdir / "pruned_gene_trees" / f"gene_tree{i}.json")

        return self.true_gene_trees, self.pruned_gene_trees

    def _check_outdir(self) -> None:
        """Check the output directory and create the required subdirectories.

        Raises:
            FileExistsError: If the specified output path exists but is not a directory.
        """
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True, exist_ok=True)
        elif not self.outdir.is_dir():
            raise FileExistsError(f"'{self.outdir}' is not a directory")

        for directory in (
            "true_symbiont_trees",
            "pruned_symbiont_trees",
            "auxiliary_trees",
            "true_gene_trees",
            "pruned_gene_trees",
        ):
            path = self.outdir / directory
            if not path.exists():
                path.mkdir(parents=True, exist_ok=True)


def _validate_tree(tree: Tree, tree_name: str) -> None:
    """Validate that the object is a non-empty tree."""
    if not isinstance(tree, Tree) or not tree.root:
        raise TypeError(f"{tree_name} must be a non-empty tree of type 'Tree'")


def _annotate_host_tree(tree: Tree) -> dict[object, str]:
    """Prefix host labels and initialize the host component of ``mu_A``."""
    label_map = {}

    for node in tree.preorder():
        source_label = getattr(node, "label", None)
        prefixed_label = _prefix_label("H", source_label)

        node.source_label = source_label
        node.level = "host"
        node.label = prefixed_label
        node.reconc = prefixed_label

        label_map[source_label] = prefixed_label

    return label_map


def _annotate_symbiont_tree(tree: Tree, host_label_map: dict[object, str]) -> None:
    """Prefix symbiont labels and re-map their host reconciliation."""
    for node in tree.preorder():
        source_label = getattr(node, "label", None)
        source_reconc = getattr(node, "reconc", None)

        node.source_label = source_label
        node.source_reconc = source_reconc
        node.level = "symbiont"
        node.label = _prefix_label("S", source_label)
        node.reconc = _map_host_reconciliation(source_reconc, host_label_map)


def _map_host_reconciliation(
    reconciliation: object,
    host_label_map: dict[object, str],
) -> object:
    """Map a symbiont reconciliation entry to the prefixed host labels."""
    if reconciliation is None:
        return None

    if isinstance(reconciliation, (tuple, list)):
        return tuple(_map_host_reconciliation(item, host_label_map) for item in reconciliation)

    return host_label_map.get(reconciliation, _prefix_label("H", reconciliation))


def _prefix_label(prefix: str, label: object) -> str:
    """Create a prefixed, collision-safe label."""
    if label in (None, ""):
        return f"{prefix}0"
    return f"{prefix}{label}"
