"""Example simulations for the prototype hologenome module."""

from __future__ import annotations

from collections.abc import Callable
from itertools import product
from pathlib import Path
import random as py_random

from pandas import DataFrame
from pandas import Series
import numpy.random as np_random

from asymmetree.hologenome import HologenomeSimulator
from asymmetree.hologenome import to_nhx
from asymmetree.hologenome import to_simple_newick
from asymmetree.treeevolve import species_tree_n


SimulationRates = tuple[float, float, float]


def run_example_simulations(
    host_n_species: tuple[int, ...] = (5, 10, 30),
    symbiont_dtl_rates: tuple[SimulationRates, ...] = (
        (0.133, 0.266, 0.266),
        (0.3, 0.6, 0.6),
        (0.6, 1.2, 1.2),
    ),
    replace_probabilities: tuple[float, ...] = (0.0, 0.5, 1.0),
    transfer_distance_biases: tuple[str | None, ...] = (None, "exponential"),
    n_simulations: int = 5,
    seed: int = 13042026,
    output_dir: str | Path | None = None,
) -> tuple[Series, DataFrame]:
    """Run the holobiont example simulations.

    Args:
        host_n_species: Numbers of host leaves for the host-tree simulations.
        symbiont_dtl_rates: Duplication, transfer, and loss rates for the symbiont simulations.
        replace_probabilities: Probabilities for replacing-transfer events.
        transfer_distance_biases: Transfer-recipient bias modes.
        n_simulations: Number of repetitions for each parameter combination.
        seed: Random seed shared by ``random`` and ``numpy``.
        output_dir: Optional directory for writing the resulting tables as CSV files.

    Returns:
        A pandas series containing the serialized host trees and a pandas dataframe containing the
            serialized symbiont, auxiliary, and gene trees for each scenario.
    """
    py_random.seed(seed)
    np_random.seed(seed)

    host_counts: dict[int, int] = {}
    host_records: dict[str, str] = {}
    host_records_simple: dict[str, str] = {}
    host_records_simple_dated: dict[str, str] = {}
    symbiont_rows: list[dict[str, object]] = []
    symbiont_rows_simple: list[dict[str, object]] = []
    symbiont_rows_simple_dated: list[dict[str, object]] = []
    species_rows: list[dict[str, object]] = []
    species_rows_simple: list[dict[str, object]] = []
    species_rows_simple_dated: list[dict[str, object]] = []

    for n_host in host_n_species:
        host_counts[n_host] = host_counts.get(n_host, 0) + 1
        host_tree_id = f"N{n_host}.{host_counts[n_host]}"

        host_tree = species_tree_n(n_host, innovation=True)
        host_records[host_tree_id] = to_nhx(host_tree)
        host_records_simple[host_tree_id] = to_simple_newick(host_tree)
        host_records_simple_dated[host_tree_id] = to_simple_newick(
            host_tree,
            include_distances=True,
        )
        _append_species_rows(
            species_rows,
            species_rows_simple,
            species_rows_simple_dated,
            tree_id=host_tree_id,
            tree_role="host",
            tree=host_tree,
            host_tree_id=host_tree_id,
        )

        simulator = HologenomeSimulator(host_tree)

        for rates, replace_probability, transfer_distance_bias, simulation_index in product(
            symbiont_dtl_rates,
            replace_probabilities,
            transfer_distance_biases,
            range(n_simulations),
        ):
            true_trees, pruned_trees, auxiliary_trees = simulator.simulate_symbiont_trees(
                1,
                dupl_rate=rates[0],
                hgt_rate=rates[1],
                loss_rate=rates[2],
                prohibit_extinction="per_family",
                replace_prob=replace_probability,
                transfer_distance_bias=transfer_distance_bias,
            )
            true_gene_trees, pruned_gene_trees = simulator.simulate_gene_trees(
                dupl_rate=rates[0],
                hgt_rate=rates[1],
                loss_rate=rates[2],
                prohibit_extinction="per_family",
                replace_prob=replace_probability,
                transfer_distance_bias=transfer_distance_bias,
            )

            scenario_id = _scenario_id(
                host_tree_id,
                rates,
                replace_probability,
                transfer_distance_bias,
                simulation_index,
            )

            metadata = {
                "scenario_id": scenario_id,
                "host_tree_id": host_tree_id,
                "dupl_rate": rates[0],
                "hgt_rate": rates[1],
                "loss_rate": rates[2],
                "replace_prob": replace_probability,
                "transfer_distance_bias": transfer_distance_bias,
                "simulation_index": simulation_index,
            }
            scenario_trees = {
                "T_symbiont_unpruned": true_trees[0],
                "T_symbiont_pruned": pruned_trees[0],
                "T_auxiliary": auxiliary_trees[0],
                "T_gene_unpruned": true_gene_trees[0],
                "T_gene_pruned": pruned_gene_trees[0],
                "T_auxiliary_host": simulator.host_component_trees[0],
                "T_auxiliary_symbiont": simulator.symbiont_component_trees[0],
                "T_gene_host_unpruned": simulator.true_host_gene_trees[0],
                "T_gene_symbiont_unpruned": simulator.true_symbiont_gene_trees[0],
                "T_gene_host_pruned": simulator.pruned_host_gene_trees[0],
                "T_gene_symbiont_pruned": simulator.pruned_symbiont_gene_trees[0],
            }

            symbiont_rows.append(_serialize_scenario_row(metadata, scenario_trees, to_nhx))
            symbiont_rows_simple.append(
                _serialize_scenario_row(metadata, scenario_trees, to_simple_newick)
            )
            symbiont_rows_simple_dated.append(
                _serialize_scenario_row(
                    metadata,
                    scenario_trees,
                    _to_simple_dated_newick,
                )
            )
            _append_species_rows(
                species_rows,
                species_rows_simple,
                species_rows_simple_dated,
                tree_id=f"{scenario_id}.symbiont_unpruned",
                tree_role="symbiont_unpruned",
                tree=true_trees[0],
                host_tree_id=host_tree_id,
                scenario_id=scenario_id,
            )
            _append_species_rows(
                species_rows,
                species_rows_simple,
                species_rows_simple_dated,
                tree_id=f"{scenario_id}.symbiont_pruned",
                tree_role="symbiont_pruned",
                tree=pruned_trees[0],
                host_tree_id=host_tree_id,
                scenario_id=scenario_id,
            )
            _append_species_rows(
                species_rows,
                species_rows_simple,
                species_rows_simple_dated,
                tree_id=f"{scenario_id}.auxiliary_host",
                tree_role="auxiliary_host",
                tree=simulator.host_component_trees[0],
                host_tree_id=host_tree_id,
                scenario_id=scenario_id,
            )
            _append_species_rows(
                species_rows,
                species_rows_simple,
                species_rows_simple_dated,
                tree_id=f"{scenario_id}.auxiliary_symbiont",
                tree_role="auxiliary_symbiont",
                tree=simulator.symbiont_component_trees[0],
                host_tree_id=host_tree_id,
                scenario_id=scenario_id,
            )

    host_trees = Series(host_records, name="T_host")
    host_trees_simple = Series(host_records_simple, name="T_host")
    host_trees_simple_dated = Series(host_records_simple_dated, name="T_host")
    symbiont_trees = DataFrame(symbiont_rows).set_index("scenario_id")
    symbiont_trees_simple = DataFrame(symbiont_rows_simple).set_index("scenario_id")
    symbiont_trees_simple_dated = DataFrame(symbiont_rows_simple_dated).set_index(
        "scenario_id"
    )
    species_trees = DataFrame(species_rows).set_index("tree_id")
    species_trees_simple = DataFrame(species_rows_simple).set_index("tree_id")
    species_trees_simple_dated = DataFrame(species_rows_simple_dated).set_index("tree_id")

    if output_dir is not None:
        _write_outputs(
            Path(output_dir),
            host_trees,
            host_trees_simple,
            host_trees_simple_dated,
            symbiont_trees,
            symbiont_trees_simple,
            symbiont_trees_simple_dated,
            species_trees,
            species_trees_simple,
            species_trees_simple_dated,
        )

    return host_trees, symbiont_trees


def _serialize_scenario_row(
    metadata: dict[str, object],
    trees: dict[str, object],
    serializer: Callable[..., str],
) -> dict[str, object]:
    """Serialize all tree-valued columns for one scenario."""
    row = dict(metadata)

    for column, tree in trees.items():
        row[column] = serializer(tree)

    return row


def _append_species_rows(
    rows: list[dict[str, object]],
    simple_rows: list[dict[str, object]],
    simple_dated_rows: list[dict[str, object]],
    *,
    tree_id: str,
    tree_role: str,
    tree: object,
    host_tree_id: str,
    scenario_id: str | None = None,
) -> None:
    """Append one species-tree row in the three supported output formats."""
    metadata = {
        "tree_id": tree_id,
        "tree_role": tree_role,
        "host_tree_id": host_tree_id,
        "scenario_id": scenario_id,
    }

    rows.append(_species_tree_row(metadata, to_nhx(tree)))
    simple_rows.append(_species_tree_row(metadata, to_simple_newick(tree)))
    simple_dated_rows.append(_species_tree_row(metadata, _to_simple_dated_newick(tree)))


def _species_tree_row(metadata: dict[str, object], tree: str) -> dict[str, object]:
    """Create one row of a species-tree table."""
    row = dict(metadata)
    row["tree"] = tree

    return row


def _to_simple_dated_newick(tree: object) -> str:
    """Serialize a tree in the simplified dated format."""
    return to_simple_newick(tree, include_distances=True)


def _scenario_id(
    host_tree_id: str,
    rates: SimulationRates,
    replace_probability: float,
    transfer_distance_bias: str | None,
    simulation_index: int,
) -> str:
    """Create a compact scenario identifier."""
    return "_".join(
        (
            host_tree_id,
            f"d{_scaled_rate(rates[0])}",
            f"t{_scaled_rate(rates[1])}",
            f"l{_scaled_rate(rates[2])}",
            f"r{_scaled_rate(replace_probability)}",
            f"b{transfer_distance_bias or 'none'}",
            f"i{simulation_index}",
        )
    )


def _scaled_rate(value: float) -> int:
    """Scale a rate for inclusion in the scenario identifier."""
    return int(round(100 * value))


def _write_outputs(
    output_dir: Path,
    host_trees: Series,
    host_trees_simple: Series,
    host_trees_simple_dated: Series,
    symbiont_trees: DataFrame,
    symbiont_trees_simple: DataFrame,
    symbiont_trees_simple_dated: DataFrame,
    species_trees: DataFrame,
    species_trees_simple: DataFrame,
    species_trees_simple_dated: DataFrame,
) -> None:
    """Write the example outputs to CSV files in all supported tree formats."""
    output_dir.mkdir(parents=True, exist_ok=True)
    host_trees.to_csv(output_dir / "host_trees.csv", header=True)
    host_trees_simple.to_csv(output_dir / "host_trees_simple.csv", header=True)
    host_trees_simple_dated.to_csv(
        output_dir / "host_trees_simple_dated.csv",
        header=True,
    )
    symbiont_trees.to_csv(output_dir / "symbiont_scenarios.csv")
    symbiont_trees_simple.to_csv(output_dir / "symbiont_scenarios_simple.csv")
    symbiont_trees_simple_dated.to_csv(output_dir / "symbiont_scenarios_simple_dated.csv")
    species_trees.to_csv(output_dir / "species_trees.csv")
    species_trees_simple.to_csv(output_dir / "species_trees_simple.csv")
    species_trees_simple_dated.to_csv(output_dir / "species_trees_simple_dated.csv")


if __name__ == "__main__":
    host_series, symbiont_dataframe = run_example_simulations()
    print(host_series.head())
    print(symbiont_dataframe.head())
