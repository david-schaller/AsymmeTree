"""Example simulations for the prototype hologenome module."""

from __future__ import annotations

from itertools import product
from pathlib import Path
import random as py_random

from pandas import DataFrame
from pandas import Series
import numpy.random as np_random

from asymmetree.hologenome import HologenomeSimulator
from asymmetree.hologenome import to_nhx
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
    symbiont_rows: list[dict[str, object]] = []

    for n_host in host_n_species:
        host_counts[n_host] = host_counts.get(n_host, 0) + 1
        host_tree_id = f"N{n_host}.{host_counts[n_host]}"

        host_tree = species_tree_n(n_host, innovation=True)
        host_records[host_tree_id] = to_nhx(host_tree)

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

            symbiont_rows.append(
                {
                    "scenario_id": scenario_id,
                    "host_tree_id": host_tree_id,
                    "dupl_rate": rates[0],
                    "hgt_rate": rates[1],
                    "loss_rate": rates[2],
                    "replace_prob": replace_probability,
                    "transfer_distance_bias": transfer_distance_bias,
                    "simulation_index": simulation_index,
                    "T_symbiont_unpruned": to_nhx(true_trees[0]),
                    "T_symbiont_pruned": to_nhx(pruned_trees[0]),
                    "T_auxiliary": to_nhx(auxiliary_trees[0]),
                    "T_gene_unpruned": to_nhx(true_gene_trees[0]),
                    "T_gene_pruned": to_nhx(pruned_gene_trees[0]),
                    "T_auxiliary_host": to_nhx(simulator.host_component_trees[0]),
                    "T_auxiliary_symbiont": to_nhx(simulator.symbiont_component_trees[0]),
                    "T_gene_host_unpruned": to_nhx(simulator.true_host_gene_trees[0]),
                    "T_gene_symbiont_unpruned": to_nhx(
                        simulator.true_symbiont_gene_trees[0]
                    ),
                    "T_gene_host_pruned": to_nhx(simulator.pruned_host_gene_trees[0]),
                    "T_gene_symbiont_pruned": to_nhx(
                        simulator.pruned_symbiont_gene_trees[0]
                    ),
                }
            )

    host_trees = Series(host_records, name="T_host")
    symbiont_trees = DataFrame(symbiont_rows).set_index("scenario_id")

    if output_dir is not None:
        _write_outputs(Path(output_dir), host_trees, symbiont_trees)

    return host_trees, symbiont_trees


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


def _write_outputs(output_dir: Path, host_trees: Series, symbiont_trees: DataFrame) -> None:
    """Write the example outputs to CSV files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    host_trees.to_csv(output_dir / "host_trees.csv", header=True)
    symbiont_trees.to_csv(output_dir / "symbiont_scenarios.csv")


if __name__ == "__main__":
    host_series, symbiont_dataframe = run_example_simulations()
    print(host_series.head())
    print(symbiont_dataframe.head())
