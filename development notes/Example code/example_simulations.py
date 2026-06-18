"""Run a small, reproducible collection of hologenome example simulations.

This script is the main driver for the prototype examples. It first simulates
host trees, then for each host tree it simulates one symbiont history and the
corresponding auxiliary tree, and finally it simulates gene trees inside that
auxiliary tree.

In this file, one holobiont scenario means one host tree together with one
simulated symbiont tree and the single auxiliary tree that merges them. That
auxiliary tree is the shared host-symbiont setting in which gene trees are
generated. The current script produces one gene-tree realization per parameter
combination, but conceptually the same auxiliary tree can be reused to produce
many gene trees.

Each auxiliary tree is also exported in decomposed form: its host component and
its symbiont component. This makes it possible to inspect the same scenario both
as one merged three-level object and as the collection of embedded
host-symbiont scenarios induced by the auxiliary tree.

The input parameters of ``run_example_simulations()`` control a Cartesian
product of scenarios:

- ``host_n_species`` gives the host leaf counts. Each entry generates one host
  tree. Repeated values are allowed and still create independent host trees.
- ``symbiont_dtl_rates`` provides candidate ``(duplication, transfer, loss)``
  triples for the symbiont-tree simulations.
- ``gene_dtl_rates`` provides candidate ``(duplication, transfer, loss)``
  triples for the gene-tree simulations inside the resulting auxiliary tree.
- ``alpha_values`` and ``beta_values`` provide candidate host-system loss-rate
  normalization factors for the gene-tree simulations.
- ``replace_probabilities`` provides candidate probabilities for replacing HGT
  events.
- ``transfer_distance_biases`` provides candidate recipient-bias modes for
  transfer events.
- ``n_simulations`` gives the number of repeated draws for every combination of
  symbiont rates, gene rates, alpha, beta, replacement probability, and
  transfer bias.
- ``seed`` initializes both Python's ``random`` module and NumPy's random
  generator so the full run is reproducible.
- ``output_dir`` optionally selects where the generated tables are written.

The total number of simulated holobiont scenarios is:

``len(host_n_species)``
times
``len(symbiont_dtl_rates) * len(gene_dtl_rates) * len(alpha_values)``
times
``len(beta_values) * len(replace_probabilities) * len(transfer_distance_biases) * n_simulations``

That is, each host tree is combined with every symbiont-rate triple, every
gene-rate triple, every alpha value, every beta value, every replacement
probability, every transfer-bias mode, and every repetition index.

For each holobiont scenario, the outputs store:

- the unpruned and pruned symbiont trees,
- the merged auxiliary tree,
- the host and symbiont components of that auxiliary tree,
- the unpruned and pruned gene trees, and
- the host-side and symbiont-side projections of those gene trees.

When run as a script, pass a config file with ``> section`` headers. The
required sections are ``host_n_species``, ``symbiont_dtl_rates``,
``gene_dtl_rates``, ``replace_probabilities``, and ``n_simulations``. Optional
sections are ``alpha_values``, ``beta_values``, ``transfer_distance_biases``,
``seed``, and ``output_dir``.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from collections.abc import Sequence
from datetime import date
from itertools import product
from pathlib import Path
import random as py_random

import numpy.random as np_random


SimulationRates = tuple[float, float, float]
OptionalRate = float | None
DEFAULT_TRANSFER_DISTANCE_BIASES = (None, "exponential")
DEFAULT_ALPHA_VALUES = (None,)
DEFAULT_BETA_VALUES = (None,)
DEFAULT_SEED = 13042026
REQUIRED_CONFIG_KEYS = {
    "host_n_species",
    "symbiont_dtl_rates",
    "gene_dtl_rates",
    "replace_probabilities",
    "n_simulations",
}


def run_example_simulations(
    host_n_species: tuple[int, ...],
    symbiont_dtl_rates: tuple[SimulationRates, ...],
    gene_dtl_rates: tuple[SimulationRates, ...],
    replace_probabilities: tuple[float, ...],
    n_simulations: int,
    alpha_values: tuple[OptionalRate, ...] = DEFAULT_ALPHA_VALUES,
    beta_values: tuple[OptionalRate, ...] = DEFAULT_BETA_VALUES,
    transfer_distance_biases: tuple[str | None, ...] = DEFAULT_TRANSFER_DISTANCE_BIASES,
    seed: int = DEFAULT_SEED,
    output_dir: str | Path | None = None,
) -> tuple[Series, DataFrame]:
    """Run the holobiont example simulations.

    Args:
        host_n_species: Numbers of host leaves for the host-tree simulations.
        symbiont_dtl_rates: Duplication, transfer, and loss rates for the symbiont simulations.
        gene_dtl_rates: Duplication, transfer, and loss rates for the gene-tree simulations.
        alpha_values: Symbiont-side normalization factors for host-system loss refresh.
        beta_values: Host-side normalization factors for host-system loss refresh.
        replace_probabilities: Probabilities for replacing-transfer events.
        transfer_distance_biases: Transfer-recipient bias modes.
        n_simulations: Number of repetitions for each parameter combination.
        seed: Random seed shared by ``random`` and ``numpy``.
        output_dir: Optional directory for writing the resulting tables as CSV files.

    Returns:
        A pandas series containing the serialized host trees and a pandas dataframe containing the
            serialized symbiont, auxiliary, and gene trees for each scenario.
    """
    _load_simulation_dependencies()
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
    total_iterations = _simulation_count(
        host_n_species,
        symbiont_dtl_rates,
        gene_dtl_rates,
        alpha_values,
        beta_values,
        replace_probabilities,
        transfer_distance_biases,
        n_simulations,
    )

    with tqdm(total=total_iterations, desc="Simulating scenarios", unit="scenario") as progress:
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

            for (
                symbiont_rates,
                gene_rates,
                alpha,
                beta,
                replace_probability,
                transfer_distance_bias,
                simulation_index,
            ) in product(
                symbiont_dtl_rates,
                gene_dtl_rates,
                alpha_values,
                beta_values,
                replace_probabilities,
                transfer_distance_biases,
                range(n_simulations),
            ):
                true_trees, pruned_trees, auxiliary_trees = simulator.simulate_symbiont_trees(
                    1,
                    dupl_rate=symbiont_rates[0],
                    hgt_rate=symbiont_rates[1],
                    loss_rate=symbiont_rates[2],
                    prohibit_extinction="per_family",
                    replace_prob=replace_probability,
                    transfer_distance_bias=transfer_distance_bias,
                )
                true_gene_trees, pruned_gene_trees = simulator.simulate_gene_trees(
                    dupl_rate=gene_rates[0],
                    hgt_rate=gene_rates[1],
                    loss_rate=gene_rates[2],
                    prohibit_extinction="per_family",
                    replace_prob=replace_probability,
                    transfer_distance_bias=transfer_distance_bias,
                    alpha=alpha,
                    beta=beta,
                )

                scenario_id = _scenario_id(
                    host_tree_id,
                    symbiont_rates,
                    gene_rates,
                    alpha,
                    beta,
                    replace_probability,
                    transfer_distance_bias,
                    simulation_index,
                )

                metadata = {
                    "scenario_id": scenario_id,
                    "host_tree_id": host_tree_id,
                    "symbiont_dupl_rate": symbiont_rates[0],
                    "symbiont_hgt_rate": symbiont_rates[1],
                    "symbiont_loss_rate": symbiont_rates[2],
                    "gene_dupl_rate": gene_rates[0],
                    "gene_hgt_rate": gene_rates[1],
                    "gene_loss_rate": gene_rates[2],
                    "alpha": alpha,
                    "beta": beta,
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
                progress.update(1)

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


def _simulation_count(
    host_n_species: Sequence[int],
    symbiont_dtl_rates: Sequence[SimulationRates],
    gene_dtl_rates: Sequence[SimulationRates],
    alpha_values: Sequence[OptionalRate],
    beta_values: Sequence[OptionalRate],
    replace_probabilities: Sequence[float],
    transfer_distance_biases: Sequence[str | None],
    n_simulations: int,
) -> int:
    """Return the number of scenario iterations in the Cartesian product."""
    return (
        len(host_n_species)
        * len(symbiont_dtl_rates)
        * len(gene_dtl_rates)
        * len(alpha_values)
        * len(beta_values)
        * len(replace_probabilities)
        * len(transfer_distance_biases)
        * n_simulations
    )


def _load_simulation_dependencies() -> None:
    """Load simulation dependencies after CLI argument parsing."""
    global DataFrame
    global HologenomeSimulator
    global Series
    global species_tree_n
    global tqdm
    global to_nhx
    global to_simple_newick

    from pandas import DataFrame as pandas_dataframe
    from pandas import Series as pandas_series
    from tqdm import tqdm as tqdm_progress

    from asymmetree.hologenome import HologenomeSimulator as hologenome_simulator
    from asymmetree.hologenome import to_nhx as serialize_nhx
    from asymmetree.hologenome import to_simple_newick as serialize_simple_newick
    from asymmetree.treeevolve import species_tree_n as simulate_species_tree_n

    DataFrame = pandas_dataframe
    HologenomeSimulator = hologenome_simulator
    Series = pandas_series
    species_tree_n = simulate_species_tree_n
    tqdm = tqdm_progress
    to_nhx = serialize_nhx
    to_simple_newick = serialize_simple_newick


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
    symbiont_rates: SimulationRates,
    gene_rates: SimulationRates,
    alpha: OptionalRate,
    beta: OptionalRate,
    replace_probability: float,
    transfer_distance_bias: str | None,
    simulation_index: int,
) -> str:
    """Create a compact scenario identifier."""
    return "_".join(
        (
            host_tree_id,
            f"sd{_scaled_rate(symbiont_rates[0])}",
            f"st{_scaled_rate(symbiont_rates[1])}",
            f"sl{_scaled_rate(symbiont_rates[2])}",
            f"gd{_scaled_rate(gene_rates[0])}",
            f"gt{_scaled_rate(gene_rates[1])}",
            f"gl{_scaled_rate(gene_rates[2])}",
            f"a{_scaled_optional_rate(alpha)}",
            f"be{_scaled_optional_rate(beta)}",
            f"r{_scaled_rate(replace_probability)}",
            f"b{transfer_distance_bias or 'none'}",
            f"i{simulation_index}",
        )
    )


def _scaled_rate(value: float) -> int:
    """Scale a rate for inclusion in the scenario identifier."""
    return int(round(100 * value))


def _scaled_optional_rate(value: OptionalRate) -> str:
    """Scale an optional rate for inclusion in the scenario identifier."""
    return "default" if value is None else str(_scaled_rate(value))


def read_simulation_config(config_path: str | Path) -> dict[str, object]:
    """Read a human-editable simulation config file.

    Sections start with ``> section_name`` and continue until the next section. Comma-separated
    values are accepted within each section, and blank lines or comment lines starting with ``#``
    are ignored.
    """
    sections = _read_config_sections(Path(config_path))
    missing_keys = sorted(REQUIRED_CONFIG_KEYS.difference(sections))
    if missing_keys:
        raise ValueError(f"missing required config sections: {', '.join(missing_keys)}")

    config: dict[str, object] = {
        "host_n_species": _parse_int_list(sections["host_n_species"]),
        "symbiont_dtl_rates": _parse_rate_triples(sections["symbiont_dtl_rates"]),
        "gene_dtl_rates": _parse_rate_triples(sections["gene_dtl_rates"]),
        "replace_probabilities": _parse_float_list(sections["replace_probabilities"]),
        "n_simulations": _parse_single_int(sections["n_simulations"], "n_simulations"),
        "alpha_values": DEFAULT_ALPHA_VALUES,
        "beta_values": DEFAULT_BETA_VALUES,
        "transfer_distance_biases": DEFAULT_TRANSFER_DISTANCE_BIASES,
        "seed": DEFAULT_SEED,
        "output_dir": _default_output_dir(),
    }

    if "alpha_values" in sections:
        config["alpha_values"] = _parse_optional_float_list(sections["alpha_values"])
    if "beta_values" in sections:
        config["beta_values"] = _parse_optional_float_list(sections["beta_values"])
    if "transfer_distance_biases" in sections:
        config["transfer_distance_biases"] = _parse_optional_string_list(
            sections["transfer_distance_biases"]
        )
    if "seed" in sections:
        config["seed"] = _parse_single_int(sections["seed"], "seed")
    if "output_dir" in sections:
        config["output_dir"] = _parse_single_string(sections["output_dir"], "output_dir")

    return config


def _read_config_sections(config_path: Path) -> dict[str, list[str]]:
    """Collect config values keyed by their ``> section`` headers."""
    sections: dict[str, list[str]] = {}
    current_section: str | None = None

    with config_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                current_section = line[1:].strip()
                if not current_section:
                    raise ValueError(f"empty section name at line {line_number}")
                sections.setdefault(current_section, [])
                continue
            if current_section is None:
                raise ValueError(f"value before first section header at line {line_number}")
            sections[current_section].append(line)

    return sections


def _parse_items(lines: list[str]) -> list[str]:
    """Split comma-separated config lines into stripped items."""
    items = []
    for line in lines:
        items.extend(item.strip() for item in line.split(",") if item.strip())

    return items


def _parse_int_list(lines: list[str]) -> tuple[int, ...]:
    """Parse a section as a tuple of integers."""
    values = tuple(int(item) for item in _parse_items(lines))
    if not values:
        raise ValueError("expected at least one integer value")

    return values


def _parse_float_list(lines: list[str]) -> tuple[float, ...]:
    """Parse a section as a tuple of floats."""
    values = tuple(float(item) for item in _parse_items(lines))
    if not values:
        raise ValueError("expected at least one float value")

    return values


def _parse_optional_float_list(lines: list[str]) -> tuple[OptionalRate, ...]:
    """Parse floats, accepting ``default``/``None`` for simulator defaults."""
    values = []
    for item in _parse_items(lines):
        values.append(None if item.lower() in ("default", "none", "null") else float(item))
    if not values:
        raise ValueError("expected at least one float value or default")

    return tuple(values)


def _parse_rate_triples(lines: list[str]) -> tuple[SimulationRates, ...]:
    """Parse a DTL-rate section as one triple per non-empty line."""
    triples = []
    for line in lines:
        values = tuple(float(item.strip()) for item in line.split(",") if item.strip())
        if len(values) != 3:
            raise ValueError(f"expected a DTL triple, got {line!r}")
        triples.append(values)
    if not triples:
        raise ValueError("expected at least one DTL-rate triple")

    return tuple(triples)


def _parse_optional_string_list(lines: list[str]) -> tuple[str | None, ...]:
    """Parse transfer-bias modes, accepting ``None``/``null`` for unbiased transfer."""
    values = []
    for item in _parse_items(lines):
        values.append(None if item.lower() in ("none", "null") else item)
    if not values:
        raise ValueError("expected at least one transfer-distance bias mode")

    return tuple(values)


def _parse_single_int(lines: list[str], section_name: str) -> int:
    """Parse a section containing exactly one integer value."""
    values = _parse_int_list(lines)
    if len(values) != 1:
        raise ValueError(f"{section_name} expects exactly one integer")

    return values[0]


def _parse_single_string(lines: list[str], section_name: str) -> str:
    """Parse a section containing exactly one string value."""
    values = _parse_items(lines)
    if len(values) != 1:
        raise ValueError(f"{section_name} expects exactly one value")

    return values[0]


def _default_output_dir() -> Path:
    """Return a dated output directory path that does not currently exist."""
    stem = f"simulations_{date.today().isoformat()}"
    index = 1
    candidate = Path(f"{stem}_{index}")
    while candidate.exists():
        index += 1
        candidate = Path(f"{stem}_{index}")

    return candidate


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
    host_trees.to_csv(output_dir / "host_trees.tsv", header=True, sep="\t")
    host_trees_simple.to_csv(output_dir / "host_trees_simple.tsv", header=True, sep="\t")
    host_trees_simple_dated.to_csv(
        output_dir / "host_trees_simple_dated.tsv",
        header=True,
        sep="\t",
    )
    symbiont_trees.to_csv(output_dir / "symbiont_scenarios.tsv", sep="\t")
    symbiont_trees_simple.to_csv(output_dir / "symbiont_scenarios_simple.tsv", sep="\t")
    symbiont_trees_simple_dated.to_csv(
        output_dir / "symbiont_scenarios_simple_dated.tsv",
        sep="\t",
    )
    species_trees.to_csv(output_dir / "species_trees.tsv", sep="\t")
    species_trees_simple.to_csv(output_dir / "species_trees_simple.tsv", sep="\t")
    species_trees_simple_dated.to_csv(output_dir / "species_trees_simple_dated.tsv", sep="\t")


def main(argv: Sequence[str] | None = None) -> int:
    """Run simulations from a config file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "config_file",
        type=Path,
        help="Path to a simulation config file using '> section' headers.",
    )
    args = parser.parse_args(argv)

    config = read_simulation_config(args.config_file)
    host_series, symbiont_dataframe = run_example_simulations(**config)
    print(host_series.head())
    print(symbiont_dataframe.head())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
