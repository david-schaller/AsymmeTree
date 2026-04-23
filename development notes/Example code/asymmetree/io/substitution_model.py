"""Module for reading and writing substitution models."""

from __future__ import annotations

from pathlib import Path


def parse_paml(
    file_path: Path | str,
    model_type: str = "a",
) -> tuple[list[list[float]], list[float]]:
    """Parse a PAML substitution model file.

    Args:
        file_path: The path to the PAML file containing the substitution model.
        model_type: The type of the substitution model, either "a" for amino acid or "n" for
            nucleotide.

    Returns:
        A tuple containing the exchangeability matrix and the equilibrium frequencies.

    Raises:
        ValueError: If the model type is not supported.
        ValueError: If the number of equilibrium frequencies is not correct.
    """
    if model_type == "a":
        # amino acid model
        alphabet_size = 20
    elif model_type == "n":
        # nucleotide model
        alphabet_size = 4
    else:
        raise ValueError(f"model type '{model_type}' not supported")

    file_path = Path(file_path)

    with file_path.open("r") as f:
        exchangeability_matrix = [[0.0 for j in range(alphabet_size)] for i in range(alphabet_size)]

        line = f.readline().strip()
        while line == "":
            line = f.readline().strip()

        for i in range(1, alphabet_size):
            line = [float(item) for item in line.split()]

            for j in range(len(line)):
                exchangeability_matrix[i][j] = line[j]
                exchangeability_matrix[j][i] = line[j]

            line = f.readline().strip()

        line = f.readline().strip()
        while line == "":
            line = f.readline().strip()

        stat_freqs = [float(item) for item in line.split()]

        if len(stat_freqs) != alphabet_size:
            raise ValueError("wrong no. of equilibrium frequencies, check paml file")

        return exchangeability_matrix, stat_freqs
