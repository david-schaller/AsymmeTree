"""Module for reading and writing sequence files.

Read and write fasta files and multiple sequence alignments in various formats.
"""

from __future__ import annotations

from pathlib import Path
from typing import TextIO

from tralda.datastructures import TreeNode


# --------------------------------------------------------------------------------------------------
#                                    Auxiliary functions
# --------------------------------------------------------------------------------------------------


def labeled_sequences(sequences: dict | list[tuple]) -> list[tuple[str, str]]:
    """Converts dictionaries to lists and keys into strings.

    If the keys are `TreeNode` instances, then their label is taken. Also works if `sequences`
    already is a list (of tuples).

    Args:
        sequences: A dictionary mapping keys to sequences, or a list of (key, sequence) tuples.

    Returns:
        A list of (label, sequence) tuples, where labels are strings.
    """
    result = []

    if isinstance(sequences, dict):
        sequences = sequences.items()

    for key, sequence in sequences:
        if isinstance(key, TreeNode):
            result.append((str(key.label), sequence))
        else:
            result.append((str(key), sequence))

    return result


# --------------------------------------------------------------------------------------------------
#                                       Fasta files
# --------------------------------------------------------------------------------------------------


def write_fasta(file_path: Path | str, sequences: dict | list[tuple]) -> None:
    """Write sequences to a fasta file.

    Args:
        file_path: The path to the output fasta file.
        sequences: A dictionary mapping keys to sequences, or a list of (key, sequence) tuples. If
            the keys are `TreeNode` instances, then their label is taken as the sequence label.
    """
    file_path = Path(file_path)
    sequences = labeled_sequences(sequences)

    with file_path.open("w") as f:
        for label, seq in sequences:
            f.write(f">{label}\n")
            pos = 0
            while pos < len(seq):
                f.write(seq[pos : min(pos + 80, len(seq))])
                pos += 80
                f.write("\n")


# --------------------------------------------------------------------------------------------------
#                               Multiple Sequence Alignments
# --------------------------------------------------------------------------------------------------


def write_alignment(
    file_path: Path | str,
    alignment: dict | list[tuple],
    alignment_format: str = "phylip",
) -> None:
    """Write a multiple sequence alignment to a file.

    Args:
        file_path: The path to the output alignment file.
        alignment: A dictionary mapping keys to sequences, or a list of (key, sequence) tuples. If
            the keys are `TreeNode` instances, then their label is taken as the sequence label.
        alignment_format: The format of the alignment file. Options are "phylip", "clustal", and
            "pretty".

    Raises:
        ValueError: If the specified alignment format is not available.
    """
    alignment = labeled_sequences(alignment)
    file_path = Path(file_path)

    with file_path.open("w") as f:
        if alignment_format == "phylip":
            _write_phylip(f, alignment)
        elif alignment_format == "clustal":
            _write_clustal(f, alignment)
        elif alignment_format == "pretty":
            _write_pretty(f, alignment)
        else:
            raise ValueError(f"alignment format '{alignment_format}' is not available")


def _check_alignment(alignment: list[tuple]) -> tuple[int, int]:
    """Check the alignment and return the maximal label length and the sequence length.

    Args:
        alignment: A list of (label, sequence) tuples.

    Returns:
        A tuple containing the maximal label length and the sequence length.

    Raises:
        ValueError: If the aligned sequences do not have the same length.
    """
    max_length = 0  # maximal label length
    seq_length = None  # length of the aligned sequences

    for label, seq in alignment:
        if len(str(label)) > max_length:
            max_length = len(str(label))

        if seq_length is None:
            seq_length = len(seq)
        elif seq_length != len(seq):
            raise ValueError("aligned sequences must have the same length")

    return max_length, seq_length


def _write_phylip(f: TextIO, alignment: list[tuple]) -> None:
    """Write the alignment in phylip format.

    Args:
        f: The file object to which the alignment should be written.
        alignment: A list of (label, sequence) tuples.
    """
    max_length, seq_length = _check_alignment(alignment)
    max_length = max(max_length, 9)

    f.write(f"  {len(alignment)} {seq_length} i")

    format_str = f"\n{{:<{max_length + 1}}}"
    current = 0

    while current < seq_length:
        end = min(seq_length, current + 50)

        for label, seq in alignment:
            if current == 0:
                f.write(format_str.format(label))
            else:
                f.write(format_str.format(""))

            f.write(" ".join(seq[i : min(i + 10, end)] for i in range(current, end, 10)))

        if end != seq_length:
            f.write("\n")

        current += 50


def _write_clustal(f: TextIO, alignment: list[tuple]) -> None:
    """Write the alignment in clustal format.

    Args:
        f: The file object to which the alignment should be written.
        alignment: A list of (label, sequence) tuples.
    """
    max_length, seq_length = _check_alignment(alignment)

    f.write("CLUSTAL W (1.8) multiple sequence alignment\n\n")

    format_str = f"\n{{:<{max_length + 4}}}{{}}"
    current = 0

    while current < seq_length:
        end = min(seq_length, current + 60)

        for label, seq in alignment:
            f.write(format_str.format(label, seq[current:end]))

        if end != seq_length:
            f.write("\n\n")

        current += 60


def _write_pretty(f: TextIO, alignment: list[tuple]) -> None:
    """Write the alignment in a pretty format.

    Args:
        f: The file object to which the alignment should be written.
        alignment: A list of (label, sequence) tuples.
    """
    _, seq_length = _check_alignment(alignment)

    f.write("  1          11         21         31         41       50\n")
    f.write("  |          |          |          |          |        |")

    current = 0

    while current < seq_length:
        end = min(seq_length, current + 50)

        for label, seq in alignment:
            count = end - current - seq[current:end].count("-")
            seq_string = " ".join(seq[i : min(i + 10, end)] for i in range(current, end, 10))
            f.write(f"\n  {seq_string:54}{count:>6} {label}")

        if end != seq_length:
            f.write("\n\n\n")

        current += 50
