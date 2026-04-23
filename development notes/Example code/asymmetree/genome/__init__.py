"""Genome simulation.

The module provides functions that combine the simulation of phylogenetic trees and sequences. In
particular, the class GenomeSimulator combines multiple steps described in the previous sections in
order to conveniently simulate whole genomes/proteomes. The (optional) output directory contains
serialized trees, fasta files, and the true alignments.
"""

from __future__ import annotations

from asymmetree.genome.genome_simulation import GenomeSimulator as GenomeSimulator
