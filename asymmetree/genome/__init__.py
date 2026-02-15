# -*- coding: utf-8 -*-

"""Genome simulation.

The module asymmetree.genome.GenomeSimulation provides functions that combine
the simulation of phylogenetic trees and sequences. In particular, the class
GenomeSimulator combines multiple steps described in the previous sections in
order to conveniently simulate whole genomes/proteomes. The (optional) output
directory contains serialized trees, fasta files, and the true alignments.
"""

from asymmetree.genome.GenomeSimulation import GenomeSimulator