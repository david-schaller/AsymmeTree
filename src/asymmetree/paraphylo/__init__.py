# -*- coding: utf-8 -*-

"""Estimation of species trees from orthology/paralogy.

The subpackage asymmetree.paraphylo contains a method to compute rooted species
tree from orthology/paralogy relations. This is a reimplementation of ParaPhylo
which uses heuristics for the NP-hard optimization steps instead of exact ILP
solutions.
"""

from asymmetree.paraphylo.SpeciesTreeFromParalogs import TreeReconstructor
from asymmetree.paraphylo.SpeciesTreeFromProteinOrtho import reconstruct_from_proteinortho