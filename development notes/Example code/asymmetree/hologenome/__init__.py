"""Hologenome simulation.

The module provides helpers for host-symbiont simulations and auxiliary-tree construction. The
current implementation covers the holobiont layer and prepares the auxiliary tree for the later
gene-level extension.
"""

from __future__ import annotations

from .hologenome_simulation import (
    HologenomeSimulator as HologenomeSimulator,
)
from .hologenome_simulation import create_auxiliary_tree as create_auxiliary_tree
from .hologenome_simulation import to_nhx as to_nhx
