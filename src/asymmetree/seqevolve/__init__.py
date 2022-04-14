# -*- coding: utf-8 -*-

"""Evolving sequences along trees.

The subpackage asymmetree.seqevolve contains modules for the simulation of
nucleotide or amino acid sequences along a phylogenetic tree using
time-continuous Markov chains, as usually applied for this purpose (for
textbooks, see e.g. [1], [2], and [3]).

These models typically take a substitution-rate matrix and the equilibrium
frequencies of the states (i.e. the nucleotides or amino acids as input).
Moreover, insertions and deletions (indels) and heterogeneity among the sites 
can be simulated.

References
----------
.. [1] J. Felsenstein.
   Inferring Phylogenies. Sinauer Associates, Sunderland, Mass, 2004.
   ISBN 978-0-87893-177-4.
.. [2] Z. Yang.
   Computational molecular evolution.
   Oxford series in ecology and evolution. Oxford University Press, 2006.
   ISBN 978-0-19-856699-1 978-0-19-856702-8.
.. [3] Z. Yang.
   Molecular Evolution: A Statistical Approach.
   Oxford University Press, Oxford, United Kingdom; New York, NY, United
   States of America,  rst edition edition, 2014.
   ISBN 978-0-19-960261-2 978-0-19-960260-5.
"""

from asymmetree.seqevolve.Evolver import Evolver
from asymmetree.seqevolve.SubstModel import SubstModel
from asymmetree.seqevolve.IndelModel import IndelModel
from asymmetree.seqevolve.HetModel import HetModel
from asymmetree.seqevolve.Alignment import AlignmentBuilder
from asymmetree.seqevolve.EvolvingSequence import EvoSeq