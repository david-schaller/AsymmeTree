# -*- coding: utf-8 -*-

import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import (to_newick,)


__author__ = 'David Schaller'


print('Yule ------------------------')
tree = te.simulate_species_tree(10, model='yule', birth_rate=1.0)
print(to_newick(tree))

print('EBDP ------------------------')
tree2 = te.simulate_species_tree(10,
            episodes=[(1.0, 0.3, 0.8, 0.0), (0.9, 0.4, 0.6, 0.3)])
print(to_newick(tree2))

print('Yule age ------------------------')
tree3 = te.simulate_species_tree_age(2.0, model='yule', birth_rate=1.0)
print(to_newick(tree3))

print('EBDP age ------------------------')
tree4 = te.simulate_species_tree_age(2.0, model='EBDP', birth_rate=1.0,
            episodes=[(1.0, 0.3, 0.8, 0.0), (0.9, 0.4, 0.6, 0.3)])
print(to_newick(tree4))