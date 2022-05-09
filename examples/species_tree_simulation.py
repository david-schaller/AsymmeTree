# -*- coding: utf-8 -*-

import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import (to_newick,)


__author__ = 'David Schaller'


print('Yule n ------------------------')
tree = te.species_tree_n(10, model='yule', birth_rate=1.0)
print(to_newick(tree))

print('EBDP n ------------------------')
tree2 = te.species_tree_n(10, model='EBDP',
                          episodes=[(1.0, 0.3, 0.8, 0.0),
                                    (0.9, 0.4, 0.6, 0.3)])
print(to_newick(tree2))



print('Yule age ------------------------')
tree3 = te.species_tree_age(2.0, model='yule', birth_rate=1.0)
print(to_newick(tree3))

print('EBDP age ------------------------')
tree4 = te.species_tree_age(2.0, model='EBDP', birth_rate=1.0,
                            episodes=[(1.0, 0.3, 0.8, 0.0),
                                      (0.9, 0.4, 0.7, 0.3)])
print(to_newick(tree4))



print('Yule n age ------------------------')
tree5 = te.species_tree_n_age(10, 1.0, model='yule', birth_rate=1.0)
print(to_newick(tree5))

print('EBDP n age ------------------------')
tree6 = te.species_tree_n_age(10, 1.0, model='BDP',
                              birth_rate=1.0,
                              death_rate=0.5)
print(to_newick(tree6))
