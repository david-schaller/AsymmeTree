# -*- coding: utf-8 -*-

import asymmetree.treeevolve as te


__author__ = 'David Schaller'


print('Yule ------------------------')
tree = te.simulate_species_tree(10, model='yule', birth_rate=1.0)
print(tree.to_newick())

print('EBDP ------------------------')
tree2 = te.simulate_species_tree(10,
            episodes=[(1.0, 0.3, 0.8, 0.0), (0.9, 0.4, 0.6, 0.3)])
print(tree2.to_newick())

print('Yule age ------------------------')
tree3 = te.simulate_species_tree_age(2.0, model='yule', birth_rate=1.0)
print(tree3.to_newick())

print('EBDP age ------------------------')
tree4 = te.simulate_species_tree_age(2.0, model='EBDP', birth_rate=1.0,
            episodes=[(1.0, 0.3, 0.8, 0.0), (0.9, 0.4, 0.6, 0.3)])
print(tree4.to_newick())