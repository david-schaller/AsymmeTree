# -*- coding: utf-8 -*-

import os
import asymmetree.treeevolve as te


directory = 'results'
repeats = 10000
n = 20


if not os.path.exists(directory):
    os.mkdir(directory)
    

death_rates = [0.0, 0.2, 0.5, 0.8]

for d in death_rates:
    
    ages = []
    
    for i in range(repeats):
        
        if i % 1000 == 0:
            print(d,i)
        
        model = 'yule' if d == 0.0 else 'BDP'
        
        tree = te.species_tree_n(n, model=model,
                                 birth_rate=1.0, death_rate=d)
        
        ages.append(tree.root.tstamp)
    
    filename = os.path.join(directory, f'py_age{d}.txt')
    with open(filename, 'w') as f:
        for age in ages:
            f.write(f'{age}\n')
