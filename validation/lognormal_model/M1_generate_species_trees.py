# -*- coding: utf-8 -*-


import asymmetree.treeevolve as te


filename = 'testfile_scenarios.csv'
number_of_trees = 10000
variances = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]


with open(filename, 'w') as f:
    
    f.write('tree_id,n,variance,rate')

    for i in range(number_of_trees):
        print(i)
        
        S = te.species_tree_age(1.0, birth_rate=3.0, model='yule')
        n = sum(1 for _ in S.leaves())
        
        for variance in variances:
            
            rates, _ = te.autocorrelation_factors(S, variance)
            
            for v in S.leaves():
                f.write(f'\n{i},{n},{variance},{rates[v.label]}')
