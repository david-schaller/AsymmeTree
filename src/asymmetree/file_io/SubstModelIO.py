# -*- coding: utf-8 -*-


__author__ = 'David Schaller'


def parse_paml(filename, model_type='a'):
    
    if model_type == 'a':
        N = 20
    elif model_type == 'n':
        N = 4
    else:
        raise ValueError("model type '{}' not supported".format(model_type))
    
    with open(filename, 'r') as f:
        
        exchangeability_matrix = [[0.0 for j in range(N)] for i in range(N)]
        
        line = f.readline().strip()
        while line == '': line = f.readline().strip()
        
        for i in range(1,N):
            
            line = [float(item) for item in line.split()]
            
            for j in range(len(line)):
                exchangeability_matrix[i][j] = line[j]
                exchangeability_matrix[j][i] = line[j]
                
            line = f.readline().strip()
                
        line = f.readline().strip()
        while line == '': line = f.readline().strip()
            
        stat_freqs = [float(item) for item in line.split()]
        
        if len(stat_freqs) != N:
            raise RuntimeError('wrong no. of equilibrium frequencies, '\
                               'check paml file')
        
        return exchangeability_matrix, stat_freqs