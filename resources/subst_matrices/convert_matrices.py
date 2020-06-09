# -*- coding: utf-8 -*-

from glob import glob

from asymmetree.file_io.SubstModelIO import parse_paml
    

def width_max(a):
    
    maximum = 0
    
    for i in range(len(a)):
        
        if isinstance(a[i], list):
            for j in range(len(a)):
                length = len(str(a[i][j]))
                if length > maximum:
                    maximum = length
        else:
            length = len(str(a[i]))
            if length > maximum:
                maximum = length
    
    return maximum


def format_list(a, max_length):
    
    bracket = ', {:' + str(max_length) + '}'
    format_string = '[' + bracket[2:] + (len(a)-1) * bracket + ']'
    
    return format_string.format(*a)


def write_python_function(f, model, exchangeability_matrix, stat_freqs):
    
    em_max = width_max(exchangeability_matrix)
    sf_max = width_max(stat_freqs)
    
    f.write('def _{}_aa():\n    '.format(model))
    
    f.write('\n    S = np.array([' + format_list(exchangeability_matrix[0], em_max))
    for i in range(1, len(exchangeability_matrix)):
        f.write(',\n                  ' + format_list(exchangeability_matrix[i], em_max))
        
    f.write('])\n    \n    freqs = np.array(' + format_list(stat_freqs, sf_max) + ')\n    ')
    f.write('\n    return S, freqs')
    

outfile = '../../asymmetree/seqevolve/EmpiricalModels.py'

with open(outfile, 'w') as f:
    
    f.write('# -*- coding: utf-8 -*-\n\nimport numpy as np')
            
    models = []
    
    for filename in glob('*.paml'):
        
        model = filename[:filename.find('.')].upper()
        print(model)
        models.append(model)
        exchangeability_matrix, stat_freqs = parse_paml(filename, model_type='a')
        
        f.write('\n\n\n')
        write_python_function(f, model, exchangeability_matrix, stat_freqs)
        
    f.write('\n\n\n')
    f.write('empirical_models = {')
    for model in models:
        f.write("\n                    '{}': _{}_aa,".format(model, model))
    f.write('\n                   }')
    