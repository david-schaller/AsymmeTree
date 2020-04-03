# -*- coding: utf-8 -*-

from glob import glob

def parse_paml(filename):
    
    with open(filename, "r") as f:
        
        exchangeability_matrix = [[0.0 for j in range(20)] for i in range(20)]
        
        line = f.readline().strip()
        while line == "": line = f.readline().strip()
        
        for i in range(1,20):
            
            line = [float(item) for item in line.split()]
            
            for j in range(len(line)):
                exchangeability_matrix[i][j] = line[j]
                exchangeability_matrix[j][i] = line[j]
                
            line = f.readline().strip()
                
        line = f.readline().strip()
        while line == "": line = f.readline().strip()
            
        stat_freqs = [float(item) for item in line.split()]
        
        return exchangeability_matrix, stat_freqs
    

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
    
    bracket = ", {:" + str(max_length) + "}"
    format_string = "[" + bracket[2:] + (len(a)-1) * bracket + "]"
    
    return format_string.format(*a)


def write_python_function(f, model, exchangeability_matrix, stat_freqs):
    
    em_max = width_max(exchangeability_matrix)
    sf_max = width_max(stat_freqs)
    
    f.write("def _{}_aa():\n    ".format(model))
    
    f.write("\n    S = np.array([" + format_list(exchangeability_matrix[0], em_max))
    for i in range(1, len(exchangeability_matrix)):
        f.write(",\n                  " + format_list(exchangeability_matrix[i], em_max))
        
    f.write("])\n    \n    freqs = np.array(" + format_list(stat_freqs, sf_max) + ")\n    ")
    f.write("\n    return S, freqs")
    

outfile = "../../asymmetree/seqevolve/EmpiricalModels.py"

with open(outfile, "w") as f:
    
    f.write("# -*- coding: utf-8 -*-\n\nimport numpy as np")
            
    models = []
    
    for filename in glob("*.paml"):
        
        model = filename[:filename.find('.')].upper()
        print(model)
        models.append(model)
        exchangeability_matrix, stat_freqs = parse_paml(filename)
        
        f.write("\n\n\n")
        write_python_function(f, model, exchangeability_matrix, stat_freqs)
        
    f.write("\n\n\n")
    f.write("empirical_models = {")
    for model in models:
        f.write("\n                    '{}': _{}_aa,".format(model, model))
    f.write("\n                   }")
    