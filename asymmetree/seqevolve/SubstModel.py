# -*- coding: utf-8 -*-

import Matrices


class SubstModel:
    
    
    nuc_models = {"JC69", "K80"}
    aa_models = {"JC69"}
    
    nucs = "ACGT"
    nuc_dict = {nuc: index for index, nuc in enumerate(nucs)}
    
    aas = "ACDEFGHIKLMNPQRSTVWY"
    aa_dict = {aa: index for index, aa in enumerate(aas)}
    
    
    def __init__(self, model_type, model_name):
        
        model_type = model_type.lower()
        model_name = model_name.upper()
        
        if (model_type in ("nuc", "nucleotide") and
            model_name in SubstModel.nuc_models):
            
            self.model = ("nuc", model_name)
            
        elif (model_type in ("aa", "amino", "aminoacid", "protein") and
              model_name in SubstModel.aa_models):
            
            self.model = ("aa", model_name)
                
        else:
            raise ValueError("Model type '{}' is not available!".format(model_type))
            
        self._compute_Q()
        
    
    def _compute_Q(self):
        
        if self.model == ("nuc", "JC69"):
            S, self.freqs = Matrices._JC69_nuc()
        
        elif self.model == ("aa", "JC69"):
            S, self.freqs = Matrices._JC69_aa()
            
        self.Q = Matrices.build_matrix_Q(S, self.freqs)
        
    
    def to_indices(self, sequence):
        
        if self.model[0] == "nuc":
            alphabet_dict = SubstModel.nuc_dict
        else:
            alphabet_dict = SubstModel.aa_dict
        
        try:
            result = [alphabet_dict[letter] for letter in sequence]
        except KeyError:
            raise ValueError("Invalid sequence for the specified model!")
            
        return result
    
    
    def to_sequence(self, evoseq):
        
        if self.model[0] == "nuc":
            alphabet = SubstModel.nucs
        else:
            alphabet = SubstModel.aas
        
        return "".join(alphabet[x._value] for x in evoseq)
            
        