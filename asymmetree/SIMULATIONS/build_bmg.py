# -*- coding: utf-8 -*-

import os, glob

from proteinortho.POBestMatches import POBestMatches, print_BMG
    
directory = "proteinortho/test_files_2"

fasta_files = glob.glob(os.path.join(directory, "*.faa"))
candidate_dir = os.path.join(directory, "test.bm_candidates")
tree_file = os.path.join(directory, "true_species_tree")

po_bm = POBestMatches(fasta_files, candidate_dir, tree_file)
BMG = po_bm()

print_BMG(BMG, os.path.join(directory, "bmg"))
