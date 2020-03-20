#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:08:28 2020

@author: fabian
""" 

import unittest
import re
import os

from Bio import PDB
from Bio import SeqIO

from builder import superimpose as sp
from builder import build_macrocomplex as bm

def create_class():
    pdb_file = 'examples/3e0d/pairs/3e0d_AC.pdb'
    pdb_file2 = 'examples/3e0d/pairs/3e0d_AD.pdb'
    fa_file1 = 'examples/3e0d/pairs/3e0d_AD.fa'
    fa_file2 = 'examples/3e0d/pairs/3e0d_AD.fa'
    parser = PDB.PDBParser(QUIET=1)
    structure = parser.get_structure(pdb_file[:-3], pdb_file) # Parses the PDB
    structure2 = parser.get_structure(pdb_file2[:-3], pdb_file2) # Parses the PDB
    return sp.Ensemble(structure, structure2), structure, structure2

class test(unittest.TestCase):
           
    def test_recursive(self):
        parser = PDB.PDBParser(QUIET=1)  
        structure = parser.get_structure('good','examples/5nss.pdb')
        print(len(list(structure.get_chains())))
        structures = list()
        path = 'examples/5nss/pairs'
        fasta = '5nss.fa'
        for file in os.listdir(path):
            if(not os.path.isfile(os.path.join(path, file))): continue
            if(not file.endswith('.pdb')): continue
            path_file = os.path.join(path, file)
            pdb = parser.get_structure(path_file[:-3], path_file)
            structures.append(pdb)
        threshold=0.9
        distance=1
        stoichiometry=''
        model = bm.build_complex(threshold, distance, stoichiometry, 
                         structures=structures, 
                         sequences=list(SeqIO.parse(fasta, 'fasta')), 
                         verbose=False)
        io = PDB.PDBIO()
        io.set_structure(model)
        io.save('final_model.pdb')
        pass
    
if __name__ == '__main__':
    unittest.main()
