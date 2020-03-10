#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:08:28 2020

@author: fabian
"""

import unittest

from Bio import PDB
from Bio import SeqIO

from builder import superimpose as sp

def create_class():
    pdb_file = 'examples/3e0d/pairs/3e0d_AB.pdb'
    fa_file = 'examples/3e0d/pairs/3e0d_AB.fa'
    pdb_file2 = 'examples/3e0d/pairs/3e0d_AC.pdb'
    fa_file2 = 'examples/3e0d/pairs/3e0d_AC.fa'
    parser = PDB.PDBParser(QUIET=1)
    structure = parser.get_structure(pdb_file[:-3], pdb_file) # Parses the PDB
    structure2 = parser.get_structure(pdb_file2[:-3], pdb_file2) # Parses the PDB
    fasta_1 = SeqIO.parse(fa_file, 'fasta')
    fasta_2 = SeqIO.parse(fa_file2, 'fasta')
    return sp.Ensemble(structure, structure2, fasta_1, fasta_2), structure, structure2

class test(unittest.TestCase):
    
    def test_class(self):
        create_class()

    def test_align(self):
        ensembl, _, _ = create_class()
        alignment = ensembl.get_best_alignment()

    def test_chain_superimpose(self):
        ensembl, structure, structure2 = create_class()
        #alignment = ensembl.get_best_alignment()
        ensembl.superimpose('A')
        
        
    def test_save(self):
        ensembl, structure, structure2 = create_class()
        #alignment = ensembl.get_best_alignment()
        ensembl.superimpose('A')
        ensembl.save_model('final_model.pdb', 'A')
    
if __name__ == '__main__':
    unittest.main()