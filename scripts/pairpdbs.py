#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 12:11:19 2020

@author: fabian

@version: 0.9

Joins all the PDBs in a folder into pairs 
"""

from Bio import PDB
import os
import os.path
import sys

if __name__ == '__main__':
    directory = sys.argv[1]
    PDBparser = PDB.PDBParser()
    done = []
    for pdb_1 in os.listdir(directory):
        for pdb_2 in os.listdir(directory):
            if (pdb_1 == pdb_2): continue
            if not (pdb_1, pdb_2) in done or (pdb_2, pdb_1) in done:
                done.append((pdb_1, pdb_2))
                print('Joining chain {a} and {b}'.format(a=pdb_1[-5:-4], b=pdb_2[-5:-4]))
                structure1 = PDBparser.get_structure(pdb_1[:-4], os.path.join(directory,pdb_1))
                structure2 = PDBparser.get_structure(pdb_2[:-4], os.path.join(directory,pdb_2))
                io = PDB.PDBIO()
                structure1[0].child_list += [chain for chain in structure2.get_chains()]
                io.set_structure(structure1)
                io.save(os.path.join(directory,pdb_1[:-4]+pdb_2[-5:-4]+'.pdb'))