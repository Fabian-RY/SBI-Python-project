#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:11:47 2020

@author: fabian
"""

from Bio import PDB
from argparse import ArgumentParser

def _parse_args():
    parser = ArgumentParser
    

if __name__ == '__main__':
    parser = PDB.PDBParser()
    struc1 = parser.get_structure('6om3_A', 'examples/6om3/split_chain/6om3_A.pdb')
    atom1 = []
    for model in struc1:
        for chain in struc1:
            for residue in struc1:
                for atom in residue:
                    atom1.append(atom)
    struc2 = parser.get_structure('6om3_B', 'examples/6om3/split_chain/6om3_B.pdb')
    atom2 = []
    for model in struc2:
        print(model)
        for chain in struc1:
            print(chain)
            for residue in struc1:
                print(residue)
                for atom in residue:
                    atom1.append(atom)
    sup = PDB.Superimposer()
    sup.set_atoms(atom1, atom2)
    sup.apply(atom2)
    