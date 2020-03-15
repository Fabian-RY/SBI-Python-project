#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:11:47 2020

@author: fabian

Executable script of our SBI-PYT project
"""

from Bio import PDB
from argparse import ArgumentParser
import re
import os


from builder import build_macrocomplex

def _parse_args():
    '''
        Parses the arguments whenever is called using command-line interface. Mandatory arguments are:
            -i: the folder with the input pdbs
            -o: the folder where the output pdbs will be stored
            
        Optional arguments are:
            -v: Shows progress log
            -s: indicates the stoichiometry of the proteins
            -d: maximun distance (in Armstrongs) to consider that 2 atoms clash. Default 5 Armstrong
            -
    '''
    parser = ArgumentParser('Build a macromolecular complex using interacting subcomponents')
    parser.add_argument('-i', '--input-folder',
                        action='store',
                        dest='input_folder',
                        required=True)
    parser.add_argument('-o', '--output-folder',
                        action='store',
                        dest='output_folder',
                        required=True)
    parser.add_argument('-v','--verbose',
                        action='store_true',
                        dest='log')
    parser.add_argument('-s','--stoichiometry',
                        action='store',
                        dest='stoichiometry',
                        default=None)
    parser.add_argument('-d','--depth',
                        action='store',
                        dest='depth',
                        default=-1)
    return parser.parse_args()

def _parse_stoichiometry(input_arg):
    '''
        Input separated by ; and:
        A2;B3;C2
    '''
    splitted = input_arg.split(';')
    for element in splitted:
        result = re.search('^[A-Za-Z]+[0-9]+$', element)
        if (not result): 
            raise Exception('Non valid stoichimetry')
            

if __name__ == '__main__':
    arguments = _parse_args()
    parser = PDB.PDBParser(QUIET=1)  
    structures = list()
    for file in os.listdir(arguments.input_folder):
        if file.endswith('.pdb':)
            path = os.path.join(arguments.input_folder, file[:-4])
            pdb = parser.get_structure(path, path+'.pdb')
            pdb = parser.get_structure(path, path+'.pdb')
        structures.append(pdb)
    structures.sort()
    model = build_macrocomplex.build_complex(arguments.stoichiometry, *structures)
    