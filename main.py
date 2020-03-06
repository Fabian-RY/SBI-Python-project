#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:11:47 2020

@author: fabian

Executable script of our SBI-PYT project
"""

from Bio import PDB
from argparse import ArgumentParser

from . import builder

def _parse_args():
    '''
        Parses the arguments whenever is called using command-line interface. Mandatory arguments are:
            -i: the folder with the input pdbs
            -o: the folder where the output pdbs will be stored
            
        Optional arguments are:
            -v: Shows progress log
            -s: indicates the stoichiometry of the proteins
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
                        dest='stoichiometry')
    return parser.parse_args()

if __name__ == '__main__':
    arguments = _parse_args()
    pass
    