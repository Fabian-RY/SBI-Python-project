#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:11:47 2020

@author: fabian

Executable script of our SBI-PYT project
"""

from Bio import PDB
from Bio import SeqIO
from argparse import ArgumentParser

import sys
import re
import os


from builder import build_macrocomplex
from Bio import PDB

def _parse_args():
    '''
        Parses the arguments whenever is called using command-line interface. Mandatory arguments are:
            -i: the folder with the input pdbs
            -o: the folder where the output pdbs will be stored
            -f: the fasta file with the sequences
            
        Optional arguments are:
            -v: Shows progress log
            -s: indicates the stoichiometry of the proteins
            -d: maximun distance (in Armstrongs) to consider that 2 atoms clash. Default 5 Armstrong
            -t: threshold to consider whether two sequences are homologous or not
            -o: After finishing, calls modeller to optimize energies of the model
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
    parser.add_argument('-f', '--fasta',
                        action='store',
                        dest='fasta',
                        required=True)
    parser.add_argument('-v','--verbose',
                        action='store_true',
                        dest='log')
    parser.add_argument('-s','--stoichiometry',
                        action='store',
                        dest='stoichiometry',
                        help='The stoichiometry of the final molecule the builder will try to forme',
                        default='')
    parser.add_argument('-d','--distance',
                        action='store',
                        dest='distance',
                        help='Maximun distance to consider that two atoms clash, in Armstrongs. By default, 0.1 armstrong',
                        type=float,
                        default=0.1)
    parser.add_argument('-t','--threshold',
                        action='store',
                        help='Minimum score to consider two sequences homologus, and thus, the same to be build. Default 0.95',
                        dest='threshold',
                        type=float,
                        default=0.95)
    parser.add_argument('-oo','--optimize',
                        action='store_true',
                        help='Whether the model will be optimized with MODELLER after building or not. Default False',
                        dest='optimize',
                        default=False)
    return parser.parse_args()

def _parse_stoichiometry(input_arg):
    '''
        Input separated by ;
        A2;B3;C2 
        
        Last element must not have ;
    '''
    if(not input_arg): return input_arg
    splitted = input_arg.split(';')
    for element in splitted:
        result = re.search('^([A-Za-Z]+)([0-9]+)$', element)
        if (not result): 
            raise Exception('Non valid stoichimetry')
            

if __name__ == '__main__':
    
    
    ##########################################
    #     Parse input and stoichiometry      #
    ##########################################
    arguments = _parse_args()
    parser = PDB.PDBParser(QUIET=1)  
    structures = list()
    print('Reading PDB files. This may take a while...')
    for file in os.listdir(arguments.input_folder):
        if file.endswith('.pdb'):
            path = os.path.join(arguments.input_folder, file[:-4])
            pdb = parser.get_structure(path, path+'.pdb')
            structures.append(pdb)
    sequences = SeqIO.parse(arguments.fasta, 'fasta')
    stoic = _parse_stoichiometry(arguments.stoichiometry) 
    
    print('Building complex...')
    #############################################
    #          Build the complex                #
    #############################################
    try:
        model = build_macrocomplex.build_complex(threshold=arguments.threshold,
                                                 stoichiometry=stoic, 
                                                 sequences=list(sequences),
                                                 structures=structures,
                                                 distance=arguments.distance,
                                                 verbose=arguments.log)
    except ValueError:
        print('There are not enough pairs of pdbs to superimpose. You need at least two pdb files')
    
    ######################################################
    # Saving the model into a file in the indicated path #
    ######################################################
    print('Saving complex...')
    io = PDB.PDBIO()
    io.set_structure(model)
    io.save(os.path.join(arguments.output_folder, 'final_model.pdb'))
    
    ###############################
    # Try to optimize the build   #
    ###############################
    if(arguments.optimize):
        try:
            import modeller
        except ImportError:
            print('Modeller could not be found, so no optimization will be done. Please, install it before using the --optimize option', 
                  file=sys.stderr)
        else:
            print('Optimizing...')
    print('Model completed')