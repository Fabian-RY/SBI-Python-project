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

import builder
from builder import build_macrocomplex, errors

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
                        default=1)
    parser.add_argument('-t','--threshold',
                        action='store',
                        help='Minimum score to consider two sequences homologus, and thus, the same to be build. Default 0.95',
                        dest='threshold',
                        type=float,
                        default=0.95)
    parser.add_argument('-optimize','--optimize',
                        action='store_true',
                        help='Whether the model will be optimized with MODELLER after building or not. Default False',
                        dest='optimize',
                        default=False)
    parser.add_argument('-start','--start',
                        action='store',
                        help='Indicate the initial pdb from which the protein will be assembled',
                        dest='start',
                        default='')
    return parser.parse_args()

def _parse_stoichiometry(input_file):
    '''
        Input file with one element in the fasta per line, stoichiometry indicated by commas (,)
    '''
    if(not input_file): return input_file
    stoic = {}
    for line in open(input_file):
        result = re.search('^(.+),([0-9]+)$', line.strip())
        if (not result): 
            raise errors.stoichiometry_error(line.strip())
        stoic[result.group(1)] = int(result.group(2))
    return stoic

if __name__ == '__main__':
    
    
    ##########################################
    #     Parse input and stoichiometry      #
    ##########################################
    arguments = _parse_args()
    parser = PDB.PDBParser(QUIET=1)  
    structures = list()
    try:
        stoic = _parse_stoichiometry(arguments.stoichiometry)
    except errors.stoichiometry_error as e:
        print('Error: Stoichiometry file has an invalid format in line \'%s\'.' % e.error, file=sys.stderr)
        print('Aborting', file=sys.stderr)
        sys.exit(-1)
    except FileNotFoundError as e:
        print('Error: stoichiometry file does not exist')
        print('Aborting')
        sys.exit(-1)
    if(not os.path.exists(arguments.input_folder)):
        print('Error: Input folder doesn\'t exist', file=sys.stderr)
        print('Aborting', file=sys.stderr)
        sys.exit(1)
    print('Reading PDB files. This may take a while...')
    for file in os.listdir(arguments.input_folder):
        if file.endswith('.pdb'):
            if(arguments.start and os.path.join(arguments.input_folder, file) == arguments.start): continue 
            path = os.path.join(arguments.input_folder, file[:-4])
            pdb = parser.get_structure(path, path+'.pdb')
            structures.append(pdb)
    sequences = SeqIO.parse(arguments.fasta, 'fasta')
    print('Building complex...')
    if(arguments.start):
        start = parser.get_structure('initial', arguments.start)
        print('Established %s as the initial structure' % arguments.start)
        structures = [start] + structures
    initial = ''
    #############################################
    #          Build the complex                #
    #############################################
    try:
        model = build_macrocomplex.build_complex(threshold=arguments.threshold,
                                                 stoichiometry=stoic, 
                                                 sequences=list(sequences),
                                                 structures=structures,
                                                 distance=arguments.distance,
                                                 verbose=arguments.log,
                                                 initial=initial)
    except builder.errors.PDB_disagrees_fasta as e:
        print('The sequences in the PDB %s are not present in the given fasta' % (e.pdb))
        sys.exit(-1)
    except ValueError:
        print('There are not enough pairs of pdbs to superimpose. You need at least two pdb files')
        print('Aborting')
        sys.exit(-1)
    except builder.errors.chain_in_stoic_not_in_fasta as e:
        print('Conflict in stoichiometry and fasta file. %s is not present in one of those files' % e.seq)
        print('Aborting')
        sys.exit(-1)
    
    
    ######################################################
    # Saving the model into a file in the indicated path #
    ######################################################
    print('Saving complex...')
    io = PDB.PDBIO()
    io.set_structure(model)
    if(not os.path.exists(arguments.output_folder)):
        os.mkdir(arguments.output_folder)
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
        else: #Optimize
            from builder import optimize
            print('Optimizing...')
            sys.stdout = open(os.devnull, 'w')
            optimize.optimize(os.path.join(arguments.output_folder, 'final_model.pdb'), arguments.output_folder)
            sys.stdout = sys.__stdout__
            print('Done')
    print('Model completed')
    sys.exit(0)