#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Feb 25 19:39:13 2020

@author: fabian

@version: 1.0

Splits a PDB into its component chains, each in one file. If indicated, 
it will also create fasta files for each chain

"""

from Bio import PDB

import os, os.path
import argparse
import sys

AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
      'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
      'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
      'MET':'M', 'PHE':'F','PRO':'P','SER':'S',
      'THR':'T', 'TRP':'W','VAL':'V','TYR':'Y', 'TERM':'',
      'UNK':'X',' DA':'A',' DG':' G', ' DC':'C', 
      ' DT':'T','  U':'U','  A':'A','  G':'G', '  C':'C'}

def _parse_args():
    '''
        Parses the arguments needed for this script to work
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True, dest='input', action='store')
    parser.add_argument('-p','--prefix', required=True, dest='prefix', action='store')
    parser.add_argument('-o','--output', required=True, dest='output', action='store')
    parser.add_argument('-f','--fasta', action='store_true', dest='fasta', default=False)
    return parser.parse_args()


if(__name__ == '__main__'):
    '''
        Main routine of the program
        
        Checks the arguments to know if the input is correct. Then tries to parse the pdb structure indicated in arguments.input
        (may throw a Bio.PDB warning though). If output folder is a file, exits with error -1. If not, create the folder in case it doesn't exists and continues.
        
        For every chain of the complex, creates a new empty structure object, adds the chain and saves it with the estructure {prefix}_{chain_id}.pdb.
        In case the flag -f was indicate, also saves the fasta sequence in {prefix}_{chain_id}.fa. Both files are saved in the output directory
    '''
    arguments = _parse_args()
    if(os.path.isfile(arguments.input)):  # If the input is a file, tries to read it
        pdb_file = arguments.input
        PDBparser = PDB.PDBParser(QUIET=True)
        structure = PDBparser.get_structure(pdb_file[:-3], pdb_file) # Parses the PDB
        io = PDB.PDBIO()
        if(os.path.exists(arguments.output) and os.path.isfile(arguments.output)): #Checks the output is correct and either doesn't exists or it's a folder
            print('File {out} already exists!'.format(out=arguments.output), file=sys.stderr)
            sys.exit(-1) # Exits because output is a file and already exists
        elif(not os.path.exists(arguments.output)): # Folder does not exist, so let's create it
            os.mkdir(arguments.output)
        # For each chain save it in a separate file. Log is provided for very long proteins
        for chain in structure.get_chains():
            splitted = PDB.Structure.Structure(pdb_file[:-3]+chain.id)
            splitted.add(chain)
            io.set_structure(chain)
            print('Saving chain {chain_id} in {prefix}_{chain_id}.pdb'.format(prefix=arguments.prefix, chain_id=chain.id))
            # Saves also the fasta file if indicated with -f flag
            if (arguments.fasta):
                print('Saving sequence of chain {chain_id} in {prefix}_{chain_id}.fa'.format(prefix=arguments.prefix, chain_id=chain.id))
                filename = os.path.join(arguments.output, '{prefix}_{chain_id}.fa'.format(prefix=arguments.prefix, chain_id=chain.id))
                fhand = open(filename, 'w')
                fhand.write('>{prefix}_{chain_id}\n'.format(prefix=arguments.prefix, chain_id=chain.id))
                for residue in chain:
                    fhand.write(AA.get(residue.resname, '').replace(' ',''))
                fhand.write('\n')
                fhand.close()
            filename = os.path.join(arguments.output, '{prefix}_{chain_id}.pdb'.format(prefix=arguments.prefix, chain_id=chain.id))
            io.save(filename)
    elif(not os.path.exists(arguments.input)): # If there's no input file
        print('File {file} does not exist'.format(file=arguments.input), 
              file=sys.stderr)
    else: # User gave a folder as input
        print('Unknown input file or is not a file', 
              file=sys.stderr)
        