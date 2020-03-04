#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Feb 25 19:39:13 2020

@author: fabian

@version: 0.1

"""

from Bio import PDB

import os, os.path
import argparse
import sys

AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
      'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
      'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
      'MET':'M', 'PHE':'F','PRO':'P','SER':'S',
      'THR':'T', 'TRP':'Y','VAL':'V','TERM':''}

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
    arguments = _parse_args()
    if(os.path.isfile(arguments.input)):
        pdb_file = arguments.input
        PDBparser = PDB.PDBParser()
        structure = PDBparser.get_structure(pdb_file[:-3], pdb_file)
        io = PDB.PDBIO()
        if(os.path.exists(arguments.output) and os.path.isfile(arguments.output)):
            print('File {out} already exists!'.format(out=arguments.output), file=sys.stderr)
            sys.exit(-1)
        elif(not os.path.exists(arguments.output)): 
            os.mkdir(arguments.output)
        for chain in structure.get_chains():
            splitted = PDB.Structure.Structure(pdb_file[:-3]+chain.id)    
            splitted.add(chain)
            io.set_structure(chain)
            print('Saving chain {chain_id} in {prefix}_{chain_id}.pdb'.format(prefix=arguments.prefix, chain_id=chain.id))
            if (arguments.fasta):
                print('Saving sequence of chain {chain_id} in {prefix}_{chain_id}.fa'.format(prefix=arguments.prefix, chain_id=chain.id))
                filename = os.path.join(arguments.output, '{prefix}_{chain_id}.fa'.format(prefix=arguments.prefix, chain_id=chain.id))
                fhand = open(filename, 'w')
                fhand.write('>{prefix}_{chain_id}\n'.format(prefix=arguments.prefix, chain_id=chain.id))
                for residue in chain:
                    fhand.write(AA.get(residue.resname, ''))
                fhand.close()
            filename = os.path.join(arguments.output, '{prefix}_{chain_id}.pdb'.format(prefix=arguments.prefix, chain_id=chain.id))
            io.save(filename)
    elif(not os.path.exists(arguments.input)):
        print('File {file} does not exist'.format(file=arguments.input), 
              file=sys.stderr)
    else:
        print('Unknown input file or is not a file', 
              file=sys.stderr)
        