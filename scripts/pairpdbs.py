#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 12:11:19 2020

@author: fabian

@version: 0.9

Joins all the PDBs in a folder into pairs, but only if they interact
"""

from Bio import PDB
import os
import os.path
import sys
import re

def join_fasta(*fastas):
    text = ''
    for fhand in fastas:
        with open(fhand) as f1:
            for line in f1:
                text += line 
    return text

if __name__ == '__main__':
    directory = sys.argv[1] # Folder with the chain pdbs
    distance = float(sys.argv[2])
    directory_out = os.path.join(directory, 'pairs') # Folder where pairs will be saved
    PDBparser = PDB.PDBParser(QUIET=True)
    done = []
    for pdb_1 in list(filter(lambda x: x.endswith('.pdb'), os.listdir(directory))):
        for pdb_2 in list(filter(lambda x: x.endswith('.pdb'), os.listdir(directory))):
            if (pdb_1 == pdb_2):   
                continue # Not to duplicate itself
            elif not (pdb_1, pdb_2) in done or not (pdb_2, pdb_1) in done: # Not to duplicate pairs
                done.append((pdb_1, pdb_2))
                structure1 = PDBparser.get_structure(pdb_1[:-4], os.path.join(directory,pdb_1))
                structure2 = PDBparser.get_structure(pdb_2[:-4], os.path.join(directory,pdb_2))
                for chain in structure2.get_chains():
                    atoms = list(chain.get_atoms()) 
                    ns = PDB.NeighborSearch(atoms) # An object to search chains near an atom 
                    for target_atom in structure1.get_atoms():
                        near = ns.search(target_atom.coord, distance)
                        if(near and (pdb_1[-5:-4],pdb_2[-5:-4]) not in done and (pdb_2[-5:-4],pdb_1[-5:-4]) not in done): # If there's an atom near, they interact
                            print('Joining chain {a} and {b}'.format(a=pdb_1[-5:-4], b=pdb_2[-5:-4]))
                            done.append((pdb_1[-5:-4],pdb_2[-5:-4]))
                            io = PDB.PDBIO()
                            structure1[0].child_list += [chain for chain in structure2.get_chains()] # Joins all the chains
                            io.set_structure(structure1)
                            if(not os.path.exists(os.path.join(directory_out))):
                                os.mkdir(os.path.join(directory_out))
                            p = re.search('^[A-Za-z0-9]+_([A-Za-z0-9]+).pdb$', pdb_1 ).group(1)
                            q = re.search('^[A-Za-z0-9]+_([A-Za-z0-9]+).pdb$', pdb_2 ).group(1)
                            io.save(os.path.join(directory_out,pdb_1[:-5]+p+q+'.pdb'))
                            print(os.path.join(directory,pdb_1[:-4]+'.fa'))
                            if(os.path.exists(os.path.join(directory,pdb_1[:-4]+'.fa'))) and os.path.exists(os.path.join(directory,pdb_2[:-4]+'.fa')):
                                path = os.path.join(directory_out,pdb_1[:-5]+p+q+'.fa')
                                file = open(path, 'w')
                                file.write(join_fasta(os.path.join(directory,pdb_1[:-4]+'.fa'), os.path.join(directory,pdb_2[:-4]+'.fa')))
                                file.close()
