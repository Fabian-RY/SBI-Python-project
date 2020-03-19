# -*- coding: utf-8 -*-

"""
Created on Tue Feb 25 19:39:13 2020

@author: fabian

@version: 1.0

"""

from . import superimpose
from .superimpose import Ensemble

import copy

from Bio import PDB

def _get_total_chains(*structures):
    '''
        Obtains the total number of different chains that are in the different structures
    '''
    chains = set()
    for sc in structures:
        for chain in sc.get_chains():
            chains.add(chain.id)
    return chains

def get_fastas_from_structs(structs, fastas):
    chain_fasta = {}
    for struct in structs:
        chain_fasta[struct.id] = {}
        for chain in struct.get_chains():
            for seq in fastas:
                if chain.id in seq.id:
                    chain_fasta[struct.id][chain.id] = seq
    return chain_fasta

def build_complex(threshold, distance, stoichiometry, sequences, structures, verbose):
    seqs = get_fastas_from_structs(structures, sequences)
    full_structure = PDB.Structure.Structure('full')
    full_structure.add(PDB.Model.Model('model'))
    iters = 0
    pair_num = len(structures)
    if(pair_num < 2):
        raise ValueError('Needed at least 2 pairs to superpose')
    comparisons = {}
    failed = []
    current = 0
    done = []
    while True:
        '''
        print(1, len(done), len(structures))
        print(2, failed)
        print(3, current)
        '''
        if(len(done) == len(structures)): break
        structure = structures[current % len(structures)]
        if(structure.id in done): 
            current += 1
            print('Skiping')
            continue # If is already in, continue to the next one
        if verbose:
            print('Chains joined', len(list(full_structure.get_chains()))) 
            print('Pairs left', len(structures))
            print('Current pair:', structure.id)
        count = list(failed.count(element) for element in failed)
        if 2 in count: # We are repeating structures
            if verbose: print('Some pdbs could not be joined')
            break
        if len(list(full_structure.get_chains())) == 0:
            print('Initializing structure')
            seqs[full_structure.id] = {}
            for idx, chain in enumerate(structure.get_chains()):
                print(idx, chain.id, end=' ')
                chain_id = 'ABCDEFGHIJKLMNÑOPQRSTUVWXYZ0123456789|@#~½¬{[]}'[idx]
                seqs[full_structure.id][chain_id] = seqs[structure.id][chain.id]
                chain_new = PDB.Chain.Chain(chain_id)
                chain_new.child_list = list(chain.get_residues())
                model = next(full_structure.get_models())
                model.child_list.append(chain_new)
            done.append(structure.id)
            current += 1
            continue
        print('Ensembling')
        pair = Ensemble(full_structure, structure)
        pair_seqs = (seqs[full_structure.id], seqs[structure.id])
        align = pair.get_best_alignment(pair_seqs)
        if verbose: print(align[0][0], align[0][1])
        # If score is > 0.95, they should be homologs
        print(align)
        for alignment in align:
            if(align[1] > threshold):
                atoms_of_chains = pair.superimpose(align[0][0], align[0][1])
                iters += 1
                neighbor = PDB.NeighborSearch(list(full_structure.get_atoms()))
                clashes = 0
                for chain in atoms_of_chains:
                    for atom in atoms_of_chains:
                        for atom in chain.get_atoms():
                            close_atoms = neighbor.search(atom.get_coord(), distance)
                            # If there are atoms within 2 angstroms, consider a clash
                            if len(close_atoms) > 0:
                                clashes += 1
                print(5, clashes)
                if (clashes < 10):
                    for chain in full_structure.get_chains():
                        print(chain.id, end=' ')
                    for chain in atoms_of_chains:
                        chain_id = 'ABCDEFGHIJKLMNÑOPQRSTUVWXYZ0123456789|@#~½¬{[]}'[len(list(full_structure.get_chains()))]
                        chain2 = PDB.Chain.Chain(chain_id)
                        chain2.child_list += list(chain.get_residues())
                        model = next(full_structure.get_models())
                        model.child_list.append(chain2)
                        failed = []
                    done.append(structure.id)
                    print('breaking')
                    break
                else:
                    print('Failed')
                    failed.append(structure.id)
                    comparisons[structure.id] = comparisons.get(structure.id, 0) +1
                    # Pdbs clash
        else: #If there are no good alignments, is count as a failed
            failed.append(structure.id)
        current += 1
    return full_structure