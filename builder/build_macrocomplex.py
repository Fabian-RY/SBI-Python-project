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

def build_complex(threshold, repeated_limit, structures, sequences):
    seqs = get_fastas_from_structs(structures, sequences)
    full_structure = PDB.Structure.Structure('full')
    full_structure.add(PDB.Model.Model('model'))
    iters = 0
    pair_num = len(structures)
    if(pair_num < 2):
        raise ValueError('Needed at least 2 pairs to superpose')
    comparisons = {}
    failed_in_a_row = 0
    print(structures)
    failed = []
    for structure_bak in structures:
        print('Chains joined', len(list(full_structure.get_chains())))
        print('Pairs left', len(structures))
        print('Current pair:', structure_bak.id)
        structure = copy.deepcopy(structure_bak)
        count = list(failed.count(element) for element in failed)
        if 2 in count: # We are repeating structures
            print('Some pdbs could not be joined')
            break
        if len(list(full_structure.get_chains())) == 0:
            print('Initializing structure')
            seqs[full_structure.id] = {}
            for idx, chain in enumerate(structure.get_chains()):
                chain_id = 'ABCDEFGHIJKLMNO'[idx]
                seqs[full_structure.id][chain_id] = seqs[structure.id][chain.id]
                chain.id = chain_id
                model = next(full_structure.get_models())
                chain.id = chain_id
                model.child_list.append(chain)
            continue
        pair = Ensemble(full_structure, structure)
        pair_seqs = (seqs[full_structure.id], seqs[structure.id])
        align = pair.get_best_alignment(pair_seqs)
        print(align[0][0], align[0][1])
        # If score is > 0.95, they should be homologs
        if(align[1] > threshold):
            atoms_of_chains = pair.superimpose(align[0][0], align[0][1])
            iters += 1
            neighbor = PDB.NeighborSearch(list(full_structure.get_atoms()))
            clashes = 0
            for chain in atoms_of_chains:
                for atom in atoms_of_chains:
                    for atom in chain.get_atoms():
                        close_atoms = neighbor.search(atom.get_coord(), 1)
                        # If there are atoms within 2 angstroms, consider a clash
                        if len(close_atoms) > 0:
                            clashes += 1
            if (clashes < 10):
                for chain in atoms_of_chains:
                    failed_in_a_row = 0
                    chain.id = 'ABCDEFGHIJKLMNÑOPQRSTUVWXYZ0123456789|@#~½¬{[]}'[len(list(full_structure.get_chains()))]
                    model = next(full_structure.get_models())
                    model.child_list += [chain]
                    failed = []
            else:
                print('Failed')
                failed.append(structure.id)
                comparisons[structure.id] = comparisons.get(structure.id, 0) +1
                # Pdbs clash
                structures += [structure_bak]
    else:
        print('Ended for')
    io = PDB.PDBIO()

    io.set_structure(full_structure)
    io.save('final_model.pdb')