# -*- coding: utf-8 -*-

"""
Created on Tue Feb 25 19:39:13 2020

@author: fabian

@version: 1.0

"""

from superimpose import Ensemble

from Bio import PDB
from Bio import pairwise2

AA_EQUIVALENCE = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
      'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
      'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
      'MET':'M', 'PHE':'F','PRO':'P','SER':'S',
      'THR':'T', 'TRP':'Y','VAL':'V', 'TYR':'Y','TERM':'',
      'UNK':'X','DA':'A','DG':'G', 'DC':'C', 
      'DT':'T','U':'U','A':'A','G':'G', 'C':'C'}

def _get_total_chains(*structures):
    '''
        Obtains the total number of different chains that are in the different structures
    '''
    chains = set()
    for sc in structures:
        for chain in sc.get_chains():
            chains.add(chain.id)
    return chains

def get_fastas_from_structs(structs, fastas, threshold=0.9):
    chain_fasta = {}
    for struct in structs:
        chain_fasta[struct.id] = {}
        for chain in struct.get_chains():
            chain_sequence = ''
            for residue in chain.get_residues():
                residue = residue.resname.lstrip().strip()
                print(residue, end='')
                chain_sequence += AA_EQUIVALENCE[residue]
            #print(chain_sequence)
            for seq in fastas:
                #print(seq.seq)
                print(chain_sequence, seq.seq)
                alignment = pairwise2.align.globalxx(chain_sequence, seq.seq, one_alignment_only=True)
                score = alignment[0][2]/len(alignment[0][0])
                print(score)
                if score > threshold:
                    chain_fasta[struct.id][chain.id] = seq
                    break
            else:
                raise Exception('PDB does not correspond to fasta')
                pass
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
        if(len(done) == len(structures)): 
            print('done')
            break
        structure = structures[current % len(structures)]
        if(structure.id in done): 
            current += 1
            print('Skiping')
            print(failed)
            continue # If is already in, continue to the next one
        if verbose:
            print('Chains joined', len(list(full_structure.get_chains()))) 
            print('Pairs left', len(structures))
            print('Current pair:', structure.id)
        count = list(failed.count(element) for element in failed)
        print(count)
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
        alignment = pair.get_best_alignment(pair_seqs)
        # If score is > 0.95, they should be homologs
        for align in alignment:
            if(align[1] > threshold):
                #try:
                atoms_of_chains = pair.superimpose(align[0][0], align[0][1])
                #except PDB.PDBExceptions.PDBException:
                #    current += 1
                #    continue
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
                    break
        else:
            print('Failed')
            failed.append(structure.id)
            comparisons[structure.id] = comparisons.get(structure.id, 0) +1
            # Pdbs clash
        current += 1
        print('Loop')
    return full_structure