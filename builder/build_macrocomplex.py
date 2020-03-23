# -*- coding: utf-8 -*-

"""
Created on Tue Feb 25 19:39:13 2020

@author: fabian

@version: 1.0

"""

from .superimpose import Ensemble

from Bio import PDB
from Bio import pairwise2

def _chain_id():
    for i in range(65, 91, 1):
        yield chr(i)
    for i in range(97, 123, 1):
        yield i
    for i in range(48, 58, 1):
        yield i

AA_EQUIVALENCE = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
      'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
      'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
      'MET':'M', 'PHE':'F','PRO':'P','SER':'S',
      'THR':'T', 'TRP':'W','VAL':'V', 'TYR':'Y','TERM':'',
      'UNK':'X','DA':'A','DG':'G', 'DC':'C', 
      'DT':'T','U':'U','A':'A','G':'G', 'C':'C'}

def _get_chain_sequence(chain):
    '''
        Obtains the total number of different chains that are in the different structures
    '''
    chain_sequence = ''
    for residue in chain.get_residues():
        residue_val = residue.resname.lstrip().strip()
        chain_sequence += AA_EQUIVALENCE.get(residue_val,'')
    return chain_sequence            
    

def get_fastas_from_structs(structs, fastas, threshold=0.9):
    chain_fasta = {}
    for struct in structs:
        chain_fasta[struct.id] = {}
        for chain in struct.get_chains():
            chain_sequence = _get_chain_sequence(chain)
            #print(chain_sequence)
            for seq in fastas:
                #print(seq.seq)
                alignment = pairwise2.align.globalxx(chain_sequence, seq.seq, one_alignment_only=True)
                score = alignment[0][2]/len(alignment[0][0])
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
    ids = _chain_id()
    pair_num = len(structures)
    if(stoichiometry):
        current_number_of_chains = {chain_id:0 for chain_id in set(stoichiometry)}
    if(pair_num < 2):
        raise ValueError('Needed at least 2 pairs to superpose')
    comparisons = {}
    failed = []
    current = 0
    done = []
    while True:
        if(len(done) == len(structures)): 
            break
        structure = structures[current % len(structures)]
        if(structure.id in done): 
            current += 1
            continue # If is already in, continue to the next one
        count = list(failed.count(element) for element in failed)
        if 2 in count: # We are repeating structures
            if verbose: print('Some pdbs could not be joined')
            break
        if len(list(full_structure.get_chains())) == 0:
            seqs[full_structure.id] = {}
            for idx, chain in enumerate(structure.get_chains()):
                chain_id = next(ids)
                seqs[full_structure.id][chain_id] = seqs[structure.id][chain.id]
                chain_new = PDB.Chain.Chain(chain_id)
                chain_new.child_list = list(chain.get_residues())
                model = next(full_structure.get_models())
                if(stoichiometry):
                    chain_sequence = _get_chain_sequence(chain_new)
                    for seq in sequences:
                        alignment = pairwise2.align.globalxx(chain_sequence, seq.seq, one_alignment_only=True) 
                        # This chain corresponds to this stoichiometry chain
                        if(alignment[0][2]/len(alignment[0][0])> threshold and
                           current_number_of_chains[seq.id] < stoichiometry[seq.id]):
                            model.child_list.append(chain_new)
                            current_number_of_chains[seq.id] += 1
                            if(stoichiometry):
                                count = 0
                                for chain in current_number_of_chains:
                                    if current_number_of_chains[chain] == stoichiometry[chain]:
                                        count += 1
                                    if count == len(stoichiometry):
                                        print('Stoichiometry fullfilled!')
                                        return full_structure
                                    print(current_number_of_chains)
                model.child_list.append(chain_new)

            done.append(structure.id)
            current += 1
            continue
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
                if (clashes < 10):
                    for chain in atoms_of_chains:
                        chain_id = next(ids)
                        chain2 = PDB.Chain.Chain(chain_id)
                        chain2.child_list += list(chain.get_residues())
                        model = next(full_structure.get_models())
                        if(stoichiometry):
                            chain_sequence = _get_chain_sequence(chain)
                            for seq in sequences:
                                alignment = pairwise2.align.globalxx(chain_sequence, seq.seq, one_alignment_only=True) 
                                # This chain corresponds to this stoichiometry chain
                                if(alignment[0][2]/len(alignment[0][0])> threshold and
                                           current_number_of_chains[seq.id] < stoichiometry[seq.id]):
                                    model.child_list.append(chain2)
                                    current_number_of_chains[seq.id] += 1
                                    print('Stoichiometry 2')
                                    count = 0
                                    for chain in current_number_of_chains:
                                        if current_number_of_chains[chain] == stoichiometry[chain]:
                                            count += 1
                                    if count == len(stoichiometry):
                                        print('Stoichiometry fullfilled!')
                                        return full_structure
                                    print(current_number_of_chains)
                        else:
                            model.child_list.append(chain2)
                    failed = []
                    done.append(structure.id)
                    break
        else:
            failed.append(structure.id)
            comparisons[structure.id] = comparisons.get(structure.id, 0) +1
            # Pdbs clash
        current += 1
    return full_structure