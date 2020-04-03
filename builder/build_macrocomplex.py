# -*- coding: utf-8 -*-

"""
Created on Tue Feb 25 19:39:13 2020

@author: fabian

@version: 1.0

"""

from .superimpose import Ensemble
from . import errors

from Bio import PDB
from Bio import pairwise2

def _chain_id():
    '''
        The differents id the chains will have inside the model. Now supports a model up to 60 different ids:
        25 for uppercase letters, 25 for lowercase letters and 10 more from numbers 0 to 9.
    '''
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
    '''
        Determines which of the sequences of the fastas correspond to each of the pdb chains
        Two similar chains can refer to the same fasta sequence. 
        
        It will raise an exception if a schain in the pdb does not have an homologous in the fasta. 
    '''
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
                raise errors.PDB_disagrees_fasta(struct.id)
                pass
    return chain_fasta

def build_complex(threshold, distance, stoichiometry, sequences, structures, verbose, initial):
    '''
        This is the core function of the program. 
        
        Takes a list of structures and tries to join them in a single model, taking into account
        diferent parameters, such as an intial model, a threshold for the pairwise alignments and 
        the distance to consider that two chains are in the same position and it must be discarded 
        (different atoms cannot occupy the same space)
    '''
    # Initiate empty structure and check some errors
    seqs = get_fastas_from_structs(structures, sequences)
    
    # Initializing an empty structure to save the model
    full_structure = PDB.Structure.Structure('full')
    full_structure.add(PDB.Model.Model('model'))
    # We will need to assign ids. For now is limited to up to 64 different chains, but more single characters can be added
    ids = _chain_id()
    iters = 0
    pair_num = len(structures)
    # Prepare a dict that saves how many chains of every one in the stoichiometry has the model in construction
    if(stoichiometry):
        current_number_of_chains = {chain_id:0 for chain_id in set(stoichiometry)}
    # We need at least 2 pdbs with an interaction. Can be the same one repeated. This is done to avoid trouble
    if(pair_num < 2):
        raise ValueError('Needed at least 2 pdbs to superpose')
    failed = []
    current = 0
    done = []
    
    
    # If an initial structure was given, introduce it into the model
    if(initial):
        if verbose: print('Initializing complex')
        # Obtaining the sequences of the chains
        init_seqs = get_fastas_from_structs([initial], sequences)
        seqs[full_structure.id] = {}
        # For each chain, add it to the model with sequential ids: A->B-> etc
        for idx, chain in enumerate(initial.get_chains()):
            chain_id = next(ids)
            seqs[full_structure.id][chain_id] = init_seqs[initial.id][chain.id]
            chain_new = PDB.Chain.Chain(chain_id)
            chain_new.child_list = list(chain.get_residues())
            model = next(full_structure.get_models())
            # If stoichiometry is selected, then count which sequences will be added
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
            model.child_list.append(chain_new)
        done.append(initial.id)
    if verbose: print('Starting to build')
    # The main loop of the function
    while True:
        if verbose: 
            print('Loop #%i' % current)
            print('Current number of chains: %i' % len(list(full_structure.get_chains())))
        # All the structures where added correctly (should only be for non stoichiometric uses)
        if(not stoichiometry and len(done) == len(structures)): 
            break
        structure = structures[current % len(structures)]
        if(structure.id in done): 
            current += 1
            continue # If is already in, continue to the next one
        count = list(failed.count(element) for element in failed)
        # If failed has an element twice, that means that no more available chains could be added
        # Therefor, it will start an endless loop, as the structure remain equal no matter which 
        # of the left structures is trying to be added. We stop it here
        if 3 in count: # We are repeating structures
            if verbose: print('Some pdbs could not be joined.\nAnd endless loop started and no more chains could be added\nFinishing the model in the current state')
            break
        if len(list(full_structure.get_chains())) == 0:
            if verbose: print('Initializing complex')
            seqs[full_structure.id] = {}
            for idx, chain in enumerate(structure.get_chains()):
                chain_id = next(ids) # Selecting the next id available
                seqs[full_structure.id][chain_id] = seqs[structure.id][chain.id] # Assigning the new chains its sequence
                # Empty new chain with the new id
                chain_new = PDB.Chain.Chain(chain_id)
                # Adding the residues to the new chain and adding the chain to the new model
                chain_new.child_list = list(chain.get_residues())
                model = next(full_structure.get_models())
                # Check stoichiometry if needed
                if(stoichiometry):
                    # Need the exact sequence to align
                    chain_sequence = _get_chain_sequence(chain_new)
                    # Align with every sequence in the fasta to know which sequence is 
                    for seq in sequences:
                        alignment = pairwise2.align.globalxx(chain_sequence, seq.seq, one_alignment_only=True) 
                        try:
                            current_number_of_chains[seq.id]
                        except KeyError:
                            # The sequence chain is not in the stoichiometry
                            raise errors.chain_in_stoic_not_in_fasta(seq.id)
                        # If homologous and does not surpass the stoichiometry, add it
                        if(alignment[0][2]/len(alignment[0][0])> threshold and
                            current_number_of_chains[seq.id] < stoichiometry[seq.id]):
                            model.child_list.append(chain_new)
                            current_number_of_chains[seq.id] += 1
                            count = 0
                            for chain in current_number_of_chains:
                                if current_number_of_chains[chain] == stoichiometry[chain]:
                                    count += 1
                                if count == len(stoichiometry):
                                    print('Stoichiometry fullfilled!')
                                    return full_structure
                model.child_list.append(chain_new)

            done.append(structure.id)
            current += 1
            continue
        # The Enseble class keeps code ordered and abstracts the alignment
        # and superimposition of the code
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
                if (clashes < 10 and pair.rms < 0.05):
                    # The chain can be added, so let's add it
                    for chain in atoms_of_chains:
                        chain_id = next(ids) # Chosing an id
                        chain2 = PDB.Chain.Chain(chain_id) # An empty chain which will add the residues to join
                        chain2.child_list += list(chain.get_residues())
                        model = next(full_structure.get_models())
                        if(stoichiometry): # Taking the stoichiometry into account
                            chain_sequence = _get_chain_sequence(chain)
                            for seq in sequences:
                                alignment = pairwise2.align.globalxx(chain_sequence, seq.seq, one_alignment_only=True) 
                                # This chain corresponds to this stoichiometry chain
                                if(alignment[0][2]/len(alignment[0][0])> threshold and
                                           current_number_of_chains[seq.id] < stoichiometry[seq.id]):
                                    model.child_list.append(chain2)
                                    current_number_of_chains[seq.id] += 1
                                    count = 0
                                    for chain in current_number_of_chains:
                                        if current_number_of_chains[chain] == stoichiometry[chain]:
                                            count += 1
                                    if verbose: print('Current number of chains: %i' % sum(current_number_of_chains.values()))
                                    if count == len(stoichiometry):
                                        if verbose: print('Stoichiometry fullfilled!')
                                        return full_structure
                        else:
                            model.child_list.append(chain2)
                            if verbose and stoichiometry: print('Current number of chains: %i' % sum(current_number_of_chains.values()))
                            done.append(structure.id)
                    failed = []
                    break
        else:
            failed.append(structure.id)
            # Pdbs clash
        current += 1
    return full_structure