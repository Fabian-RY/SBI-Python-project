#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 18:22:06 2020

@author: fabian

This module contains a class that makes the alignments and superimpositions of 
the chains in order to assembly them. However, it does not join the chains, as some 
extra conditions may be required to do so. For example, add only necessary chain
for a desire stoichiometry
"""

import copy

from Bio import PDB
from Bio import pairwise2

class Ensemble():
    
    io = PDB.PDBIO()
    
    '''
        This class takes structures of proteins (Biopyton PDBs) and superimposes them
        in order to make one bigger structure with all the chains given
    '''
    
    def __init__(self, structure_1, structure_2):
        '''
            
        '''
        self.structure_A, self.structure_B = structure_1, structure_2
        self.superimposer = PDB.Superimposer()
        self._rms = None
    
    def get_best_alignment(self, fasta_seqs):
        '''
            Obtains the best alignment between 2 fasta files sequences.
            
            The best alignment shows which chain is common between those files. 
            To know which two sequences are more similar, we make an alignment to determine if they are homologus
            
            Returns:
                - The alignment with best score, which is supposed to be the repeated sequence

        '''
        
        aligned = []
        alignments = []
        # Each possible pair of chains is aligned and all of them are returned
        for chain_A in fasta_seqs[0]:
            seq1 = fasta_seqs[0][chain_A]
            for chain_B in fasta_seqs[1]:
                seq2 = fasta_seqs[1][chain_B]
                #if (seq1 == seq2): continue We make all posibilities so we don't skip e
                if (chain_B, chain_A) in aligned: continue
                aligned.append((chain_A, chain_B))
                # The alignment algorithms does several of them, so we take only the best one
                aligned_pair = pairwise2.align.globalxx(seq1, seq2)
                max_score = max(aligned_pair, key=lambda x: x[2])
                alignments.append(((chain_A, chain_B), max_score[2]/len(max_score[0])))
        # We return the chain ids of the best alignment and its score
        return alignments
    
    def superimpose(self, chain_A, chain_B):
        '''
            Superimpose two complexes given two chains. This chains should be homologous, 
            and they are used to calculate the matrices and rotations needed to be applied 
            to the rest of the chains
            
            input: The chain id to superimpose of structure A and the chain to superimpose of chain B
            returns a list with the atoms with their coordinates changed
        '''
        # Obtaining the chains
        chain_A = next(filter(lambda x: x.id == chain_A, self.structure_A.get_chains()))
        chain_B = next(filter(lambda x: x.id == chain_B, self.structure_B.get_chains()))
        # Obtaining the chains atoms
        moving_chains = copy.copy(self.structure_B)
        model = next(moving_chains.get_models())
        moving_chains = model.get_chains()
        superimposed_atoms = list(chain_B.get_atoms())
        chain_atoms = list(chain_A.get_atoms())
        # If both sequences have different number of atoms, get only the common part
        # This can be improved using alignments
        if (len(list(chain_A.get_atoms())) != len(list(chain_B.get_atoms()))):
            length = min(len(list(chain_A)), len(list(chain_B)))
            chain_atoms = chain_atoms[:length]
            superimposed_atoms = superimposed_atoms[:length]
        # Doing the superimpositions and returning a list with the atoms
        # With the changed positions
        self.superimposer.set_atoms(chain_atoms, superimposed_atoms)
        chains_changed = list()
        for chain in moving_chains:
            if chain.id == chain_B.id: continue # We skip the homologous one, it's present in both structures
            atoms = chain.get_atoms()
            self.superimposer.apply(atoms)
            chains_changed.append(chain)
        self._rms = self.superimposer.rms
        return chains_changed
    
    @property
    def chains(self):
        return self.structure_A, self.structure_B
       
    @property
    def rms(self):
        '''
            Returns the root-mean-square distance of the superimposed proteins.
            If superimposition wasn't done raises an error.
        '''
        if(self._rms is None):
            raise RuntimeError('Cannot determine rms unless superimposition is done')
        else:
            return self._rms
    