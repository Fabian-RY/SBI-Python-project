#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 18:22:06 2020

@author: fabian
"""

import copy
import itertools

from Bio import PDB
from Bio import Align
from Bio import SeqIO
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
        for chain_A in fasta_seqs[0]:
            seq1 = fasta_seqs[0][chain_A]
            for chain_B in fasta_seqs[1]:
                seq2 = fasta_seqs[1][chain_B]
                #if (seq1 == seq2): continue We make all posibilities so we don't skip e
                if (chain_B, chain_A) in aligned: continue
                aligned.append((chain_A, chain_B))
                # print(self.fasta_1[seq1], self.fasta_2[seq2])
                aligned_pair = pairwise2.align.globalxx(seq1, seq2)
                max_score = max(aligned_pair, key=lambda x: x[2])
                alignments.append(((chain_A, chain_B), max_score[2]/len(max_score[0])))
        return max(alignments, key=lambda x: x[1])
    
    def superimpose(self, chain_A, chain_B):
        '''
            Superimpose two complexes given 
        '''
        print(list(self.structure_A.get_chains()), list(self.structure_B.get_chains()))
        chain_A = next(filter(lambda x: x.id == chain_A, self.structure_A.get_chains()))
        chain_B = next(filter(lambda x: x.id == chain_B, self.structure_B.get_chains()))
        moving_chains = copy.copy(self.structure_B)
        model = next(moving_chains.get_models())
        moving_chains = model.get_chains()
        superimposed_atoms = list(chain_B.get_atoms())
        chain_atoms = list(chain_A.get_atoms())
        print(self.structure_A.id, self.structure_B.id)
        print(len(chain_atoms), len(superimposed_atoms))
        self.superimposer.set_atoms(chain_atoms, superimposed_atoms)
        chains_changed = list()
        for chain in moving_chains:
            if chain.id == chain_B.id: continue
            atoms = chain.get_atoms()
            self.superimposer.apply(atoms)
            chains_changed.append(chain)
        self._rms = self.superimposer.rms
        return chains_changed
    
    @property
    def chains(self):
        return self.structure_A, self.structure_B
    
    @property
    def fastas(self):
        return self.fasta_1, self.fasta_2
    
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
    