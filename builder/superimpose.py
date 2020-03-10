#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 18:22:06 2020

@author: fabian
"""
import copy

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
    
    def __init__(self, structure_1, structure_2, fasta_1, fasta_2):
        self.structure_A, self.structure_B = structure_1, structure_2
        self.superimposer = PDB.Superimposer()
        self.fasta_1 = list(fasta_1)
        self.fasta_2 = list(fasta_2)
        self._rms = None
    
    def get_best_alignment(self):
        '''
            Obtains the best alignment between 2 fasta files sequences.
            
            The best alignment shows which chain is common between those files. 
            To know which two sequences are more similar, we make an alignment to determine if they are homologus

        '''
        aligned = []
        alignments = []
        for seq1 in self.fasta_1:
            for seq2 in self.fasta_2:
                #if (seq1 == seq2): continue
                #if (seq2, seq1) in aligned: continue
                aligned.append((seq1, seq2))
                aligned_pair = pairwise2.align.globalxx(seq1, seq2)
                aligned_pair.sort()
                max_score = max(aligned_pair, key=lambda x: x[2])
                alignments.append(max_score)
                #(alignments.append(pair) for pair in aligned_pair)
        return max(alignments, key=lambda x: x[2])
    
    def superimpose(self, fixed_chain_id):
        chain_A = next(filter(lambda x: x.id == fixed_chain_id, self.structure_A.get_chains()))
        chain_B = next(filter(lambda x: x.id == fixed_chain_id, self.structure_B.get_chains()))
        chain_atoms = list(chain_A.get_atoms())
        superimposed_atoms = list(chain_B.get_atoms())
        self.superimposer.set_atoms(chain_atoms, superimposed_atoms)
        self.superimposer.apply(list(self.structure_B.get_atoms()))
        self._rms = self.superimposer.rms
    
    def save_model(self, path, chain_duplicated):
        final_structure = copy.copy(self.structure_A)        
        for chain in self.structure_B.get_chains():
            if(chain.id == chain_duplicated): continue
            final_structure[0].add(chain)
        self.io.set_structure(final_structure)
        self.io.save(path)
            
    
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
    
    def model(self, *chains):
        '''
            Takes a list of structures an superimposes them in order to make a bigger complex
        '''
        if (len(chains) > 1):
            self.initial_structure = self.model(*chains[1:])
        structure = chains[0]
        model_atoms = self.initial_structure.get_atoms()
        structure_atoms = structure.get_atoms()
        self.superimposer.set_atoms(model_atoms, structure_atoms)
        self.superimposer.apply(structure_atoms)
        self._rms = self.superimposer.rms
        return self.initial_structure
    