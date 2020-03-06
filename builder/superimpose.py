#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 18:22:06 2020

@author: fabian
"""

from Bio import PDB

class Ensemble():
    '''
        This class takes structures of proteins (Biopyton PDBs) and superimposes them
        in order to make one bigger structure with all the chains given
    '''
    
    def __init__(self, initial_structure):
        self.initial_structure
        self.superimposer = PDB.Superimposer()
        self.rms = None
        
    @property
    def rms(self):
        '''
            Returns the root-mean-square distance of the superimposed proteins.
            If superimposition wasn't done raises an error.
        '''
        if(self.rms is None):
            raise RuntimeError('Cannot determine rms unless superimposition is done')
        else:
            return self.rms
        
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
        self.rms = self.superimposer.rms
        return self.initial_structure
    
if __name__ == '__main__':
    parser = PDB.PDBParser