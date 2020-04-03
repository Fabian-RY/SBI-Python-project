# -*- coding: utf-8 -*-

'''
    Different errors that are recognised and the program must stop are written here
'''

class seq_not_in_pdb(Exception):
    '''
        The sequences are not in pdb and thus the protein complex cannot be build
    '''
    def __init__(self, seq, pdb):
        self.seq = seq
        self.pdb = pdb
        pass
    
    def __str__(self):
        pass
    
class non_valid_input(Exception):
    '''
        The input used is not correct and cannot be used
    '''
    def __init__(self):
        pass
    
    def __str__(self):
        return 'Invalid input. Check it twice'
    
class stoichiometry_error(Exception):
    '''
        The stoichiometry file is not correct and thus can not be used
    '''
    def __init__(self, line):
        self.line = line
        pass
    
    def __str__(self):
        return 'Not a valid stoichiometry line %s' % self.line
    
    @property
    def error(self):
        '''
            The stoichiometry line that is not correctly formated
        '''
        return self.line
    
class chain_in_stoic_not_in_fasta(Exception):
    '''
        A sequence asked in the stoichiometry is not available in the fasta file
    '''
    def __init__(self, sequence):
        self._seq = sequence
        
    def __str__(self):
        return 'Seq %s not found' % self.seq
    
    @property
    def seq(self):
        return self._seq
    
class PDB_disagrees_fasta(Exception):
    '''
        A PDB sequence chain has not an homologous in the PDB and cannot continue 
    '''
    def __init__(self, pdb_id):
        self.id = pdb_id
        
    def __str__(self):
        return 'PDB sequences does not correspond to any fasta sequence'
    
    @property
    def pdb(self):
        return self.id