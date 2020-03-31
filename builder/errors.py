# -*- coding: utf-8 -*-

class seq_not_in_pdb(Exception):
    
    def __init__(self, seq, pdb):
        self.seq = seq
        self.pdb = pdb
        pass
    
    def __str__(self):
        pass
    
class non_valid_input(Exception):
    
    def __init__(self):
        pass
    
class stoichiometry_error(Exception):
    
    def __init__(self, line):
        self.line = line
        pass
    
    def __str__(self):
        return 'Not a valid stoichiometry line %s' % self.line
    
    @property
    def error(self):
        return self.line
    
class chain_in_stoic_not_in_fasta(Exception):
    
    def __init__(self, sequence):
        self._seq = sequence
        
    def __str__(self):
        return 'Seq %s not found' % self.seq
    
    @property
    def seq(self):
        return self._seq
    
class PDB_disagrees_fasta(Exception):
    
    def __init__(self, pdb_id):
        self.id = pdb_id
        
    def __str__(self):
        return 'PDB sequences does not correspond to any fasta sequence'
    
    @property
    def pdb(self):
        return self.id