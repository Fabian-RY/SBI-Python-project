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