#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Mar 19 10:01:26 2020

@author: fabian

A simple function that uses modeller to improve energie

"""

# Example for: conjugate_gradients(), molecular_dynamics(), model.switch_trace()

# This will optimize stereochemistry of a given model, including
# non-bonded contacts.


from modeller import environ, selection
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions

import os

def optimize(pdb, pdb_path):
    print(1, pdb_path)
    # Environ data
    env = environ(0)
    env.io.atom_files_directory = ['../atom_files']
    env.edat.dynamic_sphere = True
    
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    
    code = pdb.split('.')[0] 
    mdl = complete_pdb(env, pdb)
    mdl.write(file=code+'.ini')
    
    # Select all atoms:
    atmsel = selection(mdl)
    
    # Generate the restraints:
    mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
    mdl.restraints.write(file=code+'.rsr')
    
    mpdf_prior = atmsel.energy()
    
    # Create optimizer objects and set defaults for all further optimizations
    cg = conjugate_gradients(output='REPORT')
    md = molecular_dynamics(output='REPORT')
    
    # Open a file to get basic stats on each optimization
    trcfil = open(code+'.D00000001', 'w')
    
    # Run CG on the all-atom selection; write stats every 5 steps
    cg.optimize(atmsel, max_iterations=20, actions=actions.trace(5, trcfil))
    # Run MD; write out a PDB structure (called '1fas.D9999xxxx.pdb') every
    # 10 steps during the run, and write stats every 10 steps
    md.optimize(atmsel, temperature=300, max_iterations=50,
                actions=[actions.write_structure(10, code+'.D9999%04d.pdb'),
                         actions.trace(10, trcfil)])
    # Finish off with some more CG, and write stats every 5 steps
    cg.optimize(atmsel, max_iterations=20,
                actions=[actions.trace(5, trcfil)])
    
    mpdf_after = atmsel.energy()
    
    mdl.write(file=os.path.join(pdb_path, 'optimized.pdb'))
    return (mpdf_prior, mpdf_after)