#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 10:32:34 2020

@author: fabian
"""

from distutils.core import setup

setup(name='pro-mod-FC',
      version='0.1',
      description='A macromolecular complex builder',
      author='Fabián & Claudio',
      packages=['builder'],
      scripts=['main.py', 'scripts/pairpdbs.py', 'scripts/pdbsplit.py'],
      install_requieres = ['Biopython'],
      url='https://github.com/Fabian-RY/SBI-Python-project')