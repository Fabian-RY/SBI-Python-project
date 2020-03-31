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
      author_email='fabianry97@gmail.com',
      packages=['builder'],
      scripts=['main.py', 'main-tk.py','scripts/pairpdbs.py', 'scripts/pdbsplit.py'],
      install_requieres = ['Biopython'],
      url='https://github.com/Fabian-RY/SBI-Python-project',
      include_package_data=True)