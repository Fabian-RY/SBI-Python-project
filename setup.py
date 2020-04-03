#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 10:32:34 2020

@author: fabian
"""

from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'theory.md'), encoding='utf-8') as f:
    long_description_in = f.read()

setup(name='promod',
      version='0.1.7',
      description='A macromolecular complex builder',
      long_description=long_description_in,
      long_description_content_type='text/markdown',
      author='Fabián Robledo, Claudio Díaz',
      author_email='fabian.robledo01@estudiant.upf.edu',
      license='MIT',
      keywords='bioinformatics structural_bioinformatics protein complex builder modeler',
      packages=['builder'],
      scripts=['promod', 'promod-tk','scripts/pairpdbs.py', 'scripts/pdbsplit.py'],
      install_requieres = ['Biopython'],
      url='https://github.com/Fabian-RY/SBI-Python-project',
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
      ],
      python_requires='>=3',
      include_package_data=True)