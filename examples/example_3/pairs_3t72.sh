#! /bin/bash

pdbsplit.py -i 3t72.pdb -o chains/ -p 3t72 -f
pairpdbs.py chains/


