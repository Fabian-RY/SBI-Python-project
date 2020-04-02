#! /bin/bash

pdbsplit.py -i 6gmh.pdb -o chains/ -p 6gmh -f
pairpdbs.py chains/ 15
cat chains/*.fa > 6gmh.fa
