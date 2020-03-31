#! /bin/bash

pdbsplit.py -i 6om3.pdb -o chains/ -p 6om3 -f
pairpdbs.py chains/ 20
cat chains/*.fa > 6om3.fa
