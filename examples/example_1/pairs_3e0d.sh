#! /bin/bash

pdbsplit.py -i 3e0d.pdb -o chains/ -p 3e0d -f
pairpdbs.py chains/
cat chains/*.fa > 3e0d.fa