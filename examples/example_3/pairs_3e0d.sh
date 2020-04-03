#! /bin/bash

pdbsplit.py -i 5nss.pdb -o chains/ -p 5nss -f
pairpdbs.py chains/ 20
cat chains/*.fa > 5nss.fa
