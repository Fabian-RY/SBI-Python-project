# SBI-Python-project

## Introduction

This programs aims to construct a protein complex given several pairs of chains that interact. This chains can be proteins, 
single DNA/RNA strands or double DNA/RNA strands.

## Instalation

To install this program open a terminal and navigate to the folder where the program was downloaded,
and run 'pip3 install .' to install the program. You can also use 'python3 setup.py install', however
using the first comand you'll be able to uninstall it easily using 'pip3 uninstall'

## Examples of use

The executable script is main.py. By executing 'main.py -h' you'll se how to execute it:

The necessary arguments are:

- The folder where the interacting pdbs are stored

- The folder where the final model will be stored. A dot ('.') can be used to select the current directory

- The fasta file with the sequences of the chains in the pdbs

*Note:* a fasta file whose sequences don't correspond to the fasta file will throw an error and stop working

Thus, at least those arguments must be fullfilled. As an example:

main.py -i examples/3e0d/pairs -o . -f 3e0d.fa

However, other arguments are elegible, to change some internal behavior of the program. Those are:

- The distance (-d), in Armstrongs, to consider that two chains clash. That means, if there are atoms nearer than the indicated
distance, it will be considered an unviable protein and the chain won't be added. A very big value can produce very incomplete proteins, while
a very low one will produce proteins with very high energies and possible unviable.

- The threshold (-t), a value to determine if two sequences are same to the stoichiometry and the superposing. A value higher than threshold consider sequences as homologous.
Default value is 0.95 (95% of identity). Lesser values can be given to increase tolerance, but being too low can introduce incorrect chains.

## Dependencies

This program relies on *Biopython* packages in order to work properly. 
Thus, you need this package in order to make the program work

If you want to optimize the final model, you also need to install modeller, which can be found [here](https://salilab.org/modeller/).
There you can download, get a license and install MODELLER

## Examples

There are a few examples in the *examples* folder. With a few chains, protein and DNA strands, and big protein complexes.