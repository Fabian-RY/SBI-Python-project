# Promod - A protein complex builder

![version](https://img.shields.io/badge/version-1.0-blue) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![](https://img.shields.io/badge/platforms-linux%2C%20windows-brightgreen.svg)



## Introduction

Promod is a Python tool whose goal is to model macrocomplexes of  molecules starting from the interacting pairs of chains that will form  the complex. It is compatible with protein chains, single and double DNA strands and RNA strands. It has several parameters available to adjust  according to the users needs.

In essence, the builder approach is similar to a genomic assembler:  it seeks for similar chains in the different PDB files given as input  and overlaps them. The builder places pair by pair all the possible  chains into a single model.

## Features

1. Building of macromolecular complexes from basic input data (pairs of interactions between molecules).
2. Optimization of the final model using MODELLER.
3. Graphical user interface (GUI), with the same functionalities as the command line interface (CLI).

## Requirements

This software is written in the Python3 language programming. You need to [install the Python3](https://wiki.python.org/moin/BeginnersGuide/Download) interpreter on your computer.

The [Biopython library](https://biopython.org/) is also needed for the correct performance of this software. Independently of the operating system you are using, you can install Biopython using `pip`:

```
pip3 install biopython
```



## Installation

The recommended way to install Promod is by using `pip`, exactly with the same command as we indicated with Biopython above.

```
pip3 install promod
```

You can install Promod easily by downloading it from the Github repository and running the next command in the downloaded folder. That is the recommended installation procedure.

```
pip3 install .
```

Alternatively, you can also run the `setup` installation script, as follows:

```
python3 setup.py install
```

You can **uninstall** this software by using the following command, if it was installed using `pip3`:

```
pip3 uninstall promod
```

However, if it was installed via `setup.py`, all files must be deleted manually.

## Execution

The basic usage of Promod, providing only the mandatory arguments, would be as follows:

```
promod -i [INPUT_FOLDER] -o [OUTPUT_FOLDER] -f [FASTA_FILE]
```

Additional arguments can also be indicated. For further information, read the manual or use the package help:

```
promod -h
```

To run the Graphical User Interface, you only need to call the corresponding file:

```
promod-tk
```



## Documentation

- [Background and scientific explanation](theory.md)

- [Tutorial with examples on how to use the program](tutorial.md)

- [Analysis of examples shown in the tutorial](analysis.md)

A PDF file with all the contents is also [available](manual.pdf).



## About

This program is a project for the Structural  Bioinformatics and Introduction to Python subjects of the [MSc in  Bioinformatics for Health Sciences](https://www.upf.edu/web/bioinformatics) at [Universitat Pompeu Fabra  (Barcelona)](www.upf.edu).

#### Fabián Robledo Yagüe and Claudio Díaz García

![](https://www.upf.edu/documents/7283915/7376028/marca.png/906677bd-c708-255e-92a3-99524273d872?t=1484925837278)

