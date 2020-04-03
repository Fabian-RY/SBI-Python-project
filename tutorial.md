## Tutorial

### Installing dependencies

In order to make Promod work, there are some dependencies that must be installed before executing the software: Python3 and Biopython.

If you are using Windows, you can download Python3 from [its website](www.python.org). If you have got Linux, it is likely that you already have got it installed. 

Independently of the operating system you run, you can install Biopython using `pip`.

```
pip3 install biopython
```

### Installation

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

### Hands on: how to use Promod

Promod has a command line interface which explains briefly how to use it before execution, when using the `-h` flag.

```
promod -h
```

There are 3 mandatory commands, while the rest are optional. The mandatory ones are related to the input files and output folder, required for the program to work. They are `-i`, the input folder; `-o`, the output folder; and `-f`, the fasta file of the sequences. All of them will be covered later in this manual.

It also has a graphical interface which accepts the same parameters, but graphically with message boxes and graphical aid. It is an independent executable, so the CLI can be used for automation by itself. The command to start the GUI is simple:

```sh
promod-tk
```

#### Input

The input must be a folder with two or more files: each file containing two proteins, a protein and a DNA or RNA strand or 2 single DNA strands (that do not necessarily form a double strand), which represent an interaction between those chains. Two structures are considered to be interacting if the minimum distance between them is in a range of a few Angstroms (1-10 Å) . However, in shorter distances, the forces between atoms are strong enough to produce changes between the chains and promote a different conformation, and, thus, a realistic model must keep the residues at an adequate distance to consider it correct.

#### Output

The output is mainly one PDB file with all the possible chains joined in the selected folder. However, if the optimized option was selected, several PDB files will be created in the same directory. These are different approaches made by MODELLER to optimize the energies of the model and the files it used. The model will be saved as `final_model.pdb`.

If the model is required with minimum energy, and thus the `-optimize` flag is used, then several files will be saved. MODELLER will save some mid-step PDB files as `final-model.DXXXXX.pdb` and the fully optimized complex will be saved as `optimized.pdb`.

#### Parameters

##### Input folder

This is a mandatory argument and should be indicated with the `-i` or `--input-folder` tags. The input folder should contain, at least, 2 PDB files. This folder may contain other files or subfolders. However, the program will ignore other files and will not check subfolders to find more PDB files. If you wanted to include some other PDB, you should copy it in the folder before running the program.

##### Output folder

This is a mandatory argument and should be indicated with the `-o` or `--output-folder` tags. The output folder is the directory where the output files are written. At the moment, the computed model file is saved as `final_model.pdb` and, thus, any other file with the same filename will be overwritten.

##### Fasta file

This is a mandatory argument and should be indicated with the `-f` or `--fasta` tags. The fasta file contains the sequences of the chains and the identifiers for the stoichiometry. It is a critical file, as the sequences must be homologous to the chain in the PDB file. We can hold up to 5% of differences (default threshold value which can be modified using the `-t` parameter). However, if one of the chains (if stoichiometry is not provided) or one of the chains in the indicated stoichiometry differs largely, the program will exit. The reason for this behaviour is to avoid guessing possible chains that could randomly match, as that would produce unexpected results.

**Note:** The sequences in the fasta file should be of a similar length than their counterparts in the PDB file. For instance, the fasta files and the PDB files downloaded directly from [Protein Data Bank](https://www.rcsb.org/) may differ in the initial and ending residues, as they are difficult to model, so they are usually removed from PDB entries.

##### Distance

This parameter can be indicated with the `-d` or `--distance` tags. You can indicate the minimum distance to consider that two proteins do not clash. The model will discard interactions with too many atoms at lower distance than the value indicated in Angstroms. Being too restrictive may discard some meaningful interactions, but being too flexible can force chains to overlap and produce interactions that are not energetically favorable.

##### Threshold

This parameter can be indicated with the `-t` or `--threshold` tags. The threshold is the minimum score an alignment must reach in order to insert the chain in the model and determine its homologous protein. A very high value (0.9 or more) is usually recommended to ensure that the chains are correctly identified. However, for sequences with less identity, it might be needed to decrease this value.

##### Stoichiometry file

This parameter can be indicated with the `-s` or `--stoichiometry` tags. An additional file with the stoichiometry can be provided to make sure that the final protein will contain a specific number of chains. This file must be properly formatted, one value per line, as follows:

```
chain_id_in_fasta,number_of_columns
```

For example, this is the stoichiometry file for the example number 1:

```{pseudocode}
3e0d_A,2
3e0d_B,2
3e0d_C,2
```

##### Starting PDB

This parameter can be indicated with the `-start` or `--start` tags. Promod's result is highly dependent on the first PDB used to build the model. In fact, results can vary a lot when different starting PDB files are used. Therefore, we provide an optional argument to select the starting PDB file to account for this variation. This allows to use the same PDB file with the same parameters, to obtain the same model. By default, the PDB files are sorted alphabetically and the first one in the list is selected.

### Uninstalling

You can uninstall this software by using the following command, if it was installed using `pip3`:

```
pip3 uninstall promod
```

However, if it was installed via `setup.py`, all files must be deleted manually.
