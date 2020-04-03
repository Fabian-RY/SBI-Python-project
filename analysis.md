## Analysis of examples

### Example 1 (3e0d)

3e0d is a small complex formed by 2 protein chains and 2 double DNA strands. It is a subunit of the eubacterial DNA polimerase, but it is perfect for an initial testing of our project. We can test how the program works, the time it takes to run and if the result is similar to the original one. In this case, we have 6 different interactions, as each DNA strand is counted as a single one interacting with another one and bound to the protein chain by one of them. In fact, this complex can be thought as a dimer: two monomers formed by a protein and a double strand DNA, which interacts with some residues in their protein chains.

The first test to our program consists of joining the different PDB files into one single structure, without further information about the stoichiometry. Therefore, the program will try to build the protein with only the information provided by the PDB files and the fasta file. The command to build it is the next one:

```
promod -i examples/example_1/chains/pairs -o examples/example_1/results  \
        -f examples/example_1/3e0d.fa
```

We only need to provide the folder with the input PDB files (`-i`), the desired output folder (`-o`) and the fasta file with the sequences of the proteins (`-f`). The program will take the different PDB files on the folder (by alphabetical order, to make it reproducible) and the result is incomplete:

![3e0d structure obtained from: a) rcsb.com b) promod without any extra info c) promod with starting point and selected stoichiometry d) without stoichiometry control][3e0d_composed]

We can observe that there is a missing double strand DNA. So, we will take a more elaborate approach to complete the model. As there is a missing protein, we should force the program to include it. There are two different ways of indicating this. The easiest one is selecting the starting complex of our model (which contains those missing parts). Different starting PDB files can result in different models, so, it is important to select a suitable one. By default, PDB filenames are sorted alphabetically and the first one is chosen as the starting point. As we said before, this is an important choice that is evidenced when we compare not choosing a starting point (Figure 1.b), choosing an adequate file (Figure 1.c) or choosing another one (Figure 1.d), which in this case results in models with additional or less chains than expected.

The other available option is to indicate the stoichiometry of the complex. We will focus on selection the starting PDB file, so the running command would be:

```
promod -i examples/example_1/chains/pairs -o examples/example_1/results \
       -f examples/example_1/3e0d.fa -start examples/example_1/chains/pairs/3e0d_ZG.pdb
```

This is more similar to the original structure, very close to the initial model. However, this is a very simple protein, and bigger complexes may be harder to build. However, from this example, we have learnt that:

* The program can build a model without any indication. However, it might be incomplete or it might have additional chains.
* The starting PDB file is an important choice that can influence the final model. Thus, it is important to remember which starting point gave which result. The same starting PDB file will give the same output.


### Example 2 (6gmh)

6gmh is a very big complex formed by 24 different chains, with a double strand of DNA. There are several chains that interact between them, but there are not repeated sequences: each monomer is unique. 

![6gmh structure obtained from a) rcsb.com b) builded without assistance c) selecting starting point and stoichiometry][6gmh_composed]

Our first try is the following one:

```
promod -i examples/example_2/chains/pairs -o examples/example_2/results  \
        -f examples/example_2/6gmh.fa 
```

In this case, there are several chains that has not been added to the model, and some others are repeated. We should select a different starting point and limit the chains with the stoichiometry, so there is one copy of each chain, at most.

```
promod -i examples/example_2/chains/pairs -o examples/example_2/results  \
        -f examples/example_2/6gmh.fa -start examples/example_2/chains/pairs/6ghm_VA.pdb \
        -s examples/example_2/6gmh.stoic
```

The `6gmh.stoic` file contains the stoichiometry of the protein. Basically, it is a csv file, without header, in which each line follows the same structure: *<chain_name>,<number_of_chains>*. Each chain_name must be in the fasta as a sequence identifier.

This starting PDB file was not in the model in our first attempt, so now we force it to be included in the computation, and start building from there. This way, the result contains more chains in the model and it is  more similar to the original one.

### Example 3 (5nss)

![5nss structure obtained from rcsb.org vs composed one][5nss_composed]

Here we want to show how to use two parameters we have not shown yet: distance and threshold. The distance is a value that avoids the situation when two atoms collide, which is not allowed because two very close atoms would have high energy. So, if the separation between them is less than the indicated distance, it is considered invalid. However, to avoid some errors, up to 10 colliding atoms are allowed. The threshold makes reference to the identity percentage between the PDB structures and the fasta sequences. This allows certain flexibility in sequences with some mutations that may alter the conformation of the chain.

```
promod -i examples/example_3/chains/pairs -o examples/example_3/results  \
        -f examples/example_3/5nss.fa -d 2 -t 0.9 
        -start examples/example_3/chains/pairs/5nss_NG.pdb
```

In this case, the result is quite good, but an important issue is the speed of the program. As the number of interactions grow, the time it takes the program to complete the modelling also increases. However, there are a few considerations about that:

1. Adding a new chain that satisfies all the restrictions is the most computationally expensive situation. The number of attempts increases with the number of chains in the model. 

2. It takes less time to discard an interaction that has not an homologous molecule already in the model than adding a new chain. During the first steps of the program, until the model has grown enough, this is the most common case. The number of calculated alignments is relatively low as there are a few chains in the model, but this may become an issue at certain point.

3. Once the model has a considerable number of chains, and particularly if there are several copies of a monomer, discarding an homologous protein because it does not fit with the current parameters is a very time-consuming process. The program will try to introduce the chain using all the possible alignments whose scores are higher than the threshold. This is an important problem for proteins with a lot of copies of several monomers.

4. The stoichiometry must be carefully selected to reach the desired structure. It may happen that the program cannot construct a model with the desire stoichiometry, but it is also possible that the stoichiometry was not appropriate for such case.

### Example 4 (6om3)

In this last example, we will cover the software behaviour in relation to the number of input PDB files and the influence of stoichiometric restrictions.

![6om3 builded from: a) rcsb.com b) promod with stoichiometry][6om3_composed]

```
promod -i examples/example_4/chains/pairs -o examples/example_4/results  \
        -f examples/example_4/6om3.fa -d 0.3 -t 0.95 \
        -s examples/example_4/stoic.csv
```

We can see that there are some absent chains in the protein complex built with Promod, which can be added using custom values of distance or an appropriate starting PDB file. Nevertheless, we will focus on another topic that has not been commented yet.

If you provide a stoichiometry file but it does not exist, the program will abort. We could have changed it to a non-stoichiometry builder, but we considered a better approach to explicitly state this situation and let the user manage this type of error.

An important thing about the stoichiometry is that it is critical to make sure that all the model is built following the user preferences. For example, without any assistance, the modelling in this example is good. However, with a selected stoichiometry, the result has some missing chains. Here we want to stress that the user should carefully check the stoichiometry file to make sure the chains are correctly selected and the number of them is adequate to our purpose.

## Limitations

1. The sequences in the PDB files and the sequences in the fasta files must be similar. The program admits a certain degree of tolerance (which can be adjusted by the user using the `-t` parameter), but selecting sequences in the fasta file which are significantly different from the PDB structures will not allow the assembly of the complex. This is particularly difficult for PDB files in which protein tails are not properly modelled and, therefore, they are not present in the structure, but they are contained in the fasta. The pairwise alignment with the tail will yield a low score, resulting in an incorrectly failed homology identification or even returning a fake corresponding sequence in the fasta just by chance.

   Our recommendation is to prepare a fasta file as similar as possible to the structural sequence to avoid this kind of errors. We provide an additional script, `pdbsplit.py`, which can extract the sequences of the chains inside the pdb file, preventing this kind of issues up to certain degree. So, if you are investigating the effect of certain mutations in a target protein, it might be adequate to prepare the fasta file according to the PDB sequences.

2. The running time of the program is proportional to the number of PDB files selected and it is also affected by the order of the structures and the starting point. The most critical step is reading all the PDB sequences. However, as the number of interactions increase, more chains might be added to the complex. This is a trade-off between the number of interactions and computational speed.

3. Selecting different starting points can produce different models. It is uncertain to state which is the most adequate starting point, as this may vary according to the goal of the assembly.

4. The fasta files require to be specially prepared for this application. There should not be repeated sequences (with at least an alignment score equal or higher than selected threshold. This also applies to sequences with unknown aminoacids, which are usually marked as X in the fasta file and will cause a mismatch with the sequence in the PDB file. Those sequences will be usually missing in the final model.



[3e0d_pure]: examples/example_1/imgs/3e0d_original.png
[3e0d_fail]: examples/example_1/imgs/3e0d_fail.png
[3e0d_builded]: examples/example_1/imgs/3e0d_builded.png
[3e0d_start]: examples/example_1/imgs/3e0d_start.png
[3e0d_composed]: examples/example_1/imgs/3e0d_composed.png

[6gmh_composed]: examples/example_2/imgs/6gmh-composed.png
[6om3_composed]: examples/example_4/imgs/6om3-composed.png

[5nss_pure]: examples/example_3/imgs/5nss-original.png
[5nss_builded]: examples/example_3/imgs/5nss-builded.png
[5nss_composed]: examples/example_3/imgs/5nss-composed.png