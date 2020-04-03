# Manual

Welcome to Promod's. Here you will find all the information you need to run this software. This program is a project for the Structural Bioinformatics subject and the Introduction to Python project of the Msc in Bioinformatics of Universitat Pompeu Fabra.

## What is Promod?

Promod is a python tool whose goal is to form macrocomplexes of molecules starting from the interacting pairs of chains that will form the complex. It is compatible with protein chains, single and double DNA strands and RNA strands. It has several parameters available to play with according to the users needs.

In essence, the builder approach is similar to a genomic assembler: seeks for similar chains in the different PDB files given as an input and overlaps them. After that, joins the chains in a single model. Pair by pair the builder puts all the possible chains inside a single model

## Scientific Background

### Introduction 

Obtaining the tridimensional shape of a protein complex is a difficult and expensive process to obtain experimentally. Adding that there are thousands of different proteins in a single organism, with different splicings and mutations, which can vary even further in different spices, it's barely impossible to determine them all experimentally. Thus, finding computational ways to obtain the structure of a protein is advisable to know the structure of variations of a protein.

Here we present a computational tool to ensemble proteins and nucleic acids knowing their interactions. Currently, there are estimations that there are at least two hundred thousand different interaction between different chains of proteins that form every possible interaction in which a protein is present. Knowing those interactions, most of the complexes can be built easily knowing the stoichiometry of the protein.

### Homology of the subunits

In order to build the model, we start with several files: each file containis two chains (whether aminoacidical, DNA or RNA) representing an interaction between those chains. By looking to similar proteins in two different pdbs, we can overlap them in the space to check if two chains of those pdbs are similar. After that, we precisely use that similiarity to join them in a single model. This is known as protein superimposition, and the full process of joining is similar to a DNA sequence assembly: we have to recognise that two sequences are highly similar, and join them in the correct spot.

In order to consider two chains in diferent pdbs the same protein and overlap them, they must be homologous. To check whether two proteins are or not are homologous, a pairwise alignment is made. This is a calculation of how similar two sequences are, that takes both mutations and gaps into account to give a numeric score between 0 and 1. A score of 1 means that both sequences are identical, while high values indicate how similar the sequences are. So, to consider them homologous, a high score must be achieved between two proteins. Only the homologous proteins will be overlaped and used to build the model.

### Superimposition of the 3D structure

Two homologous sequences will probably have similar structures. However, it's also possible that two proteins with lower homology score can have similar structures in space, as result of convergent evolution. Therefore, when two sequences are quite different but they may have similar structures, the may fullfill the same rol in our model and it's advisable to consider them: for example, two homologous proteins from distant species whose alignment doesn't pass our homologous threshold, but they still conserve the structure.

To consider them in our model, we need another measure to check if they really have similar structures, the root median square deviation (RMSD now onwars). To check this value, first we need to superimpose those two proteins; that's putting the atoms of both proteins in the same place with the same orientation to check how good or bad they overlape in space; similar structures will have atoms in nearly the same position, while different structures will differ. The average distance (measured in argmstrongs) between the superimposed atoms is the RMSD. 

This way, RMSD must also be a filter we must also filter those structures that are differents even if we consider them homologous, as they can have very different 3D structures.

### Energy levels

Lastly, it's important to consider the final energy levels of the model. A good model of a complex should have minimum energy, as a functional, valid complex usually has the lowest energy among the possible foldings of the complex. Two atoms that are too close, hidrophobic residues in the external part of the model and several other considerations can increase the final energy of the complex, which must be corrected in a good model.

Two atoms too close can also indicate that two proteins are not correctly joined. For example, if several atoms of a chain might be too close of another chain (or even *inside*), those proteins should not be joined, even if they interact. Taking this into account avoids  impossible proteins, such as those which are inside other chains but breaking several laws of thermodynamics and physics. 

### Input

The input is a folder with two or more files: each file contains two proteins, a protein and DNA/RNA strand or 2 single DNA strands (which may or may not form a double strand), which represent an interaction between those chains. Two structures are considered to be interacting if the minimum distance between them is in a range of a few Argmstrongs (1-10) . However, in lesser distances, forces between atoms are strong enough to produce changes between the chains and force them to adopt a different conformation, an thus, a realistic model must keep the residues at an adequate distance to consider it correct.

### Output

The output is mainly one pdb file with all the possible chains joined in the selected folder. However, if the optimized option was selected, several pdb files will be created in the same directory. These are different approaches made by modeller to optimize the energies of the model and the files it used. The model will be saved as final_model.pdb

If the model is required with minimum energy, and thus the -optimize flag is used, then several files will be saved. Modeller will save some mid-step pdbs as 'final-model.DXXXXX.pdb' and the fully optimized will be saved as 'optimized.pdb'

