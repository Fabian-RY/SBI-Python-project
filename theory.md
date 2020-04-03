# Background and scientific explanation

## Introduction

Proteins are the executive molecules in all organisms. They perform a wide variety of functions, from small compounds transport to signals transduction or immune responses. However, proteins do not usually work individually in cells. Most proteins interact with other proteins (protein-protein interactions, PPIs) or molecules (DNA, hormones, drugs, among others), which are essential for a proper development of their activities. That is the reason why collecting PPIs can lead to a better understanding of protein functions, biological pathways and mechanisms of disease. The identification of such PPIs has been a challenging task for years. Experimental methods provide irrefutable evidence to test interactions between proteins, but they are expensive and time-consuming. Nowadays, computational approaches are being developed to reduce as much as possible the necessity of experimental data. Some advances have been made in this field, but the range of successful organisms is short and general frameworks are lacking at the moment.

Different experimental approaches for the identification of individual PPIs are available. The most accurate methods that allow the determination of the exact atom coordinates in the complexes are X-ray crystallography, Nuclear Magnetic Resonance and Electron Microscopy. These techniques are difficult to perform, they need very specialised equipments and a strong data analysis procedure. There are also high-throughput methods based on biochemical properties, such as Tandem Affinity Purification, which can detect multiple PPIs at once, but a high rate of false positive results can be expected and the information is not very precise. Protein microarrays are also being used as a screening method that can be easily automated and parallelised.  Other traditional techniques, namely co-immunoprecipitation, yeast two-hybrid or pull-down, are still being used as well. All these approaches provide complementary views of PPIs and have their advantages and problems. The main one is the cost-effectiveness of the experiments, making in silico approximations more feasible and adequate for understanding PPIs at the atomic level.

Regarding the structural properties of the proteins, there are three main categories of methods for computational modelling of PPIs: homology or template-based modelling, *ab initio* or template-free modelling, and hybrid or integrative modelling. Homology modelling is based on the fact that the evolutionary information in both the sequences and the structures is important for PPI prediction. The latter are preferred, as most existing predictors use surface patch data, and only the residues in the interface have to be analysed, instead of the whole sequence of aminoacids. A pitfall of this method is that not all the 3D structures are available, although public databases such as PDB are unstoppably increasing in size. Template-free modelling needs more computational resources, as it explores all the possible orientation between the interacting molecules. The combination of prior knowledge about the individual structures of the components may help reduce the searching space, but still this approach is more challenging. This kind of techniques also require a careful evaluation of the results by various means and refinement of the candidate models using biological information. Finally, integrative modelling combines experimental data and bioinformatics developments to narrow the possible complexes and save computational resources and time.

In this work, our aim is to develop a software that builds protein macrocomplexes using as input pairs of protein-protein interactions. Interaction with other molecules, such as nucleic acids, is also considered. To do this, our program is based on homology modelling. Therefore, the homology of the different chains is analysed, tridimensional structures are superimposed and energy levels of the final models are considered to propose the best possible solution. Another option would have been template-free modelling, but the required computational resources and the steps of evaluation and refinement excluded this possibility, because of the limited means and knowledge that our team has got.

Here we propose Promod as a first step approach to protein-protein and protein-nucleic acid complex modelling. It is distributed as both a Python package and a standalone application to be run in UNIX-based systems. Several parameters can be tuned to adjust to the user's needs. The final result is a single model with the best option according to the algorithm that is explained below. 



## The Promod method

The algorithm used by our program can be divided in three main steps. The first one is the analysis of the homology of the subunits. Then, the superimposition of the homologous structures is performed. Finally, energy levels in the models are taken into account to discard unlikely complexes.

#### Homology of the subunits

In order to build the model, the program begins with several files, each of them containing information about two interacting chains, whether peptidic, DNA or RNA. We can overlap similar proteins from two different PDB files to check if two chains in those files are alike. If they effectively are, they can be joined in the same model. This is known as protein superimposition, and the whole process can be compared to a DNA sequence assembly. The main objective is to recognise that two sequences are the same and put them together in the correct spot.

Two chains in different PDB files must be homologous to be considered as part of the same protein and overlap them. To check if two proteins are homologous or not, a pairwise alignment is performed. This alignment consists of calculating the similarity between two sequences, taking into account both mutations and gaps,  and then returning a score between 0 and 1. A score close to 1 means that both sequences are identical, therefore high score values are evidence of a large proportion of sequence identity. Only the homologous proteins will be overlapped and used to build the model later.

#### Superimposition of the 3D structure

Two homologous sequences will probably have similar structures. However, it is also possible that two proteins with lower identity score may have similar structures in the tridimensional space, as a result of convergent evolution. Therefore, when two proteins are significantly different in sequence but they have similar structures, the may have the same role in our model and it is desirable to consider them. For example, two homologous proteins from distant species whose alignment do not pass our homology threshold, but they still preserve the structure.

In these cases, to test if two proteins really have similar structures, we need a measurement, such as the root median square deviation (RMSD). First, to calculate this value, we have to superimpose these two proteins. That means placing the atoms of both proteins in the same coordinates and orientation, to check how well they overlap in the space. Similar structures will have atoms in almost the same positions, while different structures will be more distant. This similarity between the superimposed proteins can be used to calculate the RMSD, which is the average distance between the superimposed atoms in the chains. A value close to zero indicates a perfect fit of the structures. RMSD will increase when the differences between the protein structures increase.

Due to all this, RMSD must also be a filter to account for those structures that are different, even when we had considered them homologous in previous steps.

#### Analysis of the energy levels

Lastly, it is important to consider the final energy levels of the complex. A good model should have minimum energy, as functional, naturally occurring complexes are the ones with the lowest energy among the set of possible foldings. Situations such as two atoms very close to each other, aminoacids with hydrophobic residues located in the external part of the macromolecule and several other cases can increase the final energy of the complex, which should be then corrected.

After that, it is possible that the position of the atoms inside the model is not the most adequate. Sometimes, despite having evidence of interaction between two chains, collisions or clashes appear when they are joined in the structure. Taking this into account avoids impossible models, such as those with chains crossing with each other, which will be thermodynamically unstable and unlikely to happen.



## Features

The most striking features of Promod are:

1. Building of macromolecular complexes from basic input data (pairs of interactions between molecules).

2. Optimization of the final model using MODELLER.

3. Graphical user interface (GUI), with the same functionalities as the command line interface (CLI).

   

## Implementation

The Promod package is implemented in Python3. It strongly depends on the Biopython library, which is an installation requirement for the correct execution of this software. The GUI has been developed under the Tkinter framework, providing a user-friendly interface and, thus, avoiding the use of the command line. The implementation of the GUI and the CLI are independent, so each one can be executed on its own. This software works with well-established bioinformatics formats, such as FASTA and PDB files, so no atypical formatting of input data shall be conducted. Additionally, two scripts are provided to divide a given PDB file in separate pairs of chains and join them if there is an interaction between the molecules. Some examples are also included for the testing of the program. Promod is freely available from the following Github repository. URL: https://github.com/Fabian-RY/SBI-Python-project.

The formatted documentation of the package includes dependencies and installation instructions, examples of use and a full tutorial with sample code. Each function has got its own documentation that instructs on the particularities when importing them to be used independently. The users could learn all knowledge of the library by looking up the detailed documentation. The main functionalities are further explained in the tutorial section.



## Discussion

Promod is a basic tool that relies on information from both the sequence and the tridimensional structure of pairs of interacting molecules to build a final complex model. As we can see in the analysed examples, it is successful with small complexes and it can even handle nucleic acids interaction with proteins. The MODELLER final step adds an additional refinement of the proposed model, returning to the user a good quality structure. All this methodology is built in an easy-to-use package, well documented and with a GUI to save unpleasant command line work for the standard user.

However, the main objective of modelling software is to provide new knowledge without the need of huge amount of experimental data, such as the determination of all paired interactions, in this case. A limitation of our program is the necessity of input protein-protein or protein-nucleic acid interaction. This type of software is being broadly developed by means of various approaches, which are proving to be very successful. Both docking and homology modelling paradigms are being applied to achieve this goal.

One of the strategies being exploited is evidence combining methods. The core of this strategy is integrating evidence from multiple sources, including them in comprehensive databases for later integration. They are called gold standard databases and they contain information for both training and testing of the new methods. The annotation of paired protein interactions goes beyond structural features. For instance, evolutionary relationships, functional features, network topologies, sequence-based signatures, structure-based signatures, and text mining information is recorded in these databases. Finally, machine learning algorithms are fed with subsets of data and performance is measured to propose the better candidates, depending on target species, data sources, demand of accuracy and coverage. These evidence combining methods are performed repeatedly to find converging results with different input data and classifiers.

As we can see, template-based methods are leading the in silico modelling techniques. It is reasonable to believe that there are a limited number of possible interactions and that, once we have big enough databases and curated interaction catalogs, most of the PPIs should be easy to model with high confidence.     However, *ab initio* modelling has also yielded promising results, mainly from competitions such as CASP, CAPRI or CAMEO, where world-class groups bring their latest developments to test their performance with real problems. Of course, these solutions are computationally expensive and require a long time to return a final model. In addition, complexes with weak interactions where the conformational state changes upon binding are a big challenge for docking software.

A large number of information is nowadays available from high throughput experimental techniques and a lot of structural bioinformatics software has been developed to integrate these data. The best performance of in silico techniques has been achieved at the tertiary structure level of proteins. Nevertheless, quaternary structures are the ones responsible for the majority of biological functions and their knowledge is essential to disentangle protein interaction networks in both physiological and pathological scenarios. It is clear that hybrid approaches, joining atomic-level experimental structures, database information and the newest machine learning techniques, will give promising results in the near future. 



## Future perspectives

As we have already mentioned, there are diverse strategies for modelling PPIs. The development of refined algorithms to perform this task is far beyond our current knowledge in structural bioinformatics. Therefore, if we had time and resources to expand our project, we would focus on improving the user experience and the accessibility to existing data.

First, we would like to implement different ways of modelling. It would be useful for the user to choose whether to use a template-based or a template-free approach. We could also use existing Biopython packages to perform certain operations, but further research and testing would be needed to decide the better candidates. 

Next, the result of the modelling operations could return automatic reports about the process. Not only the final model would be reported, with the tridimensional structure and their characteristics, but also information about the reasons why our software has discarded certain possible models. Therefore, direct visualization of these other options and energy plots for user reference should be displayed. A good option for the format of this report would be a Jupyter notebook, so the user can interactively check all the available resources.

Another functionality that could be feasible to achieve is the automatic download of structures from public databases, such as PDB. Both sequences and structures can be obtained using Biopython modules and some keyword and ID search should be easy to implement, as well. These would be the input data for the additional scripts that build the pairs of interacting molecules to run the core Promod builder.

Finally, it would be ideal to develop integration options of biological information to help build better complexes, in line with trending research in the field. We are not aware of Python packages that would allow to directly implement this kind of work, but maybe external tools are capable of doing it. However, a broad review of the latest literature and a challenging programming development would be needed to know the most promising approximations.



## Conclusion

Although computational approaches for building macromolecular complexes are far to be perfect at the moment, a lot of effort is being made to develop strategies to overcome the pitfalls in this field. We have provided a software that, despite not using any novel strategy, achieves notable results when using already defined structures. Further steps using machine learning approaches and a broader range of training data would improve the performance and confidence of this type of programs. Finally, the integration of different types of biological data will also be a key step in the progress of in silico modelling of protein-protein interactions. The achievement of better methodologies will definitely have an impact in applied fields, such as drug discovery, biomedical research or food industry.

## Bibliography

Chang, J., Zhou, Y., Qamar, M. T. U., *et al.* Prediction of protein–protein interactions by evidence combining methods. *International Journal of Molecular Sciences* **17**, 1946 (2016). doi: 10.3390/ijms17111946

Ding, Z., Kihara, D. Computational identification of protein-protein interactions in model plant proteomes. *Scientific Reports* **9**, 8740 (2019). doi: 10.1038/s41598-019-45072-8

Hayes, S., Malacrida, B. , Kiely, M., Kiely, P. A. Studying protein–protein interactions: progress, pitfalls and solutions. *Biochemical Society Transactions* **44**, 994-1004. doi: 10.1042/BST20160092

Keskin, O., Tuncbag, N., Gursoy, A. Predicting protein−protein interactions from the molecular to the proteome level. *Chemical Reviews* **116**, 4884-4909 (2016). doi: 10.1021/acs.chemrev.5b00683

Liu, S., Liu, C., Deng, L. Machine learning approaches for protein–protein interaction hot spot prediction: progress and comparative assessment. *Molecules* **23**, 2535 (2018). doi: 10.3390/molecules23102535

Nealon, J. O., Philomina, L. S., McGuffin, L. J. Predictive and experimental approaches for elucidating protein–protein interactions and quaternary structures. *International Journal of Molecular Sciences* **18**, 2623 (2017). doi: 10.3390/ijms18122623

Sarkar, S., Gulati, K., Kairamkonda, M., *et al.* Elucidating protein-protein interactions through computational approaches and designing small molecule inhibitors against them for various diseases. *Current Topics in Medicinal Chemistry* **18**, 1-18 (2018). doi: 10.2174/1568026618666181025114903