# Multiscale topological indices for the quantitative prediction of SARS CoV-2 binding affinity change upon mutations

This manual is for the code implementation of the paper "Multiscale topological indices for the quantitative prediction of SARS CoV-2 binding affinity change upon mutations".


## Data

The benchmark dataset from S-protein and ACE2 complex 6M0J, which contains 3669 experimentally determined binding affinity changes upon single mutation at S-protein RBD. More specifically, all the mutations occur in the region 331-530 of S-protein RBD. The 3D structures of the mutant complex from 6M0J can be generated by the SCAP of Jackal. We also give the mutation domains and binding domains (within a cut-of distance 10.0 Angstrom) for each of the mutant complex by the VMD.

Relevant PDB data can be found in [6m0j-rbd1-pdb-file dataset](https://drive.google.com/drive/folders/1dUHg50WNLhfWOuAQj5Oa3HNMawTQuFeL?usp=sharing)

6M0J-RBD1-index.xlsx gives the related binding affinity changes values

name.txt gives the related names of PDB data.

## Preparation


```
Code Requirements Platform: Python>=3.6

Python Packages needed: math, numpy>=1.19.5, scipy>=1.6.2, GUDHI 3.7.1, networkX 2.5
```

## Details about each step

__Step 1__: Read mutation coordinate from PDB

Examples of mutated coordinates are in the “coordinates” folder
```
read_mutation_side_more.py is used to read atoms coordinates in the mutation domain of each mutant protein complex.

read_origin_side_more.py is used to read atoms coordinates in the ralted mutation domain of wild type protein complex.
```

__Step 2:__ Simplicial complex representation and Topological indices generation

For each protein, the coordinate or distance matrix is used to construct the simplicial complexes to generate the adjacency matrixes.

```
complex_index_A0_fil_more.py is used to compute the topological indices of a 0-dimensional adjacency matrix for the mutant type.

complex_index_A0_fil_more_origin.py is used to compute the topological indices of a 0-dimensional adjacency matrix for the wild type.

complex_index_A1_fil_more.py is used to compute topological indices of a 1-dimensional adjacency matrix for the mutant type.

complex_index_A1_fil_more_origin.py is used to compute topological indices of a 1-dimensional adjacency matrix (including upper, lower, and general) for the wild type.
```

__Step 3:__ Correlation between the topological index and binding affinity change

```
pearson_result_add.py is used for the correlation of all the datasets.

pearson_result_add_group.py is used for the correlation of different group type.
```
