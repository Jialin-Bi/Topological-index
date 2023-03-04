# Topological-index
Multiscale topological indices

Multiscale topological indices for the quantitative prediction of SARS CoV-2 binding affinity change upon mutations

This manual is for the code implementation of the paper "Multiscale topological indices for the quantitative prediction of SARS CoV-2 binding affinity change upon mutations".
Preparation

Code Requirements Platform: Python>=3.6
Details about each step

Step 1: read mutation coordinate from PDB

Examples of mutated coordinates are in the “coodinates” folder

Step 2: Simplicial complex representation and Topological indices generation

For each protein, the coordinate or distance matrix is used to construct the simplicial complexes to generate the adjacency matrixes.

complex_index_A0_fil_more.py is used to compute the topological indices of a 0-dimensional adjacency matrix for the mutant type.
complex_index_A0_fil_more_origin.py is used to compute the topological indices of a 0-dimensional adjacency matrix for the wild type.
complex_index_A1_fil_more.py is used to compute topological indices of a 1-dimensional adjacency matrix for the mutant type.
complex_index_A1_fil_more_origin.py is used to compute topological indices of a 1-dimensional adjacency matrix (including upper, lower, and general) for the wild type.

Step 3: Correlation between the topological index and binding affinity change

pearson_result_add.py is used for the correlation of all the datasets.
pearson_result_add_group.py is used for the correlation of different group type
