# GAMBA
GAMBA is a web-based application (www.dutchconnectomelab.nl/GAMBA/) that can be used to test whether the gene expression profile(s) of the input gene(s) and neuroimaging-derived brain features show overlapped spatial patterns. Source code of the website can be found at XXXX.

For details, please see:

- Wei Y. et al., Statistical testing and annotation of gene transcriptomic-neuroimaging associations, in preparation.

- Wei, Y. et al. (2019), Genetic mapping and evolutionary analysis of human-expanded cognitive networks. Nat Commun. https://doi.org/10.1038/s41467-019-12764-8

# Processing
GAMBA is a web-application with a front-end interface based on pre-processed gene transcriptomic data and imaging data. Processing steps and statistical analyses are included in "./processing"

Briefly, given an input gene expression data matrix (region by gene) and phenotypic data matrix (region by phenotype), linear regression is performed to first examine the spatial overlap between gene expression and the phenotypic profile. Null-random, null-brain, and null_coexpression models are implemented to examine the gene specificity. Null-spin model is implemented to examine the spatial specificity. 

To get started, please see "./processing/README.txt" for details. Please note that it may take few days to finish all processing, because GAMBA pre-compute results for gene sets with a large range of sizes.

# Examples
We use three simple examples that show analyses commonly performed in literature to illustrate the usage of different statistical null models. Examples include human-supragranular-enriched (HSE) genes, APOE gene, and risk genes of autism spectrum disorder (ASD). Relevant scripts are included in ./examples

# Simulation
We simulate and analyze the outcome of different statistical evaluation approaches for a wide range of real brain phenotypes and artificial gradients. Relevant scripts are included in ./simulation
