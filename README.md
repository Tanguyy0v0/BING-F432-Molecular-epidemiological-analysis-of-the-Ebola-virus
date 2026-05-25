# BING-F432-Molecular-epidemiological-analysis-of-the-Ebola-virus

This repository contains the files required to reconstruct the spatio-temporal dispersal dynamics of a simulated Ebola Virus (EBOV) outbreak. This project was developed as part of the BING-F432 Spatial and Temporal Epidemiology course.
The objective is to use genomic data associated with metadata to infer the evolutionary and geographic history of the pathogen.

# Workflow
The analysis follows a standard genomic epidemiology pipeline

1) Initial Phylogeny Estimation : Using IQ-TREE to obtain a Maximum Likelihood tree ans selectting the best-fitting substitution model via the BIC criterion

2) Temporal Signal Evaluation : Performing a root-to-tip regression analysis in TempEst to confirm tha the viral population is measuably evolving

3) Bayesien inference : Configuring the model parameters in BEAUTi and running the MCMC in BEAST X to obtain a time-scaled phylogeny and continuous geographic reconstuction

4) Spatio-temporal analysis : Using the R Package seraphim to extract dispersal trajectories and estimate dispersal statistics

# Repository Content 
- Output_rapport_2 : 
- data : Genomic sequences in FASTA format
- phylogeny : .treefile and initial ML trees from IQ-TREE
- XML : configuration dile used for the BEAST analysis 
- results : postrior distribution of trees (.tree); and log files (.log), and the MCC tree
- scripts : R scripts for data processing ans seraphim visualiations

  # Installation
  In R, seraphim can be installed with the devtools package:
  
  install.packages("devtools"); library(devtools)
install_github("sdellicour/seraphim/unix_OS") # (for a Unix OS)
install_github("sdellicour/seraphim/windows") # (for a Windows OS)

# Package references
Dellicour S, Faria N, Rose R, Lemey P, Pybus OG (2026). SERAPHIM 2.0: an extended toolbox for studying phylogenetically informed movements. Bioinformatics 42: btag093


