# Multi-cohort cerebrospinal fluid proteomics identifies robust molecular signatures across the Alzheimer disease continuum

# Table of contents
* [Introduction](#introduction)
* [Content](#content)
* [Data](#data)
* [Requirements](#requirements)
* [License](#license)
* [Instructions](#instructions)
* [Directories](#Directories)

# Introduction
This repository contains the code for bioinformatics analyses described in the article "Multi-cohort cerebrospinal fluid proteomics identifies robust molecular signatures across the Alzheimer disease continuum".

This project investigated CSF proteomics data from the SomaScan 7K platform to identify proteins associated with Alzheimer disease. Idnetified proteins were leveraged to create AD-spcific prediction model, pseudo-trajectory analyseis, biological pathway and cell type enrichment analyses to understand underlying AD biology.

# Content
The code covers the following main analysis steps:

1. Data pre-processing: Proteomics data preparation and surrogate variable (SV) computation 
2. Differential expression analysis (Discovery, Replication, and Meta-analyses)
3. AD prediction model development using Lasso regression
4. Survival analysis to identify individuals that will conert to AD
5. AD progression analysis to distinguish between slow and fast progressors
6. Pseudo trajectory analysis to group/cluster proteins based on their expressin in AT continuum (A-T-, A+T-, A+T+)
7. Network and pathway enrichment analyses
8. Cell type enrichment analysis
   
# Data
Proteomics data analysed in this study is available at:
- ADNI: http://adni.loni.usc.edu/
- Knight-ADRC: https://dss.niagads.org/ (Accession: ng00130)
- FACE and Barcelona-1 cohorts: http://www.fundacioace.com/
- PPMI: https://www.ppmi-info.org/
- Stanford-ADRC: https://web.stanford.edu/group/adrc/cgi-bin/web-proj/datareq.php

# Requirements
The code was written in R (version 4.3.0) and relies on multiple R and Bioconductor packages, including:
- sva
- clusterProfiler 
- scran
- glmnet
- nlme
- pROC
- igraph
- survminer
- mclust
- Additional packages listed at the beginning of each R script

# License
The code is available under the MIT License.

# Instructions
The code was tested on R 4.3.0 on Linux operating systems, but should be compatible with later versions of R installed on current Linux, Mac, or Windows systems.

To run the code, the correct working directory containing the input data must be specified at the beginning of the R-scripts, otherwise the scripts can be run as-is.

The scripts should be run in the following order:

    data_preparation.R

    differential_expression_analysis.R

    prediction_modeling.R

    survival_analysis.R

    progression_analysis.R

    clustering_pseudo_trajectories.R

    network_and_pathway_analysis.R

    cell_type_enrichment.R

# Directories
- Input data should be placed under "Input_Data" folder
- Scripts should be placed under "Scripts" folder, otherwise the scripts can be run as-is.
- Figures will be generated under "Figures" folder
