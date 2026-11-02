# Week 8 methylation homework

This folder contains Austin Lee's completed QBIO 490 Week 8 assignment.

- `Methylation_Homework.Rmd`: R Notebook source with the complete workflow, code, outputs, and interpretations.
- `Methylation_Homework.nb.html`: rendered notebook for review in a browser.
- `figures/`: the required volcano, box, and bar plots.
- `results/`: full compressed result tables, compact tables, and plot data.
- `carc/`: the R pipeline and Slurm scripts used for the official USC CARC run.

The final analysis was completed on USC CARC Discovery as Slurm job 13663807 using TCGA-BRCA data downloaded from the NCI Genomic Data Commons with `TCGAbiolinks`. It used 793 primary and 5 metastatic methylation samples and 782 primary and 5 metastatic matched RNA samples. The raw GDC data are not committed because they are large public files; the notebook and `carc/` scripts document the reproducible workflow.
