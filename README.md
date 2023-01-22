
# Visualizing next-generation sequencing variant calling data with Shiny

This repository contains an R Shiny app for exploring a benchmarking experiment for calling poliovirus intrahost variants with deep sequencing. This experiment is described in [Valesano et al. *Cell Host and Microbe* 2020](https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(20)30574-6).

Experiment Overview: Wild-type poliovirus type 1 was mixed into Sabin1 virus at five nucleic acid concentrations (10^5, 10^4, 10^3, 10^2, and 10^1 genome copies/uL) and six known frequencies (0%, 1%, 2%, 5%, 10%, and 100% WT). These two viruses differ by 66 single nucleotide substitutions. A subset of samples were spiked into  stool-derived nucleic acids to simulate the complex mixture present in human samples. RNA was extracted with QIAamp Viral RNA kit, genome amplified in four segments by RT-PCR, and sequencing libraries generated with Nextera DNA Flex library preparation kit along with plasmid controls. Sequencing was performed on an Illumina MiSeq with 2x250 paired-end reads. Variants were called using deepSNV using a previously-available pipeline: https://github.com/lauringlab/variant_pipeline

This is the base folder for the Shiny app for the polio benchmarking experiment. Data folder contains raw and reference files. All scripts, including app.R, as well as any processing scripts, are in this base folder.

# Run the App

To run this app, need to have the "shiny" package installed in R, as well as the dependencies below:

--------
> packages <- c("shiny", "tidyverse", "ggplot2", "pROC", "wesanderson", "gtable", "gridExtra", "plyr", "reshape2")

--------

Then the app can be downloaded and run straight from GitHub. To run the app, start an R session and run the following commands: 

--------
> runGitHub("PolioBenchmarking_Shiny", "andrewvalesano")
--------