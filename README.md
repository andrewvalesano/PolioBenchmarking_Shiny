
# Overview

This repository contains an R Shiny app for exploring the deep sequencing data for benchmarking variant calling in poliovirus. 

Wild-type poliovirus type 1 was mixed into Sabin1 virus at five nucleic acid concentrations (10^5, 10^4, 10^3, 10^2, and 10^1 genome copies/uL) and six known frequencies (0%, 1%, 2%, 5%, 10%, and 100% WT). These two viruses differ by 66 single nucleotide substitutions. RNA was extracted with QIAamp Viral RNA kit, genome amplified in four segments by RT-PCR, and sequencing libraries generated with Nextera DNA Flex library preparation kit along with plasmid controls. Sequencing was performed on an Illumina MiSeq with 2x250 paired-end reads. Variants were called using deepSNV using a previously-available pipeline: https://github.com/lauringlab/variant_pipeline

This is the base folder for the Shiny app for the polio benchmarking experiment. Data folder contains raw and reference files. All scripts, including app.R, as well as any processing scripts, are in this base folder.

# Run the App

To run this app, need to have the "shiny" package installed in R. Then the app can be downloaded and run straight from GitHub. To run the app, start an R session and run the following commands: 

--------
>install.packages("shiny")
>library(shiny)
>runGitHub("PolioBenchmarking_Shiny", "andrewvalesano")
--------