This is the base folder for the Shiny app for the polio benchmarking experiment.

Data folder contains raw and reference files.

All scripts, including app.R, as well as any processing and ROC scripts, are in this base folder.

To run this app, need to have the "shiny" package installed in R. Run the command: 

>install.packages("shiny")

Then, start an R session and type the following commands:

>library(shiny)
>runGitHub("PolioBenchmarking_Shiny", "andrewvalesano")
