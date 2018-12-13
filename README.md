# Iteratively Adjusted Surrogate Variable Analysis (IA-SVA)

Iteratively Adjusted Surrogate Variable Analysis (IA-SVA) is a statistical framework to uncover hidden sources of variation even when these sources are correlated with the biological variable of interest. IA-SVA provides a flexible methodology to i) identify a hidden factor for unwanted heterogeneity while adjusting for all known factors; ii) test the significance of the putative hidden factor for explaining the variation in the data; and iii), if significant, use the estimated factor as an additional known factor in the next iteration to uncover further hidden factors. 

## Citing IA-SVA

Our work has been published in _Scientific Reports_: [Detection of correlated hidden factors from single cell transcriptomes using Iteratively Adjusted-SVA (IA-SVA)](https://www.nature.com/articles/s41598-018-35365-9) Donghyung Lee, Anthony Cheng, Nathan Lawlor, Mohan Bolisetty, Duygu Ucar

## Author(s)

Donghyung Lee <donghyung.lee@jax.org> , Anthony Cheng <anthony.cheng@jax.org> , and Nathan Lawlor <nathan.lawlor@jax.org>

## Installation

To install IA-SVA package, start R and enter the following commands:

      library(devtools)
      devtools::install_github("UcarLab/iasva")


## Load the package

To load this package, enter the following command to the R console:

      library(iasva)


## Bioconductor

IA-SVA is also available for download from Bioconductor [(link to development version here)](https://www.bioconductor.org/packages/devel/bioc/html/iasva.html)

To install and load the package, please enter the following into the R console:

      if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
      BiocManager::install("iasva", version = "devel")
      library(iasva)

## View Vignettes

For instructions on how to use IA-SVA, please see the package vignette.

For additional tutorials, please click __Quick View__ for any of the four examples below to view in a web browser. 

##### Example 1) Detecting hidden heterogeneity in human islet alpha cells   ([Quick View](https://dleelab.github.io/iasvaExamples/detecting_hidden_heterogeneity.html)) 


##### Example 2) Detecting cell-cycle stage difference in glioblastoma cells   ([Quick View](https://dleelab.github.io/iasvaExamples/hidden_heterogeneity_glioblastoma.html))


##### Example 3) IA-SVA based feature selection improves the performance of clustering algorithms [1]  ([Quick View](https://dleelab.github.io/iasvaExamples/tSNE_post_IA-SVA_3celltypes.html))


##### Example 4) IA-SVA based feature selection improves the performance of clustering algorithms [2]  ([Quick View](https://dleelab.github.io/iasvaExamples/tSNE_post_IA-SVA_Xin_Islets.html))

##### Example 5) Compare IA-SVA to factor analyses methods in terms of their ability to detect marker genes for different cell types [2]  ([Quick View](https://dleelab.github.io/iasvaExamples/Brain_scRNASeq_neuron_vs_oligodendrocyte_single_run.html))

##### Example 6) scRNA-seq Data Simulation [1] ([Quick View](https://dleelab.github.io/iasvaExamples/scRNAseq_simulation.html))

##### Example 7) scRNA-seq Data Simulation [2] ([Quick View](https://dleelab.github.io/iasvaExamples/sim_study_KF_HF_genes_overlap.html))

## R Shiny Web Application
An interactive web version of IA-SVA is publicly hosted on [shinyapps.io here](https://nlawlor.shinyapps.io/VSVA/) 

The github page for this application is available [here](https://github.com/nlawlor/iasva_shiny/tree/master)

## Disclaimer of Warranties and Liabilities

The Jackson Laboratory provides the software “as is” without warranty of any kind, implied or expressed. You assume full responsibility and risk of loss resulting from your downloading and use of the content of the software. We expressly disclaim any warranty of merchantability, title, security, accuracy and non-infringement. In no event shall The Jackson Laboratory be liable for any claim, damages or other liability arising from the software or the use of the software. You may only use our content in academic research but not for commercial purposes. The software is provided as an information resource only, and should not be used or relied on for any diagnostic or treatment purposes.
