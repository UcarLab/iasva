# Iteratively Adjusted Surrogate Variable Analysis (IA-SVA)

Iteratively Adjusted Surrogate Variable Analysis (IA-SVA) is a statistical framework to uncover hidden sources of variation even when these sources are correlated with the biological variable of interest. IA-SVA provides a flexible methodology to i) identify a hidden factor for unwanted heterogeneity while adjusting for all known factors; ii) test the significance of the putative hidden factor for explaining the variation in the data; and iii), if significant, use the estimated factor as an additional known factor in the next iteration to uncover further hidden factors. 

## Citing IA-SVA

__A robust statistical framework to detect multiple sources of hidden variation in single-cell transcriptomes__, Donghyung Lee, Anthony Cheng, Duygu Ucar, bioRxiv. 2017; doi: https://doi.org/10.1101/151217

## Installation

To install IA-SVA package, start R and enter the following commands:

      library(devtools)
      devtools::install_github("UcarLab/IA-SVA")


## Load the package

To load this package, enter the following command to the R console:

      library(iasva)


## View Vignettes

Click __Quick View__ to view each vignette in a web browser. 


##### Example 1) Detecting hidden heterogeneity in human islet alpha cells   ([Quick View](https://cdn.rawgit.com/dleelab/iasvaExamples/d4d63ce3/inst/doc/detecting_hidden_heterogeneity_iasvaV0.95.html)) 


##### Example 2) Detecting cell-cycle stage difference in glioblastoma cells   ([Quick View](https://cdn.rawgit.com/dleelab/iasvaExamples/d4d63ce3/inst/doc/hidden_heterogeneity_glioblastoma_iasvaV0.95.html))


##### Example 3) IA-SVA based feature selection improves the performance of clustering algorithms [1]  ([Quick View](https://cdn.rawgit.com/dleelab/iasvaExamples/d4d63ce3/inst/doc/tSNE_post_IA-SVA_3celltypes_iasvaV0.95.html))


##### Example 4) IA-SVA based feature selection improves the performance of clustering algorithms [2]  ([Quick View](https://cdn.rawgit.com/dleelab/iasvaExamples/d4d63ce3/inst/doc/tSNE_post_IA-SVA_Xin_Islets_iasvaV0.95.html))


