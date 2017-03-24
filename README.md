# Iteratively Adjusted Surrogate Variable Analysis (IA-SVA)

Iteratively Adjusted Surrogate Variable Analysis (IA-SVA) is a statistical framework to uncover hidden sources of variation even when these sources are correlated. IA-SVA provides a flexible methodology to i) identify a hidden factor for unwanted heterogeneity while adjusting for all known factors; ii) test the significance of the putative hidden factor for explaining the variation in the data; and iii), if significant, use the estimated factor as an additional known factor in the next iteration to uncover further hidden factors.

## Installation

To install this package, start R and enter the following commands:

      if(!require(devtools)){
            install.packages("devtools")
            library(devtools)
      }
      install_github("UcarLab/IA-SVA")

## Load the package

To load this package, enter the following command to the R console:

      library(iasva)

## View Vignettes

##### Example 1) Detecting hidden heterogeneity in single cell RNA seq data

      vignette("detecting_hidden_heterogeneity", package="iasva")

##### Example 2) tSNE post IA-SVA [1]

      vignette("tSNE_post_IA-SVA_3celltypes", package="iasva")  

##### Example 3) tSNE post IA-SVA [2]

      vignette("tSNE_post_IA-SVA_7celltypes", package="iasva")  
  