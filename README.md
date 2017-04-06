# Iteratively Adjusted Surrogate Variable Analysis (IA-SVA)

Iteratively Adjusted Surrogate Variable Analysis (IA-SVA) is a statistical framework to uncover hidden sources of variation even when these sources are correlated with the biological variable of interest. IA-SVA provides a flexible methodology to i) identify a hidden factor for unwanted heterogeneity while adjusting for all known factors; ii) test the significance of the putative hidden factor for explaining the variation in the data; and iii), if significant, use the estimated factor as an additional known factor in the next iteration to uncover further hidden factors. 

## Installation

To install IA-SVA package and the package with example datasets (iasvaExamples), start R and enter the following commands:

      if(!require(devtools)){
            install.packages("devtools")
            library(devtools)
      }
      install_github("UcarLab/IA-SVA")
      
      if(!require(devtools)){
        install.packages("devtools")
        library(devtools)
        }
     install_github("dleelab/iasvaExamples")


## Load the package

To load this package, enter the following command to the R console:

      library(iasva)
      library(iasvaExamples)

## View Vignettes

After installation, enter each command to the R console to display the vignette rendered in a viewer or click __Quick View__ to view each vignette in a web browser. 

##### Example 1) Detecting hidden heterogeneity in single cell RNA seq data   ([Quick View](http://htmlpreview.github.io/?https://github.com/UcarLab/IA-SVA/blob/master/inst/doc/detecting_hidden_heterogeneity.html))

      vignette("detecting_hidden_heterogeneity", package="iasva")

##### Example 2) tSNE post IA-SVA [1]   ([Quick View](http://htmlpreview.github.io/?https://github.com/UcarLab/IA-SVA/blob/master/inst/doc/tSNE_post_IA-SVA_3celltypes.html))

      vignette("tSNE_post_IA-SVA_3celltypes", package="iasva")  

##### Example 3) tSNE post IA-SVA [2]   ([Quick View](http://htmlpreview.github.io/?https://github.com/UcarLab/IA-SVA/blob/master/inst/doc/tSNE_post_IA-SVA_7celltypes.html))

      vignette("tSNE_post_IA-SVA_7celltypes", package="iasva")  
