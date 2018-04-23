# Iteratively Adjusted Surrogate Variable Analysis (IA-SVA)

Iteratively Adjusted Surrogate Variable Analysis (IA-SVA) is a statistical framework to uncover hidden sources of variation even when these sources are correlated with the biological variable of interest. IA-SVA provides a flexible methodology to i) identify a hidden factor for unwanted heterogeneity while adjusting for all known factors; ii) test the significance of the putative hidden factor for explaining the variation in the data; and iii), if significant, use the estimated factor as an additional known factor in the next iteration to uncover further hidden factors. 

## Citing IA-SVA

__A robust statistical framework to detect multiple sources of hidden variation in single-cell transcriptomes__, Donghyung Lee, Anthony Cheng, Duygu Ucar, bioRxiv. 2017; doi: https://doi.org/10.1101/151217

## Author

Donghyung Lee <donghyung.lee@jax.org> and Anthony Cheng <anthony.cheng@jax.org>

## Installation

To install IA-SVA package, start R and enter the following commands:

      library(devtools)
      devtools::install_github("UcarLab/IA-SVA")


## Load the package

To load this package, enter the following command to the R console:

      library(iasva)


## View Vignettes

For instructions on how to use IA-SVA, please see the package vignette.

## Disclaimer of Warranties and Liabilities

The Jackson Laboratory provides the software “as is” without warranty of any kind, implied or expressed. You assume full responsibility and risk of loss resulting from your downloading and use of the content of the software. We expressly disclaim any warranty of merchantability, title, security, accuracy and non-infringement. In no event shall The Jackson Laboratory be liable for any claim, damages or other liability arising from the software or the use of the software. You may only use our content in academic research but not for commercial purposes. The software is provided as an information resource only, and should not be used or relied on for any diagnostic or treatment purposes.
