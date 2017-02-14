# IA-SVA
Iteratively Adjusted Surrogate Variable Analysis (IA-SVA) is a statistical framework to uncover hidden sources of variation
even when these sources are correlated. IA-SVA provides a flexible methodology to i) identify a hidden factor for unwanted heterogeneity
while adjusting for all known factors; ii) test the significance of the putative hidden factor for explaining the variation in the data;
and iii), if significant, use the estimated factor as an additional known factor in the next iteration to uncover further hidden factors.

## Installation
1. Create a tarball of this repository after cloning it `tar zcvf iasva.tar.gz IA-SVA.git`
2. Install the R package using `R CMD install iasva.tar.gz`

## Using library functions
1. Load the library using `library(iasva)`