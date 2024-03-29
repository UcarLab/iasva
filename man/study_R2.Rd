% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/study_R2.R
\name{study_R2}
\alias{study_R2}
\title{A function to study different values of R2}
\usage{
study_R2(Y, iasva.sv, selected.svs = 2, no.clusters = 2, verbose = FALSE)
}
\arguments{
\item{Y}{A SummarizedExperiment or SingleCellExperiment class containing
read counts where rows represent genes and columns represent samples.}

\item{iasva.sv}{matrix of estimated surrogate variables,
one column for each surrogate variable.}

\item{selected.svs}{list of SVs that are selected for the
analyses. Default is SV2}

\item{no.clusters}{No of clusters to be used in the analyses.
Default is 2.}

\item{verbose}{If verbose = TRUE, the function outputs detailed messages.}
}
\value{
a summary plot that represents silhoutte index and marker gene counts
 as a function of R2 and corresponding matrices.
}
\description{
study_R2() studies how different R2 thresholds is changing:
 1) number of marker genes; 
 2) clustering quality (assuming number of clusters is known). 
It generated diagnostic plots that shows how selected genes and
 clustering quality changes as a function of R2 threshold.
}
\examples{
counts_file <- system.file("extdata", "iasva_counts_test.Rds",
 package = "iasva")
counts <- readRDS(counts_file)
anns_file <- system.file("extdata", "iasva_anns_test.Rds",
 package = "iasva")
 anns <- readRDS(anns_file)
Geo_Lib_Size <- colSums(log(counts + 1))
Patient_ID <- anns$Patient_ID
mod <- model.matrix(~Patient_ID + Geo_Lib_Size)
summ_exp <- SummarizedExperiment::SummarizedExperiment(assays = counts)
iasva.res<- iasva(summ_exp, mod[, -1],verbose = FALSE, 
permute = FALSE, num.sv = 5)
iasva.sv <- iasva.res$sv
study_res <- study_R2(summ_exp, iasva.sv)

}
