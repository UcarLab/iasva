\name{iasva}
\alias{iasva}
\title{Iteratively Adjusted Surrogate Variable Analysis}
\description{IA-SVA can be used to uncover hidden heterogeneity in bulk and single cell transcriptomic data}

\usage{
iasva(Y, X, permute=TRUE, num.p=100, num.sv=NULL, sig.cutoff=0.05, intercept=TRUE, verbose=FALSE)
}

\arguments{
\item{Y}{Samples (m) by genes (n) read counts matrix}
\item{X}{Primary variable of interest}
\item{permute}{Allow permutation}
\item{num.p}{Number of permutations allowed}
\item{num.sv}{Number of surrogate variables to estimate}
\item{sig.cutoff}{Significance threshold}
\item{intercept}{Include the linear intercept}
\item{verbose}{Output all error and non-essential messages}
}

\author{Donghyun Lee, Anthony Cheng, Duygu Ucar}
\examples{
library(iasva)
data(iasva_test)
iasva(iasvaY, iasvaX)
}
