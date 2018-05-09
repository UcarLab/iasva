#' A function for iteratively adjusted surrogate variable analysis (IA-SVA)
#'
#' The iterative procedure of IA-SVA is implemented in this function (iasva).
#' iasva() function iteratively runs iasva_unit() function 
#' to identify a hidden factor for unwanted variation while accounting for
#'  all known factors and test the significance of its contribution on the
#'  unmodeled variation in the data. If the test statistic of detected factor
#'  is significant, iasva() includes the factor as a known variable in the 
#'  next iteration to find further hidden factors.  
#'  
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @param Y A SummarizedExperiment class containing read counts where
#' rows represent genes and columns represent samples.
#' @param X  A model matrix of known variables
#'  including the primary variables of interest. 
#' @param intercept If intercept = FALSE, the linear intercept
#'  is not included in the model.
#' @param num.sv number of surrogate variables to estimate. 
#' @param permute If permute = TRUE, a permutation test (Buja and Eyuboglu 1992,
#'  Leek and Storey 2008) is conducted to assess the significance of 
#'  the putative hidden factor.
#' @param num.p number of permutations to be used to calculate the
#'  permuation test p-value.
#' @param sig.cutoff significance threshold for the permutation test
#' @param threads number of cores to be used in permutation test.
#' @param num.sv.permtest num of top singular values to be used in
#'  computing the permutation test statistic. If num.sv.permtest = NULL,
#'   all singular values are used.
#' @param tol stopping tolerance for the augmented implicitly restarted
#'  Lanczos bidiagonalization algorithm
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#'
#' @return sv matrix of estimated surrogate variables, one column
#'  for each surrogate variable. 
#' @return pc.stat.obs vector of PC test statistic values, 
#' one value for each surrogate variable. 
#' @return pval vector of permuation p-values, 
#' one value for each surrogate variable.
#' @return n.sv number of significant/obtained surrogate variables. 
#' 
#' @examples
#' counts_file <- system.file("extdata", "iasva_counts_test.Rds",
#'  package = "iasva")
#' counts <- readRDS(counts_file)
#' anns_file <- system.file("extdata", "iasva_anns_test.Rds",
#'  package = "iasva")
#'  anns <- readRDS(anns_file)
#' Geo_Lib_Size <- colSums(log(counts + 1))
#' Patient_ID <- anns$Patient_ID
#' mod <- model.matrix(~Patient_ID + Geo_Lib_Size)
#' summ_exp <- SummarizedExperiment::SummarizedExperiment(assays = counts)
#' iasva.res<- iasva(summ_exp, mod[, -1],verbose = FALSE, 
#'  permute = FALSE, num.sv = 5)
#' @export

iasva <- function(Y, X, intercept = TRUE, num.sv = NULL, permute = TRUE,
                  num.p = 100, sig.cutoff = 0.05, threads = 1,
                  num.sv.permtest = NULL, tol = 1e-10, verbose = FALSE) {
  # error handling
  stopifnot(class(Y)[1] == "SummarizedExperiment", is.numeric(tol),
            is.matrix(X), is.numeric(num.p), is.numeric(sig.cutoff),
            is.numeric(threads))
  message("IA-SVA running...")
  sv <- NULL
  pc.stat.obs <- NULL
  pval <- NULL
  isv <- 0
  while (TRUE) {
    if (!is.null(num.sv)) {
      if (isv == num.sv) {
        break
      }
    }
    iasva.res <- iasva_unit(Y, X, intercept, permute, num.p,
                            threads, num.sv.permtest, tol, verbose)
    if (iasva.res$pval < sig.cutoff) {
      sv <- cbind(sv, iasva.res$sv)
      pc.stat.obs <- cbind(pc.stat.obs, iasva.res$pc.stat.obs)
      pval <- c(pval, iasva.res$pval)
      X <- cbind(X, iasva.res$sv)
    } else {
      break
    }
    isv <- isv + 1
    message("\nSV ", isv, " Detected!")
  }
  if (isv > 0) {
    colnames(sv) <- paste0("SV", seq(from = 1, to = ncol(sv)))
    message("\n# of significant surrogate variables: ", length(pval))
    return(list(sv = sv, pc.stat.obs = pc.stat.obs,
                pval = pval, n.sv = length(pval)))
  } else {
    message("\nNo significant surrogate variables")
  }
}