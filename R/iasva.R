#' A function for iteratively adjusted surrogate variable analysis (IA-SVA)
#'
#' The iterative procedure of IA-SVA is implemented in this function (iasva). iasva() iteratively runs iasva.unit() function 
#' to identify a hidden factor for unwanted variation while accounting for all known factors and test the significance 
#' of its contribution on the unmodeled variation in the data. If the test statistic of detected factor is significant, 
#' iasva() includes the factor as a known variable in the next iteration to find further hidden factors.  
#'
#'
#' @param Y read counts matrix with samples in row and genes in column.
#' @param X  known variables including the primary variables of interest. 
#' @param intercept If intercept=FALSE, the linear intercept is not included in the model.
#' @param num.sv number of surrogate variables to estimate. 
#' @param permute If permute=TRUE, a permutation test (Buja and Eyuboglu 1992, Leek and Storey 2008) is conducted to assess the significance of the putative hidden factor.
#' @param num.p number of permutations to be used to calculate the permuation test p-value.
#' @param sig.cutoff significance threshold for the permutation test
#' @param threads number of cores to be used in permutation test.
#' @param num.sv.permtest num of top singluar values to be used in computing the permutation test statistic.
#' @param tol stopping tolerance for the augmented implicitly restarted Lanczos bidiagonalization algorithm
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#'
#' @return sv matrix of estimated surrogate variables, one column for each surrogate variable. 
#' @return pc.stat.obs vector of PC test statistic values, one value for each surrogate variable. 
#' @return pval vector of permuation p-values, one value for each surrogate variable.
#' @return n.sv number of significant/obtained surrogate variables. 
#' 
#' @examples
#' library(iasva)
#' data(iasva_test)
#' iasva.res <- iasva(iasvaY, iasvaX) 
#'  
#' @export

iasva <- function(Y, X, intercept=TRUE, num.sv=NULL, permute=TRUE, num.p=100, sig.cutoff= 0.05, threads=1, num.sv.permtest=NULL, tol=1e-10, verbose=FALSE){
  cat("IA-SVA running...")
  sv <- NULL
  pc.stat.obs <- NULL
  pval <- NULL
  rsq <- NULL
  isv <- 0
  while(TRUE){
    if(!is.null(num.sv)){
      if(isv==num.sv){ 
        break
      }
    }
    iasva.res <- iasva.unit(Y, X, intercept, permute, num.p, threads, num.sv.permtest, tol, verbose)
    if(iasva.res$pval < sig.cutoff){
      sv <- cbind(sv, iasva.res$sv)
      pc.stat.obs <- cbind(pc.stat.obs, iasva.res$pc.stat.obs)
      pval <- c(pval, iasva.res$pval)
      X <- cbind(X,iasva.res$sv)
    } else { break }
    isv <- isv+1
    cat(paste0("\nSV",isv, " Detected!"))
  }
  if (isv > 0) {
    colnames(sv) <- paste0("SV", 1:ncol(sv))
    cat(paste0("\n# of significant surrogate variables: ",length(pval)))
    return(list(sv=sv, pc.stat.obs=pc.stat.obs, pval=pval, n.sv=length(pval)))
  } else {
    cat ("\nNo significant surrogate variables")
  }
}






