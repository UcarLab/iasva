#' A function for iteratively adjusted surrogate variable analysis (IA-SVA)
#'
#' The iterative procedure of IA-SVA is implemented in this function (iasva). iasva() iteratively runs iasva.unit() function 
#' to identify a hidden factor for unwanted variation while accounting for all known factors and test the significance 
#' of its contribution on the unmodeled variation in the data. If the test statistic of detected factor is significant, 
#' iasva() includes the factor as a known variable in the next iteration to find further hidden factors.  
#'
#' @param Y The read counts matrix with samples in row and genes in column.
#' @param X  The known variables including the primary variables of interest. 
#' @param permute If permute=TRUE, a permutation test (Buja and Eyuboglu 1992, Leek and Storey 2008) is conducted to assess the significance of the putative hidden factor.
#' @param num.p The number of permutations that will be used to calculate the permuation test p-value.
#' @param num.sv The number of surrogate variables to estimate. 
#' @param sig.cutoff The significance threshold for the permutation test
#' @param intercept If intercept=FALSE, the linear intercept is not included in the model.
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#'
#' @return sv The matrix of estimated surrogate variables, one column for each surrogate variable. 
#' @return sv.resid The matrix of residuals. 
#' @return pc.stat.obs The vector of PC test statistic values, one value for each surrogate variable. 
#' @return pval The vector of permuation p-values, one value for each surrogate variable.
#' @return wgt The matrix of gene weights (0-1 normalized R-squared values), one column for each surrogate variable.
#' @return rsq The matrix of R-squared values, one column for each surrogate variable.
#' @return n.sv The number of significant/obtained surrogate variables. 
#' 
#' @examples
#' library(iasva)
#' data(iasva_test)
#' iasva(iasvaY, iasvaX)  
#' @export
#'

iasva <- function(Y, X, permute=TRUE, num.p=100, num.sv=NULL, sig.cutoff= 0.05, intercept=TRUE, verbose=FALSE){
  cat("IA-SVA running...")
  wgt <- NULL
  sv <- NULL
  sv.resid <- NULL
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
    iasva.res <- iasva.unit(Y, X, permute, num.p, intercept, verbose)
    if(iasva.res$pval < sig.cutoff){
      wgt <- cbind(wgt, iasva.res$wgt)
      rsq <- cbind(rsq, iasva.res$rsq)
      sv <- cbind(sv, iasva.res$sv)
      sv.resid <- cbind(sv.resid, iasva.res$sv.resid)
      pc.stat.obs <- cbind(pc.stat.obs, iasva.res$pc.stat.obs)
      pval <- c(pval, iasva.res$pval)
      X <- cbind(X,iasva.res$sv)
    } else { break }
    isv <- isv+1
    cat(paste0("\nSV",isv, " Detected!"))
  }
  colnames(sv) <- paste0("SV", 1:ncol(sv))
  colnames(sv.resid) <- paste0("SV", 1:ncol(sv))
  
  row.names(wgt) <- NULL
  cat(paste0("\n# of significant surrogate variables: ",length(pval)))
  return(list(sv=sv, sv.resid=sv.resid, pc.stat.obs=pc.stat.obs, pval=pval, wgt=wgt, rsq=rsq, n.sv=length(pval)))
}

