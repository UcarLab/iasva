#' A function for fast IA-SVA 
#'
#' The iterative procedure of fast IA-SVA is implemented in this function (fast.iasva). fast.iasva() iteratively runs fast.iasva.unit() function 
#' to identify a hidden factor for unwanted variation while accounting for all known factors, compute its contribution 
#' (i.e., the percentage of unmodeled variation explained by the hidden factor) on the unmodeled variation in the data. 
#' If the contribution is greater than a user-defined cutoff (retention.cutoff, default = 5%), 
#' fast.iasva() includes the factor as a known variable in the next iteration to find further hidden factors.  
#'
#'
#' @param Y read counts matrix with samples in row and genes in column.
#' @param X  known variables. 
#' @param intercept If intercept=FALSE, the linear intercept is not included in the model.
#' @param num.sv number of surrogate variables to estimate. 
#' @param component.retention If component.retention=TRUE, the percentage of unmodeled variance explained by the putative hidden factor.
#' @param retention.cutoff percetage threshold for component retention.
#' @param num.sv.retention num of top singular values to be used in computing the percentage of unmodeled variation explained by the putative hidden factor. If num.sv.retention=NULL, all singular values are used.
#' @param tol stopping tolerance for the augmented implicitly restarted Lanczos bidiagonalization algorithm
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#'
#' @return sv matrix of estimated surrogate variables, one column for each surrogate variable. 
#' @return percentage vector of percentages of unmodeled variance explained by each surrogate variable, one value for each surrogate variable.
#' @return n.sv number of obtained surrogate variables. 
#' 
#' @examples
#' library(iasva)
#' data(iasva_test)
#' iasva.res <- fast.iasva(iasvaY, iasvaX) 
#'  
#' @export

fast.iasva <- function(Y, X, intercept=TRUE, num.sv=NULL, component.retention=TRUE, retention.cutoff= 5, num.sv.retention=NULL, tol=1e-10, verbose=FALSE){
  cat("fast IA-SVA running...")
  sv <- NULL
  percentage <- NULL
  rsq <- NULL
  isv <- 0
  while(TRUE){
    if(!is.null(num.sv)){
      if(isv==num.sv){ 
        break
      }
    }
    iasva.res <- fast.iasva.unit(Y, X, intercept, component.retention, num.sv.retention, tol, verbose)
    if(iasva.res$percentage < retention.cutoff){
      sv <- cbind(sv, iasva.res$sv)
      percentage <- c(percentage, iasva.res$percentage)
      X <- cbind(X,iasva.res$sv)
    } else { break }
    isv <- isv+1
    cat(paste0("\nSV",isv, " Detected!"))
  }
  if (isv > 0) {
    colnames(sv) <- paste0("SV", 1:ncol(sv))
    cat(paste0("\n# of obtained surrogate variables: ",length(percentage)))
    return(list(sv=sv, percentage=percentage, n.sv=length(percentage)))
  } else {
    cat ("\nNo surrogate variables obtained")
  }
}