#' A function for fast IA-SVA 
#'
#' The iterative procedure of fast IA-SVA is implemented in this function (fast.iasva). fast.iasva() iteratively  
#' identify a hidden factor for unwanted variation while accounting for all known factors, compute its contribution 
#' (i.e., the percentage of unmodeled variation explained by the hidden factor) on the unmodeled variation in the data. 
#' If the contribution is greater than a user-defined cutoff (pct.cutoff, default = 1%), 
#' the factor is retained and used as a known variable in the next iteration to find further hidden factors.  
#'
#'
#' @param Y read counts matrix with samples in row and genes in column.
#' @param X  known variables. 
#' @param intercept If intercept=FALSE, the linear intercept is not included in the model.
#' @param num.sv number of surrogate variables to estimate. 
#' @param pct.cutoff percetage threshold for SV retention. IA-SVA computes the percentage of unmodeled variance explained by the putative hidden factor and compare it with the user-defined threshold. If the percentage is greater than the threshold, SV is retained.
#' @param num.tsv num of top singular values to be used in computing the percentage of unmodeled variation explained by the putative hidden factor. If num.tsv=NULL, all singular values are used.
#' @param tol stopping tolerance for the augmented implicitly restarted Lanczos bidiagonalization algorithm
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#'
#' @return sv matrix of estimated surrogate variables, one column for each surrogate variable. 
#' @return pct vector of percentages of unmodeled variance explained by each surrogate variable, one value for each surrogate variable.
#' @return n.sv number of obtained surrogate variables. 
#' 
#' @examples
#' library(iasva)
#' data(iasva_test)
#' iasva.res <- fast.iasva(iasvaY, iasvaX) 
#'  
#' @export

fast.iasva <- function(Y, X, intercept=TRUE, num.sv=NULL, pct.cutoff= 1, num.tsv=NULL, tol=1e-10, verbose=FALSE){
  cat("fast IA-SVA running...")
  if(min(Y)<0){ Y <- Y + abs(min(Y)) }
  lY <- log(Y+1)
  svd.all <- NULL
  if(is.null(num.tsv)){
    svd.all <- svd(t(lY)-rowMeans(t(lY)))
  } else {
    svd.all <- irlba::irlba(t(lY)-rowMeans(t(lY)), num.tsv, tol=tol)
  }
  sv <- NULL
  pct <- NULL
  isv <- 0
  while(TRUE){
    if(!is.null(num.sv)){
      if(isv==num.sv){ 
        break
      }
    }
    if(intercept){
      fit <- .lm.fit(cbind(1,X), lY)
    } else {
      fit <- .lm.fit(X, lY)
    }
    resid <- resid(fit)
    tresid = t(resid)
    
    if(verbose) {cat("\n Perform SVD on residuals")}
    svd.resid <- irlba::irlba(tresid-rowMeans(tresid), 1, tol=tol)
    if(verbose) {cat("\n Regress residuals on PC1")}
    fit <- .lm.fit(cbind(1,svd.resid$v[,1]), resid)
    if(verbose) {cat("\n Get Rsq")}
    rsq.vec <- calc.rsq(resid, fit)
    if(verbose) {cat("\n Rsq 0-1 Normalization")}
    rsq.vec[is.na(rsq.vec)] <- min(rsq.vec, na.rm=TRUE)
    wgt <- (rsq.vec-min(rsq.vec))/(max(rsq.vec)-min(rsq.vec)) #0-1 normalization
    if(verbose) {cat("\n Obtain weighted log-transformed read counts")}
    tlY = t(lY)*wgt # weigh each row (gene) with respect to its Rsq value.
    if(verbose) {cat("\n Perform SVD on weighted log-transformed read counts")}
    sv.i <- irlba::irlba(tlY-rowMeans(tlY), 1, tol=tol)$v[,1] # estimated sv
    
    if(verbose) {cat("\n Compute the percentage of unmodeled variance explained by SV")}
    pct.i <- (svd.resid$d[1]^2/sum(svd.all$d^2))*100
    if(verbose) {cat("\n ",pct.i,"% of unmodeled variance is explained by SV")}
    
    if(pct.i >= pct.cutoff){
      sv <- cbind(sv, sv.i)
      pct <- c(pct, pct.i)
      X <- cbind(X, sv.i)
    } else { break }
    isv <- isv+1
    cat(paste0("\nSV",isv, " Detected!"))
  }
  if (isv > 0) {
    colnames(sv) <- paste0("SV", 1:ncol(sv))
    cat(paste0("\n# of obtained surrogate variables: ",length(pct)))
    return(list(sv=sv, pct=pct, n.sv=length(pct)))
  } else {
    cat ("\nNo surrogate variables obtained")
  }
}

#' calc.rsq
#'  
#' @param resid residual matrix
#' @param fit output of .lm.fit()
#'  
#' @return R squared

calc.rsq <- function(resid, fit) {
  RSS <- colSums(resid(fit)^2)
  TSS <- colSums(t(t(resid) - colSums(resid)/ncol(resid)) ^ 2)  # vectorized
  # TSS <- colSums(sweep(resid, 2, mean, "-") ^ 2)  # alt-2
  return(1-(RSS/(nrow(resid)-2))/(TSS/(nrow(resid)-1)))
}