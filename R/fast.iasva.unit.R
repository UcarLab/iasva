#' A function for estimating a surrogate variable (fast IA-SVA) 
#'
#' fast.iasva.unit() estimates a hidden factor for unwanted variation while accounting for all known factors and computes 
#' its contribution (i.e., the percentage of unmodeled variation explained by the hidden factor) on the unmodeled variation in the data. 
#' This function is called in fast.iasva() function to iteratively identify hidden factors.   
#'
#'
#' @param Y read counts matrix with samples in row and genes in column.
#' @param X  known variables. 
#' @param intercept If intercept=FALSE, the linear intercept is not included in the model.
#' @param component.retention If component.retention=TRUE, the percentage of unmodeled variance explained by the putative hidden factor.
#' @param num.sv.retention num of top singular values to be used in computing the percentage of unmodeled variation explained by the putative hidden factor. If num.sv.retention=NULL, all singular values are used.
#' @param tol stopping tolerance for the augmented implicitly restarted Lanczos bidiagonalization algorithm
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#'
#' @return sv estimated surrogate variable.
#' @return percentage percentage of unmodeled variance explained by the estimated surrogate variable.

fast.iasva.unit <- function(Y, X, intercept=TRUE, component.retention=TRUE, num.sv.retention=NULL, tol=1e-10, verbose=FALSE){
  if(min(Y)<0){ Y <- Y + abs(min(Y)) }
  lY <- log(Y+1)
  if(intercept){
    fit <- .lm.fit(cbind(1,X), lY)
  } else {
    fit <- .lm.fit(X, lY)
  }
  resid <- resid(fit)
  tresid = t(resid)
  
  if(verbose) {cat("\n Perform SVD on residuals")}
  svd_pca <- irlba::irlba(tresid-rowMeans(tresid), 1, tol=tol)
  
  if(verbose) {cat("\n Regress residuals on PC1")}
  fit <- .lm.fit(cbind(1,svd_pca$v[,1]), resid)
  
  if(verbose) {cat("\n Get Rsq")}
  rsq.vec <- calc.rsq(resid, fit)
  
  if(verbose) {cat("\n Rsq 0-1 Normalization")}
  rsq.vec[is.na(rsq.vec)] <- min(rsq.vec, na.rm=TRUE)
  wgt <- (rsq.vec-min(rsq.vec))/(max(rsq.vec)-min(rsq.vec)) #0-1 normalization
  
  if(verbose) {cat("\n Obtain weighted log-transformed read counts")}
  tlY = t(lY)*wgt # weigh each row (gene) with respect to its Rsq value.
  
  if(verbose) {cat("\n Perform SVD on weighted log-transformed read counts")}
  sv <- irlba::irlba(tlY-rowMeans(tlY), 1, tol=tol)$v[,1]
  
  percentage <- NULL
  if(component.retention==TRUE){
    if(verbose) {cat("\n Compute the percentage of unmodeled variance explained by SV")}
    if(is.null(num.sv.retention)){
      svd.res.obs <- svd(tresid - rowMeans(tresid))
    } else {
      svd.res.obs <- irlba::irlba(tresid-rowMeans(tresid), num.sv.retention, tol=tol)
    }

    percentage <- (svd.res.obs$d[1]^2/sum(svd.res.obs$d^2))*100
    if(verbose) {cat("\n ",percentage,"% of unmodeled variance is explained by SV")}
    
  } else {
    percentage <- -1
  }
  return(list(sv=sv, percentage=percentage))
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