#' A function for estimating a surrogate variable
#'
#' iasva.unit() estimates a hidden factor for unwanted variation while accounting for all known factors and tests the significance 
#' of its contribution on the unmodeled variation in the data. This function is called in iasva() function to iteratively identify hidden 
#' factors.   
#'
#'
#' @param Y read counts matrix with samples in row and genes in column.
#' @param X  known variables including the primary variables of interest. 
#' @param intercept If intercept=FALSE, the linear intercept is not included in the model.
#' @param permute If permute=TRUE, a permutation test (Buja and Eyuboglu 1992, Leek and Storey 2008) is conducted to assess the significance of the putative hidden factor.
#' @param num.p number of permutations to be used to calculate the permuation test p-value.
#' @param threads number of cores to be used in permutation test.
#' @param num.sv.permtest num of top singluar values to be used in computing the permutation test statistic.
#' @param tol stopping tolerance for the augmented implicitly restarted Lanczos bidiagonalization algorithm
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#' 
#' @return sv estimated surrogate variable.
#' @return pc.stat.obs PC test statistic value of the estimated surrogate variable. 
#' @return pval permuation p-value of the estimated surrogate variable.

iasva.unit <- function(Y, X, intercept=TRUE, permute=TRUE, num.p=100, threads=1, num.sv.permtest=NULL, tol=1e-10, verbose=FALSE){
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
  
  if(permute==TRUE){
    if(verbose) {cat("\n Assess the significance of the contribution of SV")}
    if(is.null(num.sv.permtest)){
      svd.res.obs <- svd(tresid - rowMeans(tresid))
    } else {
      svd.res.obs <- irlba::irlba(tresid-rowMeans(tresid), num.sv.permtest, tol=tol)
    }

    pc.stat.obs <- svd.res.obs$d[1]^2/sum(svd.res.obs$d^2)
    if(verbose) {cat("\n PC test statistic value:", pc.stat.obs)}
    
    # Generate an empirical null distribution of the PC test statistic.
    pc.stat.null.vec <- rep(0, num.p)
    permute.svd <- permute.svd.factory(lY, X, num.sv.permtest, tol, verbose)
    if (threads > 1) {
      threads <- min(threads, parallel::detectCores()-1)
      cl <- parallel::makeCluster(threads)
      pc.stat.null.vec <- tryCatch(parallel::parSapply(cl, 1:num.p, permute.svd), error=function(err){parallel::stopCluster(cl); stop(err)})
      parallel::stopCluster(cl)
    } else {
      pc.stat.null.vec <- sapply(1:num.p, permute.svd)
    }
    if(verbose) {cat("\n Empirical null distribution of the PC statistic:", sort(pc.stat.null.vec))}
    pval <- sum(pc.stat.obs <= pc.stat.null.vec)/(num.p+1)
    if(verbose) {cat("\n Permutation p-value:", pval)}
  } else {
    pc.stat.obs <- -1
    pval <- -1
  }
  return(list(sv=sv, pc.stat.obs=pc.stat.obs, pval=pval))
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


#' permute.svd.factory
#'  
#' @param lY log-transformed read counts matrix
#' @param X matrix of known factors
#' @param num.sv.permtest num of top singluar values to be used in computing the permutation test statistic.
#' @param tol stopping tolerance for the augmented implicitly restarted Lanczos bidiagonalization algorithm
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#'  
#' @return permute.svd()

permute.svd.factory <- function(lY, X, num.sv.permtest, tol, verbose) {
  
  permute.svd <- function(i) {
    permuted.lY <- apply(t(lY), 1, sample, replace=FALSE)
    tresid.null <- t(resid(.lm.fit(cbind(1,X), permuted.lY)))
    if(is.null(num.sv.permtest)) {
      svd.res.null <- svd(tresid.null)
    } else {
      svd.res.null <- irlba::irlba(tresid.null, num.sv.permtest, tol=tol)
    }
    return(svd.res.null$d[1]^2/sum(svd.res.null$d^2))
  }
  
  return(permute.svd)
}