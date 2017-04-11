#' A function for estimating a surrogate variable
#'
#' iasva.unit() estimates a hidden factor for unwanted variation while accounting for all known factors and tests the significance 
#' of its contribution on the unmodeled variation in the data. This function is called in iasva() function to iteratively identify hidden 
#' factors.   
#'
#' @param Y The read counts matrix with samples in row and genes in column.
#' @param X  The known variables including the primary variables of interest. 
#' @param permute If permute=TRUE, a permutation test (Buja and Eyuboglu 1992, Leek and Storey 2008) is conducted to assess the significance of the putative hidden factor.
#' @param num.p The number of permutations that will be used to calculate the permuation test p-value.
#' @param intercept If intercept=FALSE, the linear intercept is not included in the model.
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#'
#' @return sv The estimated surrogate variable.
#' @return pc.stat.obs The PC test statistic value of the estimated surrogate variable. 
#' @return pval The permuation p-value of the estimated surrogate variable.
#' 
#' @examples
#' 
#' @export
#'


iasva.unit <- function(Y, X, permute=TRUE, num.p=100, intercept=TRUE, verbose=FALSE){
  if(min(Y)<0){ Y <- Y + abs(min(Y)) }
  lY <- log(Y+1)
  if(intercept){
    fit <- lm(lY~X)
  } else {
    fit <- lm(lY~X-1)
  }
  resid <- resid(fit)
  tresid = t(resid)
  svd_pca <- svd(tresid-rowMeans(tresid))
  if(verbose) {cat("\n Run Linear Regression to Get Rsq")}
  fit <- lm(resid ~ svd_pca$v[,1])
  if(verbose) {cat("\n Get Rsq")}
  rsq.vec <- unlist(lapply(summary(fit), function(x) x$adj.r.squared))
  if(verbose) {cat("\n Rsq 0-1 Normalization")}
  rsq.vec[is.na(rsq.vec)] <- min(rsq.vec, na.rm=TRUE)
  wgt <- (rsq.vec-min(rsq.vec))/(max(rsq.vec)-min(rsq.vec)) #0-1 normalization
  if(verbose) {cat("\n Conduct SVD on weighted log-transformed read counts")}
  tlY = t(lY)*wgt # weigh each row (gene) with respect to its Rsq value.
  
  sv <- svd(tlY - rowMeans(tlY))$v[,1]
  
  if(verbose) {cat("\n Assess the significance of the contribution of SV")}
  svd.res.obs <- svd(tresid - rowMeans(tresid))
  pc.stat.obs <- svd.res.obs$d[1]^2/sum(svd.res.obs$d^2)
  
  if(permute==TRUE){
    # generate an empirical null distribution of the PC test statistic.
    pc.stat.null.vec <- rep(0, num.p)
    for(i in 1:num.p){
      if(verbose) {cat("\n For-Loop: Permute")}
      permuted.lY <- apply(t(lY), 1, sample, replace=FALSE)
      if(verbose) {cat("\n For-Loop: Get Residuals")}
      tresid.null <- t(resid(lm(permuted.lY~X)))
      if(verbose) {cat("\n For-Loop: Conduct SVD")}
      svd.res.null <- svd(tresid.null)
      if(verbose) {cat("\n For-Loop: Compute PC test statistic")}
      pc.stat.null.vec[i] <- svd.res.null$d[1]^2/sum(svd.res.null$d^2)
    }
    if(verbose) {cat("\n Compute permutation p-value")}
    pval <- sum(pc.stat.obs <= pc.stat.null.vec)/(num.p+1)
  } else {
    pval <- -1
  }
  return(list(sv=sv, pc.stat.obs=pc.stat.obs, pval=pval))
}
