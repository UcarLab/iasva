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
  #svd_pca <- svd(tresid)
  svd_pca <- svd(tresid-rowMeans(tresid))
  sv.resid <- svd_pca$v[,1]
  if(verbose) {cat("\n Run Linear Regression to Get Rsq")}
  fit <- lm(resid ~ svd_pca$v[,1])
  if(verbose) {cat("\n Get Rsq")}
  rsq.vec <- unlist(lapply(summary(fit), function(x) x$adj.r.squared))
  if(verbose) {cat("\n Rsq 0-1 Normalization")}
  rsq.vec[is.na(rsq.vec)] <- min(rsq.vec, na.rm=TRUE)
  wgt <- (rsq.vec-min(rsq.vec))/(max(rsq.vec)-min(rsq.vec)) #0-1 normalization
  ##plot(wgt)
  if(verbose) {cat("\n Conduct SVD on weighted log-transformed read counts")}
  tlY = t(lY)*wgt # weight each row using normalized rsq values
  
  sv <- svd(tlY - rowMeans(tlY))$v[,1]
  #svs <- svd(tlY - rowMeans(tlY))$v
  #max.index <- which.max(cor(svs, sv.resid))
  #sv <- svs[,max.index]
  
  if(verbose) {cat("\n Assess the significance of the contribution of SV")}
  svd.res.obs <- svd(tresid - rowMeans(tresid))
  pc.stat.obs <- svd.res.obs$d[1]^2/sum(svd.res.obs$d^2)
  
  if(permute==TRUE){
    ## test significance here
    pc.stat.null.vec <- rep(0, num.p)
    for(i in 1:num.p){
      if(verbose) {cat("\n For-Loop: Permute")}
      permuted.lY <- apply(t(lY), 1, sample, replace=FALSE)
      if(verbose) {cat("\n For-Loop: Get Residuals")}
      tresid.null <- t(resid(lm(permuted.lY~X)))
      if(verbose) {cat("\n For-Loop: Conduct SVD")}
      svd.res.null <- svd(tresid.null)
      if(verbose) {cat("\n For-Loop: Get PC test statistic")}
      pc.stat.null.vec[i] <- svd.res.null$d[1]^2/sum(svd.res.null$d^2)
    }
    if(verbose) {cat("\n Get permutation p-value")}
    pval <- sum(pc.stat.obs <= pc.stat.null.vec)/(num.p+1)
  } else {
    pval <- -1
  }
  return(list(sv=sv, sv.resid=sv.resid, pc.stat.obs=pc.stat.obs, pval=pval, wgt=wgt, rsq=rsq.vec))
}

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
  row.names(wgt) <- NULL
  cat(paste0("\n# of significant surrogate variables: ",length(pval)))
  return(list(sv=sv, sv.resid=sv.resid, pc.stat.obs=pc.stat.obs, pval=pval, wgt=wgt, rsq=rsq, n.sv=length(pval)))
}

