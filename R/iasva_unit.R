#' @importFrom stats .lm.fit cutree dist hclust lm p.adjust resid
#' @importFrom irlba irlba
#' @importFrom parallel detectCores makeCluster parSapply stopCluster
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom SummarizedExperiment SummarizedExperiment assay

iasva_unit <- function(Y, X, intercept = TRUE, permute = TRUE, num.p = 100,
                       threads = 1, num.sv.permtest = NULL,
                       tol = 1e-10, verbose = FALSE) {
  # transpose the read counts
  Y <- t(assay(Y))
  if (min(Y) < 0) {
    Y <- Y + abs(min(Y))
  }
  lY <- log(Y + 1)
  if (intercept) {
    fit <- .lm.fit(cbind(1, X), lY)
  } else {
    fit <- .lm.fit(X, lY)
  }
  resid <- resid(fit)
  tresid <- t(resid)
  if (verbose) {
    cat("\n Perform SVD on residuals")
  }
  svd_pca <- irlba(tresid - rowMeans(tresid), 1, tol = tol)
  if (verbose) {
    cat("\n Regress residuals on PC1")
  }
  fit <- .lm.fit(cbind(1, svd_pca$v[, 1]), resid)
  if (verbose) {
    cat("\n Get Rsq")
  }
  rsq.vec <- calc_rsq(resid, fit)
  if (verbose) {
    cat("\n Rsq 0-1 Normalization")
  }
  rsq.vec[is.na(rsq.vec)] <- min(rsq.vec, na.rm = TRUE)
  wgt <- (rsq.vec - min(rsq.vec)) / (max(rsq.vec) - min(rsq.vec))
  if (verbose) {
    cat("\n Obtain weighted log-transformed read counts")
  }
  tlY <- t(lY) * wgt
  if (verbose) {
    cat("\n Perform SVD on weighted log-transformed read counts")
  }
  sv <- irlba(tlY - rowMeans(tlY), 1, tol = tol)$v[, 1]
  if (permute == TRUE) {
    if (verbose) {
      cat("\n Assess the significance of the contribution of SV")
    }
    if (is.null(num.sv.permtest)) {
      svd.res.obs <- svd(tresid - rowMeans(tresid))
    } else {
      svd.res.obs <- irlba(tresid - rowMeans(tresid),
                                  num.sv.permtest, tol = tol)
    }
    pc.stat.obs <- svd.res.obs$d[1] ^ 2 / sum(svd.res.obs$d ^ 2)
    if (verbose) {
      cat("\n PC test statistic value:", pc.stat.obs)
    }
    # Generate an empirical null distribution of the PC test statistic.
    pc.stat.null.vec <- rep(0, num.p)
    permute.svd <- permute_svd_factory(lY, X, num.sv.permtest, tol, verbose)
    if (threads > 1) {
      #threads <- min(threads, detectCores() - 1)
      #cl <- makeCluster(threads)
      #clusterExport(cl, "irlba")
      #pc.stat.null.vec <- tryCatch(parSapply(cl, 
      #                            seq(from = 1, to = num.p), permute.svd),
      #                             error = function(err) {
      #                               stopCluster(cl); stop(err)
      #                               })
      #stopCluster(cl)
      pc.stat.null.vec <- tryCatch(bplapply(seq(from = 1, to = num.p), permute.svd, BPPARAM=MulticoreParam(workers=threads)),
                                   error=function(err) {invisible(err); stop(err)})
      
    } else {
      pc.stat.null.vec <- sapply(seq(from = 1, to = num.p), permute.svd)
    }
    if (verbose) {
      cat("\n Empirical null distribution of the PC statistic:",
          sort(pc.stat.null.vec))
    }
    pval <- sum(pc.stat.obs <= pc.stat.null.vec) / (num.p + 1)
    if (verbose) {
      cat("\n Permutation p-value:", pval)
    }
  } else {
    pc.stat.obs <- -1
    pval <- -1
  }
  return(list(sv = sv, pc.stat.obs = pc.stat.obs, pval = pval))
}

calc_rsq <- function(resid, fit) {
  RSS <- colSums(resid(fit) ^ 2)
  TSS <- colSums(t(t(resid) - colSums(resid) / ncol(resid)) ^ 2)
  return(1 - (RSS / (nrow(resid) - 2)) / (TSS / (nrow(resid) - 1)))
}

permute_svd_factory <- function(lY, X, num.sv.permtest, tol, verbose) {
  permute.svd <- function(i) {
    permuted.lY <- apply(t(lY), 1, sample, replace = FALSE)
    tresid.null <- t(resid(.lm.fit(cbind(1, X), permuted.lY)))
    if (is.null(num.sv.permtest)) {
      svd.res.null <- svd(tresid.null)
    } else {
      svd.res.null <- irlba(tresid.null, num.sv.permtest, tol = tol)
    }
    return(svd.res.null$d[1] ^ 2 / sum(svd.res.null$d ^ 2))
  }
  return(permute.svd)
}
