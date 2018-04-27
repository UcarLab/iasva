#' A function for finding markers for hidden factors
#' 
#' Function takes a read counts matrix of entire gene set and a matrix
#'  of surrogate variables estimated by IA-SVA as input, 
#'  identifies marker genes highly correlated with each surrogate variable
#'  and returns a read counts matrix of the markers.
#' @importFrom stats lm p.adjust
#'  
#' @param Y read counts matrix with samples in row and genes in column.
#' @param iasva.sv  matrix of estimated surrogate variables,
#'  one column for each surrogate variable. 
#' @param method multiple testing adjustment method, default = "BH".
#' @param sig.cutoff significance cutoff.
#' @param rsq.cutoff R squared cutoff.
#' @param verbose If verbose = TRUE, the function outputs detailed messages.
#'
#' @return marker.counts read counts matrix of markers,
#'  one column for each cell.
#' 
#' @examples
#' library(iasva)
#' counts_file <- system.file("extdata", "iasva_counts_test.Rds",
#'  package = "iasva")
#' counts <- readRDS(counts_file)
#' anns_file <- system.file("extdata", "iasva_anns_test.Rds",
#'  package = "iasva")
#'  anns <- readRDS(anns_file)
#' Geo_Lib_Size <- colSums(log(counts + 1))
#' Patient_ID <- anns$Patient_ID
#' mod <- model.matrix(~Patient_ID + Geo_Lib_Size)
#' iasva.res <- iasva(t(counts), mod[, -1], num.sv = 5, permute = FALSE)
#' markers <- find_markers(t(counts), iasva.res$sv)
#' @export

find_markers <- function(Y, iasva.sv, method = "BH", sig.cutoff = 0.05,
                         rsq.cutoff = 0.3, verbose = FALSE) {
  if (min(Y) < 0) {
    Y <- Y + abs(min(Y))
  }
  lY <- log(Y + 1)
  all.markers <- NULL
  num.sv <- ncol(iasva.sv)
  for (i in seq(from = 1, to = num.sv, by = 1)) {
    fit <- lm(lY ~ iasva.sv[, i])
    pval.vec <- unlist(lapply(summary(fit), function(x) x$coefficient[2, 4]))
    rsq.vec <- unlist(lapply(summary(fit), function(x) x$adj.r.squared))
    pval.vec[is.na(pval.vec)] <- 1
    rsq.vec[is.na(rsq.vec)] <- 0
    adj.pval.vec <- p.adjust(pval.vec, method = method, n = length(pval.vec))
    markers <- colnames(Y)[adj.pval.vec < sig.cutoff & rsq.vec > rsq.cutoff]
    cat(paste0("# of markers (", colnames(iasva.sv)[i], "): ",
               length(markers), "\n"))
    all.markers <- c(all.markers, markers)
  }
  all.markers <- unique(all.markers)
  cat("total # of unique markers: ", length(all.markers))
  marker.counts <- t(Y[, colnames(Y) %in% all.markers])
  return(marker.counts)
}
