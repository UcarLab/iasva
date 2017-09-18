#' A function for finding markers for hidden factors
#' 
#' find.markers() takes a read counts matrix of entire gene set and a matrix of surrogate variables estimated by IA-SVA as input, 
#' identifies marker genes highly correlated with each surrogate variable and returns a read counts matrix of the markers.
#'  
#'
#' @param Y read counts matrix with samples in row and genes in column.
#' @param iasva.sv  matrix of estimated surrogate variables, one column for each surrogate variable. 
#' @param method multiple testing adjustment method, default="BH".
#' @param sig.cutoff significance cutoff.
#' @param rsq.cutoff R squared cutoff.
#' @param verbose If verbose=TRUE, the function outputs detailed messages. 
#'
#' @return marker.counts read counts matrix of markers, one column for each cell. 
#' 
#' @examples
#' library(iasva)
#' data(iasva_test)
#' iasva.res <- iasva(iasvaY, iasvaX) 
#' markers <- find.markers(iasvaY, iasva.res$sv)
#' 
#' @export

find.markers <- function(Y, iasva.sv, method="BH", sig.cutoff=0.05, rsq.cutoff=0.3, verbose=FALSE){
  if(min(Y)<0){ Y <- Y + abs(min(Y)) }
  lY <- log(Y+1)
  all.markers <- NULL
  num.sv <- ncol(iasva.sv)
  for(i in 1:num.sv){
    fit <- lm(lY~iasva.sv[,i])
    pval.vec <- unlist(lapply(summary(fit), function(x) x$coefficient[2,4]))
    rsq.vec <- unlist(lapply(summary(fit), function(x) x$adj.r.squared))
    pval.vec[is.na(pval.vec)] <- 1
    rsq.vec[is.na(rsq.vec)] <- 0
    adj.pval.vec <- p.adjust(pval.vec, method=method, n= length(pval.vec))
    markers <- colnames(Y)[adj.pval.vec < sig.cutoff & rsq.vec > rsq.cutoff]
    cat(paste0("# of markers (",colnames(iasva.sv)[i],"): ", length(markers),"\n"))
    all.markers <- c(all.markers, markers)
  }
  all.markers <- unique(all.markers)
  cat("total # of unique markers: ",length(all.markers))
  marker.counts <- t(Y[,colnames(Y)%in%all.markers])
  return(marker.counts)
}

