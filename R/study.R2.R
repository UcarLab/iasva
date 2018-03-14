#' A function to study different values of R2 
#' 
#' study.R2() studies how different R2 thresholds is changging: 1) number of marker genes; 2) clustering quality (assuming number of clusters is known). 
#' It generated diagnostic plots that shows how selected genes and clustering quality changes as a function of R2 threshold.
#'  
#'
#' @param Y read counts matrix with samples in row and genes in column.
#' @param iasva.sv  matrix of estimated surrogate variables, one column for each surrogate variable. 
#' @param selected.svs list of SVs that are selected for the analyses. Default is SV2
#' @param no.clusters No of clusters to be used in the analyses. Default is 2. 
#'
#' @return a summary plot that represents silhoutte index and marker gene counts as a function of R2 and corresponding matrices. 
#' 
#' 
#' @export

study.R2 <- function(Y, iasva.sv, selected.svs=2, no.clusters=2, verbose=FALSE){
  C.scores <- matrix(0,0,0)
  Number.of.genes <- matrix(0,0,0)
  for (i in seq(0.1,0.9,0.05)){
    marker.counts <- find.markers(Y, as.matrix(iasva.sv[,selected.svs]), rsq.cutoff = i)
    no.genes <- dim(marker.counts)[1]
    if(no.genes == 0){
      break
    }
    else{
      my.dist <- dist(t(log(marker.counts+1)))
      my.clustering <- hclust(my.dist, method = "ward.D2")
      my.silhoutte <-silhouette(cutree(my.clustering,no.clusters),my.dist)
      C1 <- mean(my.silhoutte[my.silhoutte[,1]==1,3])
      C2 <- mean(my.silhoutte[my.silhoutte[,1]==2,3])
      average.C <- (C1+C2)/2
      C.scores <- c(C.scores, average.C)
      Number.of.genes <- c(Number.of.genes,no.genes)
    }
  }
  
  output.matrix <- rbind(C.scores, Number.of.genes)
  end.point <- (length(C.scores)-1)*0.05+0.1
  colnames(output.matrix) <- seq(0.1,end.point,0.05)
 
  par(mar=c(5,5,5,5))
  plot(Number.of.genes,  xlab = "R^2", ylab = "Number genes selected", xaxt="n", main = "Number of selected genes vs. Cluster quality", pch = 18, col ="blue",type="b", lty=2, cex=2)
  Axis(1,at=seq(1,length(Number.of.genes)), side = 1, labels=seq(0.1,end.point,0.05),las = 2)
  par(new=T)
  plot(C.scores, xlab='', ylab='', axes=F, pch = 18 , col ="red",type="b", lty=2, cex=2)
  Axis(side=4)
  mtext(side = 4, line = 2, 'Average Silhoutte Score', col = "red")
  par(new=F)

  return(output.matrix)
}