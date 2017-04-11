## ----install_packages, echo=TRUE, eval=FALSE-----------------------------
#  #devtools
#  if(!require(devtools)){
#          install.packages("devtools")
#          library(devtools)
#  }
#  #IA-SVA
#  install_github("UcarLab/IA-SVA")
#  #IAsvaexamples
#  install_github("dleelab/iasvaExamples")
#  
#  #SVA
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("sva")
#  
#  # Rtsne
#  install.packages("Rtsne",repos = "http://cran.us.r-project.org")
#  
#  #pheatmap
#  install.packages("phetmap",repos = "http://cran.us.r-project.org")

## ----load_packages, echo=TRUE, message=FALSE-----------------------------
rm(list=ls())
library(iasva)
library(iasvaExamples)
library(sva)
library(Rtsne)
library(pheatmap)

## ----load_data, echo=TRUE------------------------------------------------
data("Lawlor_Islet_scRNAseq_Read_Counts")
data("Lawlor_Islet_scRNAseq_Annotations")
ls()
counts <- Lawlor_Islet_scRNAseq_Read_Counts
anns <- Lawlor_Islet_scRNAseq_Annotations
dim(anns)
dim(counts)
# set seed for the random number generator.
set.seed(344532456)

## ----alpha_cells, echo=TRUE, results='asis'------------------------------
counts <- counts[, (anns$Phenotype!="Non-Diabetic")&(anns$Cell_Type=="GCG")] 
anns <- subset(anns, (Phenotype!="Non-Diabetic")&(Cell_Type=="GCG"))
dim(counts)
dim(anns)

anns <- droplevels(anns)

prop.zeros <- sum(counts==0)/length(counts)
prop.zeros

# filter out genes that are sparsely and lowly expressed
filter = apply(counts, 1, function(x) length(x[x>5])>=3)

counts = counts[filter,]
dim(counts)

prop.zeros <- sum(counts==0)/length(counts)
prop.zeros

## ----geno_lib_size, echo=TRUE, fig.width=6, fig.height=4-----------------
Geo_Lib_Size <- colSums(log(counts+1))
barplot(Geo_Lib_Size, xlab="Cell", ylab="Geometric Lib Size", las=2)
lcounts <- log(counts + 1)
pca.res = svd(lcounts - rowMeans(lcounts))$v
# PC1 and Geometric library size correlation
cor(Geo_Lib_Size, pca.res[,1])

## ----run_iasva, echo=TRUE, fig.width= 6, fig.height=5--------------------
Patient_ID <- anns$Patient_ID
mod <- model.matrix(~Patient_ID+Geo_Lib_Size)
iasva.res<- iasva(t(counts), mod[,-1],verbose=FALSE, permute=FALSE, num.sv=3)
iasva.sv <- iasva.res$sv
plot(iasva.sv[,1], iasva.sv[,2], xlab="SV1", ylab="SV2")

Cell_Type <- as.factor(iasva.sv[,2] > -0.2) 
levels(Cell_Type)=c("Cell1","Cell2")

# We identified 6 outlier cells based on SV2 that are marked in red
pairs(iasva.sv, main="IA-SVA", pch=21, bg=c("red","green3")[Cell_Type], oma=c(4,4,6,12)) #4,4,6,12
legend(0.8, 0.55, levels(Cell_Type), fill=c("red", "green3"), bty="n")

## ----find_markers, echo=TRUE, fig.width=8, fig.height=12-----------------
marker.counts <- find.markers(t(counts), as.matrix(iasva.sv[,2]))
nrow(marker.counts)

anno.col <- data.frame(Cell_Type=Cell_Type)
rownames(anno.col) <- colnames(marker.counts)
head(anno.col)
pheatmap(log(marker.counts+1), show_colnames =FALSE, clustering_method = "ward.D2",cutree_cols = 2,annotation_col = anno.col)

## ----run_tsne, echo=TRUE, fig.width=6, fig.height=6----------------------
set.seed(323542534)
tsne.res <- Rtsne(t(lcounts), dims = 2)

plot(tsne.res$Y, main="tSNE", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, bg=c("red","green3")[Cell_Type], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Type), fill=c("red", "green3"), bty="n")

## ----run_tsne_post_iasva, echo=TRUE, fig.width=6, fig.height=6-----------
set.seed(345233)
tsne.res <- Rtsne(unique(t(log(marker.counts+1))), dims = 2)

plot(tsne.res$Y, main="tSNE post IA-SVA", xlab="tSNE Dim1", ylab="tSNE Dim2", pch=21, bg=c("red","green3")[Cell_Type], oma=c(4,4,6,12))
legend("bottomright", levels(Cell_Type), fill=c("red", "green3"), bty="n")

## ----run_pca, echo=TRUE, fig.width=6, fig.height=5-----------------------
pca_res = svd(lcounts - rowMeans(lcounts))$v

pairs(pca.res[,1:3], main="PCA", pch=21, bg=c("red","green3")[Cell_Type], oma=c(4,4,6,12)) #4,4,6,12
legend(0.8, 0.55, levels(Cell_Type), fill=c("red", "green3"), bty="n")

## ----run_sva, echo=TRUE, fig.width=6, fig.height=5-----------------------
mod1 <- model.matrix(~Patient_ID+Geo_Lib_Size)
mod0 <- cbind(mod1[,1])

sva.res = svaseq(counts,mod1,mod0, n.sv=3)$sv
##plot(sva.res[,1], sva.res[,2], xlab="SV1", ylab="SV2", col=Cell_Type)
pairs(sva.res[,1:3], main="SVA", pch=21, bg=c("red","green3")[Cell_Type], oma=c(4,4,6,12)) #4,4,6,12
legend(0.8, 0.55, levels(Cell_Type), fill=c("red", "green3"), bty="n")

## ----session_info, echo=TRUE---------------------------------------------
sessionInfo()

