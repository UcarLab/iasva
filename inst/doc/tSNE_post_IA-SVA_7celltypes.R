## ----load_packages, echo=TRUE--------------------------------------------
rm(list=ls())
library(iasva)
library(iasvaExamples)
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
counts <- counts[, (anns$Cell_Type!="none")] 
anns <- subset(anns, (Cell_Type!="none"))
dim(counts)
dim(anns)

anns <- droplevels(anns)

prop.zeros <- sum(counts==0)/length(counts)
prop.zeros

filter = apply(counts, 1, function(x) length(x[x>5])>=3)

counts = counts[filter,]
dim(counts)

prop.zeros <- sum(counts==0)/length(counts)
prop.zeros

Patient_ID <- anns$Patient_ID
Cell_Type <- anns$Cell_Type
Phenotype <- anns$Phenotype
Batch <- anns$Batch

## ----geno_lib_size, echo=TRUE, fig.width=6, fig.height=4-----------------
Geo_Lib_Size <- colSums(log(counts+1))
barplot(Geo_Lib_Size, xlab="Cell", las =2)
lcounts <- log(counts + 1)
pca.res = svd(lcounts - rowMeans(lcounts))$v
cor(Geo_Lib_Size, pca.res[,1])

## ----run_tsne, echo=TRUE, fig.width=6, fig.height=6----------------------
tsne.res <- Rtsne(t(lcounts), dims = 2)

plot(tsne.res$Y, main="tSNE", pch=21, bg=c("red","green3","blue","grey","orange","black","purple")[Cell_Type], oma=c(4,4,6,12), cex=1.5)
legend("topleft", levels(Cell_Type), fill=c("red", "green3", "blue","grey","orange","black","purple"), bty="n")

## ----run_iasva, echo=TRUE, fig.width= 6, fig.height=5--------------------
mod <- model.matrix(~Phenotype+Patient_ID+Batch+Geo_Lib_Size)
iasva.res<- iasva(t(counts), mod[,-1],verbose=FALSE, permute=FALSE, num.sv=6)
iasva.sv <- iasva.res$sv

pairs(iasva.sv, main="tSNE", pch=21, bg=c("red","green3","blue","grey","orange","black","purple")[Cell_Type], oma=c(4,4,6,12))
legend(0.82, 0.65, levels(Cell_Type), fill=c("red", "green3", "blue","grey","orange","black","purple"), bty="n")

## ----find_markers, echo=TRUE, fig.width=6, fig.height=5------------------
marker.counts <- find.markers(t(counts), as.matrix(iasva.sv[,c(1,3,4,6)]),  rsq.cutoff = 0.4)
nrow(marker.counts)

anno.col <- data.frame(Cell_Type=Cell_Type)
rownames(anno.col) <- colnames(marker.counts)
head(anno.col)
pheatmap(log(marker.counts+1), show_colnames =FALSE, show_rownames = TRUE, clustering_method = "ward.D2",cutree_cols = 7,annotation_col = anno.col)

## ----run_tsne_post_iasva, echo=TRUE, fig.width=6, fig.height=6-----------
set.seed(75458456)
tsne.res <- Rtsne(unique(t(log(marker.counts+1))), dims = 2)
plot(tsne.res$Y, main="tSNE", pch=21, bg=c("red","green3","blue","grey","orange","black","purple")[Cell_Type], oma=c(4,4,6,12), cex=1.5)
legend("topright", levels(Cell_Type), fill=c("red", "green3", "blue","grey","orange","black","purple"), bty="n")

## ----session_info, echo=TRUE---------------------------------------------
sessionInfo()

