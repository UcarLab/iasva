library(Seurat)
library(dplyr)
library(Matrix)

# Load the PBMC dataset
pbmc.data <- Read10X("~/Downloads/filtered_gene_bc_matrices/hg19/")

#Examine the memory savings between regular and sparse matrices
dense.size <- object.size(as.matrix(pbmc.data))
dense.size


sparse.size <- object.size(pbmc.data)
sparse.size

dense.size/sparse.size

set.seed(31223133)
# Randomly select 500 cells only
ran.cell500 <- sample(1:2700, 500, replace=FALSE)
pbmc.data <- pbmc.data[,ran.cell500]
dim(pbmc.data)

# Initialize the Seurat object with the raw (non-normalized data)
# Note that this is slightly different than the older Seurat workflow, where log-normalized values were passed in directly.
# You can continue to pass in log-normalized values, just set do.logNormalize=F in the next step.
pbmc <- new("seurat", raw.data = pbmc.data)

# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, do.logNormalize = T, total.expr = 1e4, project = "10X_PBMC")


# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
# NOTE: You must have the Matrix package loaded to calculate the percent.mito values.
mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
VlnPlot(pbmc, c("nGene", "nUMI", "percent.mito"), nCol = 3)

#GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage, and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(pbmc, "nUMI", "percent.mito")
GenePlot(pbmc, "nUMI", "nGene")

#We filter out cells that have unique gene counts over 2,500
#Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nGene but on anything in the object (as in GenePlot above)
pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = 2500)
pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.05)


#note that this overwrites pbmc@scale.data. Therefore, if you intend to use RegressOut, you can set do.scale=F and do.center=F in the original object to save some time.
pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))

pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)

length(pbmc@var.genes)

pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)

VizPCA(pbmc, 1:2)

PCAPlot(pbmc, 1, 2)

PCHeatmap(pbmc, pc.use = 1, cells.use = 100, do.balanced = TRUE)

PCHeatmap(pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = FALSE)

JackStrawPlot(pbmc, PCs = 1:12)

PCElbowPlot(pbmc)

#save.SNN=T saves the SNN so that the  SLM algorithm can be rerun using the same graph, but with a different resolution value (see docs for full details)
pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = T)

pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = T)

TSNEPlot(pbmc)








dim(pbmc@raw.data)



mod1 <- model.matrix(~pbmc@data.info$percent.mito+pbmc@data.info$nUMI)[,-1]
#geo.lib.size <- colSums(log(pbmc@raw.data+1))
#mod1 <- model.matrix(~geo.lib.size)[,-1]

iasva.res<- iasva(t(as.matrix(pbmc@data)), mod1 , verbose=FALSE, permute=FALSE, num.sv=3)
iasva.sv <- iasva.res$sv
pairs(iasva.sv, main="IA-SVA")
cor(iasva.sv)

markers <- find.markers(t(as.matrix(pbmc@raw.data)), iasva.sv)
pheatmap(log(markers+1), show_colnames =FALSE, show_rownames=TRUE, clustering_method = "ward.D2",cutree_cols = 4)
tsne.res <- Rtsne(t(log(markers+1)), dims = 2)
pairs(tsne.res$Y, main="tSNE")

tsne.res <- Rtsne(t(as.matrix(pbmc@raw.data)), dims = 2)
pairs(tsne.res$Y, main="tSNE")







rm(list=ls())
library(Rtsne)
library(iasva)
library(pheatmap)


calc_cpm <- function(expr_mat, spikes = NULL){
    norm_factor<-NULL
    if(is.null(spikes)){
      norm_factor <- colSums(expr_mat)
    } else {
      norm_factor <- colSums(expr_mat[-spikes, ])
    }
    return ((t(t(expr_mat)/norm_factor))*10^6)
}

path.data <- "/Users/leed1/Desktop/My_Project/Data/Dendritic_Cell_Mohan/CD1C_CD141_Mixed/"
all.counts <- read.table(paste0(path.data,"Counts.csv"), header=TRUE, sep=",")
anno <-   read.table(paste0(path.data,"CellClusterAnnotations.csv"), header=FALSE, sep=",")
dim(all.counts)
head(all.counts[,1:5])
dim(anno)
head(anno)
gene.name <- all.counts$X
all.counts <- all.counts[,-1]
rownames(all.counts) <- gene.name
head(all.counts[,1:5])

cell.type <- anno$V2
table(cell.type)


lib.size <- colSums(all.counts)

counts2 <- all.counts[,lib.size>500]
dim(counts2)

anno2 <- anno[lib.size>500,]

filter = apply(counts2, 1, function(x) length(x[x>=1])>=3)  ##2,5
counts2 <- counts2[filter,]
dim(counts2)

counts2.cpm <- calc_cpm(counts2, spikes = NULL)

lib.size <- colSums(log(counts2.cpm+1))
lib.size.cpm <- colSums(counts2.cpm)
plot(lib.size.cpm)
lib.size2 <- colSums(counts2)
summary(lib.size)

cell.type <- anno2$V2

plot(as.factor(cell.type), lib.size.cpm)

cor(lib.size, cell.type)

mod1 <- model.matrix(~lib.size)

#iasva.res<- iasva(t(counts2.cpm), mod1[,-1],verbose=FALSE,permute=FALSE,num.sv=3)
iasva.res<- iasva2(t(counts2.cpm), NULL, verbose=FALSE,permute=FALSE,num.sv=3)

iasva.sv <- iasva.res$sv
pairs(iasva.sv, main="IA-SVA", col=cell.type)

plot(lib.size,-iasva.sv[,1], col=cell.type)
plot(lib.size2,-iasva.sv[,1], col=cell.type)

cor(iasva.sv)
corrplot(cor(batch_iasva), order="hclust")
corrplot(cor(cbind(lib.size, iasva.sv)))

tsne.res <- Rtsne(t(counts2), dims = 2)
pairs(tsne.res$Y, main="tSNE", col=cell.type)

cor(cbind(lib.size, iasva.sv))


markers.all <- find.markers(t(counts2.cpm), as.matrix(iasva.sv[,2]))
rownames(markers.all)

pheatmap(markers.all, show_colnames =FALSE, show_rownames=TRUE, clustering_method = "ward.D2",cutree_cols = 2)


tsne.res <- Rtsne(unique(t(markers.all)), dims = 2, perplexity = 7)
pairs(tsne.res$Y, main="tSNE")

## CD1C only
counts <- all.counts[,cell.type<2]
dim(counts)

lib.size <- colSums(counts)

counts <- counts[,lib.size>5000]
dim(counts)

zero.prop <- sum(counts==0)/length(as.matrix(counts))
zero.prop #95.9%

filter = apply(counts, 1, function(x) length(x[x>=1])>=3)  ##2,5
counts <- counts[filter,]
dim(counts)

counts.cpm <- calc_cpm(counts, spikes = NULL)

zero.prop <- sum(counts==0)/length(as.matrix(counts))
zero.prop

Num_Exp <- colSums(log(counts+1))
Num_Exp2 <- colSums(counts)
summary(Num_Exp)
sum(Num_Exp>2000)
ldat <- log(counts + 1)
batch_pca = svd(counts.cpm - rowMeans(counts.cpm))$v
pairs(batch_pca[,1:10], main="PCA")

cor(batch_pca[,1], Num_Exp)

batch_tsne <- Rtsne(t(counts.cpm), dims = 3)
pairs(batch_tsne$Y, main="tSNE")


batch_iasva.recur<- iasva2(t(counts.cpm), NULL, verbose=FALSE,permute=FALSE,num.sv=20, log=FALSE)
batch_iasva <- batch_iasva.recur$sv
cor(batch_iasva)[,1]<0.2
pairs(batch_iasva, main="IA-SVA")
corrplot(cor(batch_iasva), order="hclust")

cor(cbind(Num_Exp,batch_iasva))[,1]<0.2

marker.counts <- find.markers2(t(counts.cpm), as.matrix(batch_iasva), log=FALSE)
rownames(marker.counts)

marker.counts <- find.markers(t(counts), as.matrix(batch_iasva[,(cor(cbind(Num_Exp,batch_iasva))[,1]<0.2)[-1]]))
pheatmap(marker.counts, show_colnames =FALSE, show_rownames=TRUE, clustering_method = "ward.D2",cutree_cols = 4)

batch_tsne <- Rtsne(unique(t(marker.counts)), dims = 2)
pairs(batch_tsne$Y, main="tSNE")

## This code find markers for IA-SVA detected hidden factors
## This should be added in the package as a function.
lY <- log(t(counts)+1)
all.markers <- NULL
for(i in 1:4){
  fit <- lm(lY~batch_iasva[,i])
  pval.vec <- unlist(lapply(summary(fit), function(x) x$coefficient[2,4]))
  rsq.vec <- unlist(lapply(summary(fit), function(x) x$adj.r.squared))
  fdr.vec <- p.adjust(pval.vec, method="BH", n= length(pval.vec))
  markers <- rownames(counts)[fdr.vec < 0.01 & rsq.vec > 0.6]
  cat(length(markers),"\n")
  all.markers <- c(all.markers, markers)
}
length(all.markers)
all.markers <- unique(all.markers)
length(all.markers)
marker.counts <- counts[rownames(counts)%in%all.markers,]
dim(marker.counts)
##


#marker.counts[marker.counts==2] <- 1

pheatmap(marker.counts, show_colnames =FALSE, show_rownames=TRUE, clustering_method = "ward.D2",cutree_cols = 4)

mldat0 <- log(marker.counts+1)
tmldat0 <- unique(t(mldat0))
batch_tsne <- Rtsne(tmldat0, dims = 3)
pairs(batch_tsne$Y, main="tSNE")

