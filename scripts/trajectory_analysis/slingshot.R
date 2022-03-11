# BiocManager::install("slingshot")
# BiocManager::install("tradeSeq")
# BiocManager::install("TSCAN")
# ?slingshot::slingshot

library(slingshot)
library(tradeSeq)
library(mclust, quietly = TRUE)
gmt_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
deg_stat <- genes_df$morans_I
names(deg_stat) <- genes_df$gene_short_name
save_dir <- output_dir



# TO DO 
# Get cluster label: optimize number of clusters using GMM
# Partition function later
# Get PCA using Seurat function without scale data 
# Get slingshot trajectory 
# Return pseudotime 
# Run DE analysis, regress out cell cycle




BPPARAM <- BiocParallel::bpparam()
BPPARAM 

sce <- SingleCellExperiment(assays = List(counts = counts(cds)))
dim(counts(cds))

# sce <- SingleCellExperiment(assays = List(counts = counts))
# assays(sce)$logcounts <- logcounts(cds)
dim(reducedDims(cds)[['PCA']][,1:20])
dim(reducedDims(cds)[['UMAP']])
length(colData(cds)$cluster_label)
length(unique(colData(cds)$cluster_label))
colData(cds)$cluster_label[1:3]
reducedDims(sce) <- SimpleList(PCA = reducedDims(cds)[['PCA']], UMAP = reducedDims(cds)[['UMAP']])
colData(sce)$cluster_label <- colData(cds)$cluster_label
unique(colData(sce)$cluster_label)
## ----sling_sce----------------------------------------------------------------
sce <- slingshot(sce, clusterLabels = 'cluster_label', 
                 reducedDim = 'UMAP', approx_points = 150)  #

## ----plot_curve_1-------------------------------------------------------------
summary(sce$slingPseudotime_1)

sce$cluster_label <- as.factor(sce$cluster_label)
cls <- sce$cluster_label
names(cls) <- colnames(sce)
head(cls)
# plot(reducedDims(sce)$PCA, col =  brewer.pal(5,"Set1")[cls], pch=16, asp = 1)
# lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black') # discrete curve
library(RColorBrewer)
unique(cls)
col_map <- c('red','yellow','green','cyan','pink')
names(col_map) <- c("Cls_1","Cls_25","Cls_3","Cls_4","Cls_6")

png(paste0(output_dir,"slingshot_5clusters_",datatag,".png"), height = 2*450, width=2*500,res = 2*72)
plot(reducedDims(sce)$UMAP, col =  col_map[cls], pch=16, asp = 1) #brewer.pal(5,"Set1")
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black') # discrete curve
dev.off()
t <- unique(cls)
col_map <- colorRampPalette(brewer.pal(8, "Set1"))(length(t))
names(col_map) <- t
png(paste0(output_dir,"slingshot_15clusters_",datatag,".png"), height = 2*450, width=2*500,res = 2*72)
plot(reducedDims(sce)$UMAP, col =  col_map[cls], pch=16, asp = 1) #
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black') # discrete curve
dev.off()

lin1 <- getLineages(reducedDims(cds)[['UMAP']], sce$cluster_label, start.clus = 'Cls_4')
lin1
crv1 <- getCurves(lin1, approx_points=150) # default: 100 
?getCurves
crv1
class(crv1)
rd <- reducedDims(cds)[['UMAP']]
cl <- sce$cluster_label
plot(rd, col = brewer.pal(10,"Paired")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')

lin1 <- getLineages(reducedDims(cds)[['UMAP']], sce$cluster_label, start.clus = 'Cls_4')

crv1 <- getCurves(lin1, approx_points =20) # default: 100 

set.seed(6)
sce1 <- fitGAM(counts = as.matrix(counts(cds)), sds = crv1, 
               nknots = 6, verbose = FALSE, parallel=TRUE, BPPARAM = BPPARAM)
saveRDS(sce1, paste0(output_dir,'fitGAM.rds'))

BPPARAM$workers <- 5 # use 2 cores
sce1 <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE, parallel=TRUE, BPPARAM = BPPARAM)
mean(rowData(sce)$tradeSeq$converged)
# test for dynamic expression
ATres <- associationTest(sce1)
summary(ATres$pvalue)
sum(ATres$pvalue<0.05)
View(head(ATres))
saveRDS(ATres, paste0(output_dir,'ATres.rds'))
ATres <- readRDS(paste0(output_dir,'ATres.rds'))
class(ATres)
ATres <- ATres %>%
  dplyr::filter(!is.na(pvalue))

ATres <- ATres %>%
  dplyr::filter(pvalue<0.05)
dim(ATres)
View(ATres[1:10,])

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:50]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$cluster_label[pst.ord]
heatclus <- as.numeric(gsub('Cls_','',heatclus))
class(heatdata)
length(heatclus)
heatmap(log1p(as.matrix(heatdata)), Colv = NA,
        ColSideColors = brewer.pal(5,"Set1")[heatclus])


pvalLineage <- getSmootherPvalues(sce1)
statLineage <- getSmootherTestStats(sce1)

genes_df <- data.table::fread(paste0(output_dir,'pseudotime_genes_thrs_0.75qt.csv.gz')) %>% as.data.frame()
dim(genes_df)
View(head(genes_df))
summary(genes_df$p_value)

genes_df$morans_I

heatmap(log1p(heatdata), 
        ColSideColors = brewer.pal(9,"Set1")[heatclus])
## ----options, results="hide", include=FALSE, cache=FALSE, message=FALSE-------
knitr::opts_chunk$set(fig.align="center", cache=TRUE,error=FALSE, #stop on error
                      fig.width=5, fig.height=5, autodep=TRUE,
                      results="markup", echo=TRUE, eval=TRUE)
#knitr::opts_knit$set(stop_on_error = 2L) #really make it stop
#knitr::dep_auto()
options(getClass.msg=FALSE)
graphics:::par(pch = 16, las = 1)
set.seed(12345) ## for reproducibility
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)

## ----dataSetup_sim------------------------------------------------------------
# generate synthetic count data representing a single lineage
means <- rbind(
  # non-DE genes
  matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
         ncol = 300, byrow = TRUE),
  # early deactivation
  matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # late deactivation
  matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # early activation
  matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # late activation
  matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # transient
  matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
         ncol = 300, byrow = TRUE)
)
counts <- apply(means,2,function(cell_means){
  total <- rnbinom(1, mu = 7500, size = 4)
  rmultinom(1, total, cell_means)
})
rownames(counts) <- paste0('G',1:750)
colnames(counts) <- paste0('c',1:300)
sce <- SingleCellExperiment(assays = List(counts = counts))

## ----data_sling---------------------------------------------------------------
# library(slingshot, quietly = FALSE)
data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl

dim(rd) # data representing cells in a reduced dimensional space
length(cl) # vector of cluster labels
head(rd)
# cl[1:3]
## ----genefilt-----------------------------------------------------------------
# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

## ----norm---------------------------------------------------------------------
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

## ----pca, cache=TRUE----------------------------------------------------------
pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

## ----umap, cache=TRUE---------------------------------------------------------
library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

## ----add_RDs, cache=TRUE------------------------------------------------------
reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

## ----clustering_mclust--------------------------------------------------------
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
head(cl1)
colData(sce)$GMM <- cl1
?Mclust
library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

## ----clustering---------------------------------------------------------------
cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

## ----sling_sce----------------------------------------------------------------
sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')

## ----plot_curve_1-------------------------------------------------------------
summary(sce$slingPseudotime_1)


library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')  # smooth curve

## ----plot_curve_2-------------------------------------------------------------
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black') # discrete curve

sce$cluster_label <- as.factor(sce$cluster_label)
cls <- sce$cluster_label
names(cls) <- colnames(sce)
head(cls)
plot(reducedDims(sce)$PCA, col =  brewer.pal(5,"Set1")[cls], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black') # discrete curve

## ----tradeseq, eval=FALSE-----------------------------------------------------
#  library(tradeSeq)
#  
#  # fit negative binomial GAM
#  sce <- fitGAM(sce)
#  
#  # test for dynamic expression
#  ATres <- associationTest(sce)

## ----heatmaps, eval=FALSE-----------------------------------------------------
#  topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
#  pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
#  heatdata <- assays(sce)$counts[topgenes, pst.ord]
#  heatclus <- sce$GMM[pst.ord]
#  
#  heatmap(log1p(heatdata), Colv = NA,
#          ColSideColors = brewer.pal(9,"Set1")[heatclus])

## ----heatmapsREAL, echo=FALSE, fig.height=7-----------------------------------
topgenes <- paste0('G',501:750)
dim(assays(sce)$counts)
# tradeSeq has too many dependencies (174 at the time of this writing), but I 
# promise I actually ran it and got this result. This is a *very* clean example 
# dataset
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

## ----sling_lines_unsup--------------------------------------------------------
lin1 <- getLineages(rd, cl, start.clus = '1')
lin1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')

## ----lines_sup_end------------------------------------------------------------
lin2 <- getLineages(rd, cl, start.clus= '1', end.clus = '3')

plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)

## ----curves-------------------------------------------------------------------
?getCurves
crv1 <- getCurves(lin1, approx_points = 150)
crv1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')

## ----sling_approxpoints-------------------------------------------------------
sce5 <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA',
                  approx_points = 5)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce5$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce5)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce5), lwd=2, col='black')

## ----sling_omega--------------------------------------------------------------
rd2 <- rbind(rd, cbind(rd[,2]-12, rd[,1]-6))
cl2 <- c(cl, cl + 10)
pto2 <- slingshot(rd2, cl2, omega = TRUE, start.clus = c(1,11))

plot(rd2, pch=16, asp = 1,
     col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2])
lines(SlingshotDataSet(pto2), type = 'l', lwd=2, col='black')

colData(sce)$slingshot$pseudotime
## ----sling_multtraj-----------------------------------------------------------
plot(rd2, pch=16, asp = 1,
     col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2])
lines(SlingshotDataSet(pto2), lwd=2, col='black')

## ----session------------------------------------------------------------------
sessionInfo()



library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")

class(crv)
head(crv)
counts[1:3,1:3]
set.seed(5)
icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, nGenes = 200, verbose = T)

set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
View(head(pseudotime))
cellWeights <- slingCurveWeights(crv)
View(head(cellWeights))
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE, sce=TRUE)
class(sce)
table(rowData(sce)$tradeSeq)
assoRes <- associationTest(sce)
head(assoRes)
?associationTest
startRes <- startVsEndTest(sce)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts, gene = sigGeneStart)
class(sce)
plotGeneCount(crv, counts, gene = sigGeneStart)

?fitGAM
dim(counts)

BiocManager::install("TrajectoryUtils")

BiocManager::version('slingshot')
remove.packages('slingshot')
remove.packages('SingleCellExperiment')
devtools::install('/home/htran/storage/install_software/slingshot/')
devtools::install('/home/htran/storage/install_software/TrajectoryUtils/')
devtools::install('/home/htran/storage/install_software/SingleCellExperiment/')
data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl
pto <- slingshot(rd, cl, start.clus = '1')
rd2 <- cbind(rd[,2] + rnorm(nrow(rd)), -rd[,1] + rnorm(nrow(rd)))
pto.new <- embedCurves(pto, rd2)
pto.new

plot(rd2, col = cl, asp = 1)
lines(SlingshotDataSet(pto.new), lwd = 3)

rd <- cbind(rd, rnorm(nrow(rd)))
pto <- slingshot(rd, cl, start.clus = "1")
sds <- SlingshotDataSet(pto)
# plot3d.SlingshotDataSet(sds, type = 'b')


# https://bustools.github.io/BUS_notebooks_R/slingshot.html#preprocessing
library(viridis)
nc <- 3
pt <- slingPseudotime(pto)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(rd, col = colors, pch = 16, cex = 0.5, main = i)
  lines(pto, lwd = 2, col = 'black', type = 'lineages')
}
