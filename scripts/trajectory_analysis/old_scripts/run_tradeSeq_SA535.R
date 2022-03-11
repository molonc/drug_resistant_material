# tradeSeq 

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot, quietly = TRUE)
  # library(mclust, quietly = TRUE)
  library(tradeSeq)
  library(grDevices)
  library(RColorBrewer)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(clusterExperiment)
})
# ?tradeSeq::fitGAM


datatag <- 'SA535'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
output_dir <- paste0(save_dir,'tradeseq/')
nfeatures_use <- 3000
norm_sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
dim(norm_sce)

var_genes_df <- data.table::fread(paste0(save_dir, "SA535_3000_hvg_genes.csv")) %>% as.data.frame()
# View(head(var_genes_df))
rownames(norm_sce)[1]
var_genes_df$bc[1]
norm_sce <- norm_sce[var_genes_df$bc[1:500],] # get results first
counts <- as.matrix(counts(norm_sce))
dim(counts)
meta_genes <- data.frame(ens_gene_id=rownames(norm_sce), gene_symbol=rowData(norm_sce)$Symbol, stringsAsFactors=F)
rownames(meta_genes) <- meta_genes$ens_gene_id

# View(as.matrix(counts(sce)[1:5,1:5]))
# sce1 <- sce[, sce$cluster_label==6]
# unique(sce1$treat)
# table(sce$treat, sce$cluster_label)



# pseudotime <- data.table::fread(paste0(output_dir, "slingshot_SA535_7_PCA_pseudotime.csv")) %>% as.data.frame()
# # View(head(pseudo))
# pseudotime <- pseudotime %>% tibble::column_to_rownames(var='cell_id')
# 
# cellWeights <- data.table::fread(paste0(output_dir, "slingshot_SA535_7_PCA_cellWeights.csv")) %>% as.data.frame()
# cellWeights <- cellWeights %>% tibble::column_to_rownames(var='cell_id')
# View(head(cellWeights))
crv <- readRDS(paste0(save_dir, "slingshot_pseudotime_SA535_7_PCA_crv.rds"))
set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
ts_sce <- readRDS(paste0(output_dir, "fitGAM_out.rds"))

# assoRes <- associationTest(ts_sce)
# saveRDS(assoRes, paste0(output_dir, "assoRes_out.rds"))
assoRes <- readRDS(paste0(output_dir, "assoRes_out.rds"))

# startRes <- startVsEndTest(ts_sce, lineages=TRUE)
# saveRDS(startRes, paste0(output_dir, "startRes_out.rds"))
# print(head(startRes))

startRes <- readRDS(paste0(output_dir, "startRes_out.rds"))

# endRes <- diffEndTest(ts_sce, pairwise=TRUE)
# print(head(endRes))
# saveRDS(endRes, paste0(output_dir, "endRes_out.rds"))


# patternRes <- patternTest(ts_sce)
# oPat <- order(patternRes$waldStat, decreasing = TRUE)
# print(head(rownames(patternRes)[oPat]))
# saveRDS(patternRes, paste0(output_dir, "patternRes_out.rds"))

class(patternRes)
View(head(endRes))


# lin 1: treatment vs 2: drug holiday
endRes <- readRDS(paste0(output_dir, "endRes_out.rds"))

colnames(endRes)
endRes12ss <- endRes %>%
  dplyr::select(waldStat_1vs2, pvalue_1vs2, logFC1_2)

endRes12 <- endRes %>%
  dplyr::filter(pvalue_1vs2<0.05)
endRes12$gene_symb <- meta_genes[rownames(endRes12),'gene_symbol']

dim(endRes12) # 263 genes in total of 500 hvg genes
View(head(endRes12ss))
length(intersect(endRes12$gene_symb, patternRes$gene_symb))
endRes12 <- endRes12[order(endRes12$waldStat_1vs2, decreasing = TRUE),]
endRes12$gene_symb <- meta_genes[rownames(endRes12),'gene_symbol']
data.table::fwrite(endRes12, paste0(output_dir,'endRes12_Rx_RxH.csv'))
o <- order(endRes12$waldStat_1vs2, decreasing = TRUE)
plt_ls <- list()
for(i in rep(1:5,1)){
  sigGene <- names(ts_sce)[o[i]]
  gsymb <- meta_genes[sigGene,'gene_symbol']
  p <- plotSmoothers(ts_sce, counts, sigGene) + labs(title=gsymb)
  plt_ls[[gsymb]] <- p
}
# output_dir <- paste0(output_dir,'tradeseq/')
# dir.create(output_dir)
p_total <- cowplot::plot_grid(plotlist = plt_ls, ncol = 5, align='h')
png(paste0(output_dir,"tradeseq_",datatag,"_endRes12_Rx_RxH_top5genes.png"), 
    height = 2*250, width=2*1300,res = 2*72)
print(p_total)
dev.off()
rownames(endRes12)[1]
endRes12$logFC1_2
endRes12$gene_symb <- meta_genes[rownames(endRes12),'gene_symbol']
deg_df <- endRes12 %>%
  dplyr::select(logFC1_2, gene_symb)%>%
  dplyr::rename(logFC=logFC1_2)

dim(deg_df)
# TO DO: what is the pathway that related to RxH genes? 144 genes
gsea_out <- get_custom_pathway_results(deg_df,                      # named vector of statistical significance 
                                       desc='endRes12', base_name = datatag,  
                                       pathway_name=c('hallmark'), #'cosmic' or 'cisplatin_resistance', or 'metastasis'
                                       groups_use=c("Rx","RxH"),    # vector of 2 elements, 2 group name used for DE analysis
                                       output_dir,
                                       reference_genes_set=NULL,
                                       n_top=20,
                                       ref_dif=NULL)

gsea_out$pathway
class(gsea_out)
gsea_out <- gsea_out %>%
  dplyr::filter(padj < 0.05)
# "HALLMARK_G2M_CHECKPOINT"

patternRes <- readRDS(paste0(output_dir, "patternRes_out.rds"))
patternRes$gene_symb <- meta_genes[rownames(patternRes),'gene_symbol']
patternRes$gene_symb[1]
patternRes <- patternRes %>%
  dplyr::filter(pvalue<0.05)
patternRes <- patternRes %>%
  dplyr::filter(pvalue<0.05 & gene_symb=='MYC')
patternRes$pvalue
dim(patternRes)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
patternRes <- patternRes[order(patternRes$waldStat, decreasing = TRUE),]
data.table::fwrite(patternRes, paste0(output_dir,'patternRes.csv'))
patternRes$waldStat
deg_df <- patternRes %>%
  dplyr::select(waldStat, gene_symb)%>%
  dplyr::rename(logFC=waldStat)

dim(deg_df)
# TO DO: what is the pathway that related to RxH genes? 144 genes
gsea_out <- get_custom_pathway_results(deg_df,                      # named vector of statistical significance 
                                       desc='patternRes', base_name = datatag,  
                                       pathway_name=c('hallmark'), #'cosmic' or 'cisplatin_resistance', or 'metastasis'
                                       groups_use=c("Rx","RxH"),    # vector of 2 elements, 2 group name used for DE analysis
                                       output_dir,
                                       reference_genes_set=NULL,
                                       n_top=20,
                                       ref_dif=NULL)
gsea_out$pathway[gsea_out$padj<0.05]
head(rownames(patternRes)[oPat])
meta_genes[rownames(patternRes)[oPat][1:10],'gene_symbol']
plt_ls <- list()
for(i in rep(1:10,1)){
  sg <- rownames(patternRes)[oPat][i]
  gsymb <- meta_genes[sg,'gene_symbol']
  p <- plotSmoothers(ts_sce, counts, gene = sg)+ labs(title=gsymb)
  plt_ls[[gsymb]] <- p
  print(gsymb)
}
p_total <- cowplot::plot_grid(plotlist = plt_ls, ncol = 5, align='hv')
png(paste0(output_dir,"tradeseq_",datatag,"_patternRes_Rx_RxH_UnRx_top10genes.png"), 
    height = 2*300, width=2*400,res = 2*72)
print(p)
dev.off()

png(paste0(output_dir,"tradeseq_",datatag,"_patternRes_MYC_gene.png"), 
    height = 2*500, width=2*1300,res = 2*72)
print(p_total)
dev.off()

# 1: Rx, 2: RxH 3: UnRx
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(ts_sce)[oStart[1]]
p <- plotSmoothers(ts_sce, counts, gene = sigGeneStart)
p


earlyDERes <- earlyDETest(ts_sce, knots = c(1, 2))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
earlyDERes <- earlyDERes[order(earlyDERes$waldStat, decreasing = TRUE),]
earlyDERes$gene_symb <- meta_genes[rownames(earlyDERes),'gene_symbol']
earlyDERes <- earlyDERes %>%
  dplyr::filter(pvalue<0.05)
dim(earlyDERes)
data.table::fwrite(earlyDERes, paste0(output_dir,'earlyDERes.csv'))

head(rownames(earlyDERes)[oEarly])
plt_ls <- list()
for(i in rep(1:10,1)){
  sg <- rownames(earlyDERes)[oEarly][i]
  gsymb <- meta_genes[sg,'gene_symbol']
  p <- plotSmoothers(ts_sce, counts, gene = sg)+ labs(title=gsymb)
  plt_ls[[gsymb]] <- p
  print(gsymb)
}
p_total <- cowplot::plot_grid(plotlist = plt_ls, ncol = 5, align='hv')
png(paste0(output_dir,"tradeseq_",datatag,"_earlyDrive_Rx_RxH_UnRx_top10genes.png"), 
    height = 2*500, width=2*1300,res = 2*72)
print(p_total)
dev.off()
norm_sce
pcs <- reducedDims(norm_sce[rownames(ts_sce),colnames(ts_sce)])[['PCA']]
reducedDims(ts_sce) <- SimpleList(PCA = as.matrix(pcs))
dim(reducedDims(ts_sce)[['PCA']])
dim(reducedDims(norm_sce)[['PCA']])
dim(ts_sce)
library(clusterExperiment)
nPointsClus <- 20
clusPat <- tradeSeq::clusterExpressionPatterns(ts_sce, nPoints = nPointsClus,
                                               ncores = 5, nReducedDims = 30,reduceMethod = "PCA",
                                               genes = rownames(ts_sce)[1:200], minSizes = 10)
?tradeSeq::clusterExpressionPatterns
clusterLabels <- primaryCluster(clusPat$rsec)
length(clusterLabels)
summary(as.factor(clusterLabels))
cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1] # remove unclustered genes
cUniq
for (xx in cUniq[1:length(cUniq)]) {
  cId <- which(clusterLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 4),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:3, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + #guides(color = FALSE) +
    scale_color_manual(values = c("orange", "darkseagreen3","blue","yellow"),
                       breaks = c("0", "1","2","3"))  
  print(p)
}
dim(patternRes)
# Using umap clustering 
# View(head(patternRes))
preprocess_mat <- as.matrix(logcounts(norm_sce[rownames(patternRes),colnames(ts_sce)]))
dim(preprocess_mat)
gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, 
                                    resolution=0.05, max_components = 2)
summary(as.factor(gene_module_df$module))
gene_module_df$id[1]
colnames(meta_genes)
gene_module_df$gene_symb <- meta_genes[gene_module_df$id,'gene_symbol']
data.table::fwrite(gene_module_df, paste0(output_dir,'gene_modules_patternTest_',dim(patternRes)[1],'_genes.csv'))

norm_sce$treatmentSt <- norm_sce$treat
cell_group_df <- tibble::tibble(cell=colnames(preprocess_mat), 
                                cell_group=colData(norm_sce)[colnames(preprocess_mat),'treatmentSt']) #[colnames(cds)] #partitions(cds)
agg_mat <- aggregate_gene_expression_v3(preprocess_mat, gene_module_df, cell_group_df, max_agg_value = 4, min_agg_value = -4)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
# colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
print(dim(agg_mat))
# library(pheatmap)
p1 <- pheatmap::pheatmap(agg_mat, cluster_rows=T, cluster_cols=TRUE,
                         scale="column", clustering_method="ward.D2",
                         fontsize=16)
png(paste0(output_dir, "genes_modules_",datatag,".png"), height =2*600, width= 2*800,res = 2*72)
print(p1)
dev.off()


# data frame: lineage - pseudotime, cells id, genes id, expression
# color by lineages
colnames(patternRes)
patternRes1 <- patternRes
patternRes1$id <- rownames(patternRes1)
patternRes1$gene_symb <- NULL

obs_module <- 6
gene_module_df$module
obs_genes_df <- gene_module_df %>%
  dplyr::filter(module==obs_module)
dim(obs_genes_df)
# in case of drug holiday
obs_genes_df <- obs_genes_df %>%
  dplyr::filter(id %in% rownames(endRes12))

obs_genes_df <- obs_genes_df %>% dplyr::inner_join(patternRes1, by=c('id'))
obs_genes_df <- obs_genes_df[order(obs_genes_df$waldStat, decreasing = T),]
obs_genes_df <- obs_genes_df[1:5,]
View(obs_genes_df)
genes_use <- obs_genes_df$id

plt_ls <- list()
for(sg in genes_use){
  gsymb <- meta_genes[sg,'gene_symbol']
  p <- tradeSeq::plotSmoothers(ts_sce, counts, gene = sg)+ labs(title=gsymb)
  plt_ls[[gsymb]] <- p
  print(gsymb)
}
p_total <- cowplot::plot_grid(plotlist = plt_ls, ncol = 5, align='hv')
png(paste0(output_dir,"tradeseq_",datatag,"_patternRes_top5genes_module_",obs_module,".png"), 
    height = 2*250, width=2*1300,res = 2*72)
print(p_total)
dev.off()



