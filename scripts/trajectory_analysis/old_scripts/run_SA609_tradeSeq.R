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
})
datatag <- 'SA609'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
nfeatures_use <- 3000
output_dir <- paste0(save_dir,'tradeseq/')
# norm_sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
dim(norm_sce)

var_genes_df <- data.table::fread(paste0(output_dir, datatag,'_',nfeatures_use,"_hvg_genes.csv")) %>% as.data.frame()
# View(head(var_genes_df))
rownames(norm_sce)[1]
var_genes_df$bc[1]
norm_sce <- norm_sce[var_genes_df$bc[1:300],] # get results first
dim(norm_sce)
# View(as.matrix(counts(sce)[1:5,1:5]))
# sce1 <- sce[, sce$cluster_label==6]
# unique(sce1$treat)
# table(sce$treat, sce$cluster_label)

counts = as.matrix(counts(norm_sce))

meta_genes <- data.frame(ens_gene_id=rownames(norm_sce), gene_symbol=rowData(norm_sce)$Symbol, stringsAsFactors=F)
rownames(meta_genes) <- meta_genes$ens_gene_id
# pseudotime <- data.table::fread(paste0(output_dir, "slingshot_SA535_7_PCA_pseudotime.csv")) %>% as.data.frame()
# # View(head(pseudo))
# pseudotime <- pseudotime %>% tibble::column_to_rownames(var='cell_id')
# 
# cellWeights <- data.table::fread(paste0(output_dir, "slingshot_SA535_7_PCA_cellWeights.csv")) %>% as.data.frame()
# cellWeights <- cellWeights %>% tibble::column_to_rownames(var='cell_id')
# View(head(cellWeights))
crv <- readRDS(paste0(output_dir, "slingshot_pseudotime_SA609_10_PCA_crv.rds"))
set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
# ts_sce <- tradeSeq::fitGAM(counts = as.matrix(counts(norm_sce)), 
#                            pseudotime = pseudotime, cellWeights = cellWeights,
#                            nknots = 6, verbose = FALSE,
#                            parallel=T) 
# saveRDS(ts_sce, paste0(output_dir, "fitGAM_out.rds"))

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
endRes <- readRDS(paste0(output_dir, "endRes_out.rds"))

# patternRes <- patternTest(ts_sce)
# oPat <- order(patternRes$waldStat, decreasing = TRUE)
# print(head(rownames(patternRes)[oPat]))
# saveRDS(patternRes, paste0(output_dir, "patternRes_out.rds"))
patternRes <- readRDS(paste0(output_dir, "patternRes_out.rds"))
# patternRes <- data.table::fread(paste0(output_dir, "patternRes.csv")) %>% as.data.frame()
# class(patternRes)
# View(head(patternRes))


# lin 1: treatment vs 2: drug holiday
colnames(endRes)
endRes12 <- endRes %>%
  dplyr::filter(pvalue_1vs2<0.05)
dim(endRes12)
View(head(endRes12))
endRes12 <- endRes12[order(endRes12$waldStat_1vs2, decreasing = TRUE),]

data.table::fwrite(endRes12, paste0(output_dir,'endRes12_Rx_RxH.csv'))
o <- order(endRes12$waldStat_1vs2, decreasing = TRUE)
plt_ls <- list()
for(i in rep(1:5,1)){
  sigGene <- names(ts_sce)[o[i]]
  gsymb <- meta_genes[sigGene,'gene_symbol']
  p <- plotSmoothers(ts_sce, counts, sigGene) + labs(title=gsymb)
  plt_ls[[gsymb]] <- p
}
output_dir <- paste0(output_dir,'tradeseq/')
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
                           pathway_name=c('custom_pathways'), #'cosmic' or 'cisplatin_resistance', or 'metastasis'
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
patternRes$gene_symb <- meta_genes[rownames(patternRes),'gene_symbol']
patternRes$gene_symb[1]
patternRes <- patternRes %>%
  dplyr::filter(pvalue<0.05)
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