devtools::install_github("stephenturner/annotables")

library(infercnv)
library(SingleCellExperiment)
library(annotables)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

set.seed(123L)

ca_x3 <- readRDS("../../results/v1/outputs/align_clones/clonealign_fit/SA609X3XB01584.rds")
ca_x10 <- readRDS("../../results/v1/outputs/align_clones/clonealign_fit/SA609X10XB02454.rds")

load("../../data/external/recount2/rse_gene_breast.Rdata")
sce <- readRDS("../../results/v1/outputs/integrate/batch_correct/SA609.rds")
rownames(sce) <- rowData(sce)$ID

sce_x10 <- sce[, sce$id == "SA609X10XB02454"]
ca <- ca_x10

sce_tmp <- sce_x10[, sample(ncol(sce_x10), 500)]

sce_tmp <- sce_tmp[rowSums(counts(sce_tmp)) > 50,]

colnames(sce_tmp) <- sce_tmp$Barcode

gtex_counts <- assay(rse_gene, 'counts')

unique_ids <- gsub("\\.[-0-9]+", "", rownames(gtex_counts))
rownames(gtex_counts) <- unique_ids

common_genes <- intersect(rownames(sce_tmp), rownames(gtex_counts))

counts <- cbind(
  as.matrix(counts(sce_tmp)[common_genes,]),
  gtex_counts[common_genes,]
)

source <- c(
  rep("SA609", ncol(sce_tmp)),
  rep("GTEX", ncol(gtex_counts))
)

go <- grch37 %>% 
  dplyr::select(ensgene, chr, start, stop = end) 

go <- go[!duplicated(go),]

write_tsv(go, "gene_order.tsv", col_names = FALSE)

cell_names <-  colnames(counts)
tibble(cell_names, source) %>% 
  write_tsv("annotations.tsv", col_names = FALSE)



infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
                                   annotations_file="annotations.tsv",
                                    delim="\t",
                                    gene_order_file="gene_order.tsv",
                                    ref_group_names=c("GTEX")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="output", 
                             cluster_by_groups=TRUE, 
                             window_length = 51,
                             denoise=TRUE,
                             scale_data=TRUE,
                             HMM=FALSE)

norm_expr <- infercnv_obj@expr.data

quants <- quantile(as.vector(norm_expr), p = c(0.01, 0.99))

norm_expr[norm_expr < quants[1]] <- quants[1]
norm_expr[norm_expr > quants[2]] <- quants[2]

norm_expr <- norm_expr[,sce_tmp$Barcode]

clones <- ca$clone_fit$clone
names(clones) <- ca$clone_fit$Barcode
clones <- clones[colnames(norm_expr)]

df_clone <- data.frame(clone = clones)
rownames(df_clone) <- names(clones)

go <- infercnv_obj@gene_order

go$chr <- factor(go$chr, c(as.character(1:23)))

go <- rownames_to_column(go,'gene') %>% 
  dplyr::filter(!is.na(chr)) %>% 
  dplyr::arrange(chr, start) %>% 
  dplyr::select(gene, chr)

rownames(go) <- go$gene
go$gene <- NULL

norm_expr <- norm_expr[rownames(go),]

unique_clones <- unique(clones)
clone_cols <- brewer.pal(length(unique_clones), "Set1")
names(clone_cols) <- unique_clones

ann_colors = list(
  clone = clone_cols
)

norm_expr <- norm_expr[, infercnv_obj@tumor_subclusters$subclusters[[1]][[1]]]

png('test.png', width = 970, height = 800)
pheatmap(t(norm_exprs[rownames(norm_exprs) %in% ca$retained_genes,]), 
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         show_colnames = FALSE, 
         show_rownames = FALSE,
         annotation_row = df_clone,
         annotation_col = go, 
         annotation_colors = ann_colors)
dev.off()


png('test.png', width = 970, height = 800)
pheatmap(t(norm_exprs[rownames(norm_exprs) %in% ca$retained_genes,]), 
         cluster_cols = FALSE, 
         show_colnames = FALSE, 
         show_rownames = FALSE,
         annotation_row = df_clone,
         annotation_col = go, 
         annotation_colors = ann_colors)
dev.off()






