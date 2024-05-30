

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(clonealign)
})


infercnv_obj <- readRDS(snakemake@input$rds)
ca <- readRDS(snakemake@input$clonealign_fit)

ca <- clonealign:::recompute_clone_assignment(ca, 0.51)
ca$clone_fit$clone <- ca$clone

norm_expr <- infercnv_obj@expr.data

quants <- quantile(as.vector(norm_expr), p = c(0.01, 0.99))

norm_expr[norm_expr < quants[1]] <- quants[1]
norm_expr[norm_expr > quants[2]] <- quants[2]

common_cells <- intersect(colnames(norm_expr), ca$clone_fit$Barcode)

ca$clone_fit <- ca$clone_fit[match(common_cells, ca$clone_fit$Barcode),]

norm_expr <- norm_expr[,ca$clone_fit$Barcode]

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

# norm_expr <- norm_expr[, unliinfercnv_obj@tumor_subclusters$subclusters[[1]] ]

png(snakemake@output$png_all_genes, width = 800, height = 2 * ncol(norm_expr))
pheatmap(t(norm_expr), 
         cluster_cols = FALSE, 
         cluster_rows = TRUE,
         show_colnames = FALSE, 
         show_rownames = FALSE,
         annotation_row = df_clone,
         annotation_col = go, 
         annotation_colors = ann_colors)
dev.off()


png(snakemake@output$png_clonealign_only, width = 1200, height =  2 * ncol(norm_expr))
pheatmap(t(norm_expr[rownames(norm_expr) %in% ca$retained_genes,]), 
         cluster_cols = FALSE, 
         show_colnames = FALSE, 
         show_rownames = FALSE,
         annotation_row = df_clone,
         annotation_col = go, 
         annotation_colors = ann_colors)
dev.off()

