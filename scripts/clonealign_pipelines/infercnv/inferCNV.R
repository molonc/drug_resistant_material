
suppressPackageStartupMessages({
  library(infercnv)
  library(SingleCellExperiment)
  library(annotables)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(123L)

sample <- snakemake@wildcards$sample

load(snakemake@input$gtex)
sce <- readRDS(snakemake@input$sce)
infercnv_dir <- snakemake@output$infercnv_dir
output_rds <- snakemake@output$rds

#  delete inferCNV files to overwrite
unlink(infercnv_dir, recursive = TRUE)

rownames(sce) <- rowData(sce)$ID

n_cells_to_sample <- min(ncol(sce), 500)
sce <- sce[, sample(ncol(sce), n_cells_to_sample)]

sce <- sce[rowSums(counts(sce)) > 0,]

colnames(sce) <- sce$Barcode

gtex_counts <- assay(rse_gene, 'counts')

unique_ids <- gsub("\\.[-0-9]+", "", rownames(gtex_counts))
rownames(gtex_counts) <- unique_ids

common_genes <- intersect(rownames(sce), rownames(gtex_counts))

counts <- cbind(
  as.matrix(counts(sce)[common_genes,]),
  gtex_counts[common_genes,]
)

source <- c(
  rep("SA609", ncol(sce)),
  rep("GTEX", ncol(gtex_counts))
)

go <- grch37 %>% 
  dplyr::select(ensgene, chr, start, stop = end) 

go <- go[!duplicated(go),]

write_tsv(go, "gene_order.tsv", col_names = FALSE)

cell_names <-  colnames(counts)
tibble(cell_names, source) %>% 
  write_tsv("annotations.tsv", col_names = FALSE)



infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=counts,
                                   annotations_file="annotations.tsv",
                                    delim="\t",
                                    gene_order_file="gene_order.tsv",
                                    ref_group_names=c("GTEX")) 

infercnv_obj <- infercnv::run(infercnv_obj,
                             cutoff=.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=infercnv_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             scale_data = TRUE,
                             HMM=FALSE)

saveRDS(infercnv_obj, output_rds)







