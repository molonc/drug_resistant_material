
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('../../results/v2/outputs/preprocess/sce_annotated/SA906_p15b.rds', '../../data/external/recount2/rse_gene_breast.Rdata', "sce" = '../../results/v2/outputs/preprocess/sce_annotated/SA906_p15b.rds', "gtex" = '../../data/external/recount2/rse_gene_breast.Rdata'),
    output = list('../../results/v2/scratch/infercnv/SA906_p15b.rds', '../../results/v2/scratch/infercnv/SA906_p15b', "rds" = '../../results/v2/scratch/infercnv/SA906_p15b.rds', "infercnv_dir" = '../../results/v2/scratch/infercnv/SA906_p15b'),
    params = list(),
    wildcards = list('SA906_p15b', "sample" = 'SA906_p15b'),
    threads = 1,
    log = list('../../results/v2/logs/infercnv/infercnv/SA906_p15b.log'),
    resources = list(),
    config = list("outdir" = '../../results/v2/outputs', "scratchdir" = '../../results/v2/scratch', "logdir" = '../../results/v2/logs', "workspaces" = list("preprocess" = '../preprocess', "murine" = '../murine', "integrate" = '../integrate', "parse_cnv" = '../parse_cnv', "align_clones" = '../align_clones', "report" = '../report', "infercnv" = '../infercnv', "de" = '../differential_expression'), "scrnaseq_data" = list("sce_library_id" = list("SA532X2XB00147" = 'TENX065', "SA532X6XB00768" = 'TENX066', "SA532X10XB02956" = 'TENX064', "SA609X3XB01584" = 'TENX063', "SA609X10XB02454" = 'TENX062', "SA906_p50b" = 'SCRNA10X_SA_CHIP0080_002', "SA906_p15b" = 'SCRNA10X_SA_CHIP0080_001', "SA906_p11a" = 'SCRNA10X_SA_CHIP0149_001', "SA906_p57a" = 'SCRNA10X_SA_CHIP0149_002', "SA039_p23" = 'SCRNA10X_SA_CHIP0142_001', "SA039_p53" = 'SCRNA10X_SA_CHIP0142_003')), "dlp_data" = list("segments" = list("SA609X3XB01584" = c('../../data/dlp/pdx/cnvd_SA609_EQCJE7191B__868_31_X3.csv'), "SA609X6XB01899" = c('../../data/dlp/pdx/cnvd_SA609_AACIC5777L__868_31_X6.csv'), "SA609X10XB02454" = c('../../data/dlp/pdx/cnvd_SA609_YSOFQ4754R__868_31_X10.csv'), "SA532X2XB00147" = c('../../data/dlp/pdx/cnvd_SA532_VSWXW9069R__459_82_X2.csv'), "SA532X6XB00768" = c('../../data/dlp/pdx/cnvd_SA532_X5_X7_concat.csv'), "SA532X10XB02956" = c('../../data/dlp/pdx/cnvd_SA532_UXOLO8354S__459_82_X9.csv'), "SA906_p30b" = c('../../data/dlp/cl/cnvd_SA906b_IIMJY5615K__419_12_X30.csv', '../../data/dlp/cl/cnvd_SA906b_HRURK1074W__419_12_X35.csv', '../../data/dlp/cl/cnvd_SA906b_YGGKP3732G__419_12_X40.csv'), "SA906_p50b" = c('../../data/dlp/cl/cnvd_SA906b_FMZIH9314Q__419_12_X45.csv', '../../data/dlp/cl/cnvd_SA906b_CPMGU9440J__419_12_X50.csv', '../../data/dlp/cl/cnvd_SA906b_TWFTZ0341T__419_12_X55.csv'), "SA906_p15b" = c('../../data/dlp/cl/cnvd_SA906b_XYTQV2658O__419_12_X20.csv'), "SA906_p11a" = c('../../data/dlp/cl/cnvd_SA906a_AYYUR1554L__471_23_X10.csv', '../../data/dlp/cl/cnvd_SA906a_LVQLJ2710I__471_23_X15.csv', '../../data/dlp/cl/cnvd_SA906a_MZGGB1710S__471_23_X25.csv'), "SA906_p57a" = c('../../data/dlp/cl/cnvd_SA906a_UXOLO8354S__471_23_X30.csv', '../../data/dlp/cl/cnvd_SA906a_FHWHQ2109J__471_23_X40.csv', '../../data/dlp/cl/cnvd_SA906a_HWNCQ2142A__471_23_X50.csv', '../../data/dlp/cl/cnvd_SA906a_FRPKB4575M__471_23_X57.csv'), "SA039_p23" = c('../../data/dlp/cl/cnvd_SA039_AYYUR1554L__459_82_X25.csv', '../../data/dlp/cl/cnvd_SA039_UXOLO8354S__459_82_X30.csv'), "SA039_p53" = c('../../data/dlp/cl/cnvd_SA039_LVQLJ2710I__459_82_X51.csv', '../../data/dlp/cl/cnvd_SA039_FHWHQ2109J__459_82_X60.csv'))), "metadata" = list("sample_metadata" = '../../config/metadata/sample_metadata.yaml'), "conda_path" = '/home/rstudio/miniconda/bin/conda', "conda_environment" = 'r-tensorflow', "qc" = list("max_mito" = 20, "max_ribo" = 60, "min_features" = 1000, "nmads" = 3, "doublet_method" = 'scrublet'), "dimreduce" = list("umap_neighbors" = 30, "umap_min_dist" = 0.5), "integrate" = list("batch_variable" = 'id', "integration_method" = 'scanorama', "rowdata_cols_remove" = c('mean_counts', 'log10_mean_counts', 'n_cells_by_counts', 'pct_dropout_by_counts', 'total_counts', 'log10_total_counts', 'start_position', 'end_position', 'entrezgene', 'hgnc_symbol', 'chromosome_name', 'ensembl_gene_id'), "groups" = list("SA609" = c('SA609X3XB01584', 'SA609X10XB02454'), "SA532" = c('SA532X2XB00147', 'SA532X6XB00768', 'SA532X10XB02956'), "SA906b" = c('SA906_p50b', 'SA906_p15b'), "SA906a" = c('SA906_p11a', 'SA906_p57a'), "SA039" = c('SA039_p23', 'SA039_p53'))), "select_genes" = list("SA609" = list("pct_pure" = 0.4), "SA532" = list("pct_pure" = 0.4), "SA906b" = list("pct_pure" = 0.4), "SA906a" = list("pct_pure" = 0.4), "SA039" = list("pct_pure" = 0.4)), "clonealign" = list("gex_var_quantile" = 0.5, "initial_shrink" = list("SA609" = 0, "SA532" = 10, "SA906b" = 10, "SA906a" = 10, "SA039" = 10), "n_extra_genes" = 100, "clones_to_ignore" = list("SA609" = 'None', "SA532" = 'A', "SA906b" = 'None', "SA906a" = 'None', "SA039" = 'A,B,C')), "report" = list("low_qc_cluster" = list("SA609" = '0', "SA532" = 'none')), "dimred_type" = 'scanorama_TSNE', "differential_expression" = list("clones" = list("SA532X2XB00147" = c('D', 'B'), "SA532X6XB00768" = c('A', 'B'), "SA532X10XB02956" = c('B', 'C', 'D'), "SA609X3XB01584" = c('A', 'B', 'C_D'), "SA609X10XB02454" = c('E', 'F_G', 'I'), "SA906_p50b" = c('B_C', 'D'), "SA906_p15b" = c('F_H_I', 'G'), "SA906_p11a" = c('D_F_H_I', 'C'), "SA906_p57a" = c('B', 'D_F_H_I', 'A'), "SA039_p23" = c('A', 'B'), "SA039_p53" = c('E', 'D')))),
    rule = 'infercnv',
    bench_iteration = as.numeric(NA),
    scriptdir = '/cellassign/fitness-scrna/pipelines/infercnv',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)

######## Original script #########

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







