
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scran)
  library(scater)
  library(pheatmap)
  library(clonealign)
  library(annotables)
  library(tidyverse)
  library(fgsea)
})


# sample <- "SA906_p50b"
sample <- snakemake@params[['sample']]

# interesting_clones <- c("B_C", "D")
interesting_clones <- snakemake@params[['interesting_clones']]

ca <- readRDS(snakemake@input[['clonealign_fit']])

sce <- readRDS(snakemake@input[['sce']])

# save.image(paste0("/cellassign/fitness-scrna/results/dumps/", sample, '.rds'))

sce <- sce[ , match(ca$clone_fit$Barcode, sce$Barcode)]

ca <- clonealign:::recompute_clone_assignment(ca, 0.5)

sce$clone <- ca$clone

sce <- sce[, sce$clone %in% interesting_clones]

## Make sure all clones are represented

if(!all.equal(sort(interesting_clones), sort(unique(sce$clone)))) {
  stop("Not all interesting_clones present after filtering")
}

# rownames(sce) <- scater::uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
rownames(sce) <- rowData(sce)$ID

bad_genes <- grepl("^MT-|^RP[L|S]", rowData(sce)$Symbol)
sce <- sce[!bad_genes,]

sce_de <- sce[rowSums(as.matrix(counts(sce))) > 50, ]

## Load and parse CNV data

cnv_data <- read_tsv(snakemake@input[['clone_cn']])
cnv_data <- dplyr::filter(cnv_data, cluster %in% interesting_clones)

## GSEA

load(snakemake@input[['gs_hallmark']])
load(snakemake@input[['gs_go']])
load(snakemake@input[['gs_oncogenic']])

get_de_pathway_results <- function(clone) {


  fm <- findMarkers(sce_de, sce$clone == clone)


  lfcs <- as.data.frame(fm[[ 2 ]]) %>% 
    rownames_to_column('ensgene')

  lfcs <- dplyr::select(grch37, ensgene, entrez, symbol) %>% 
    inner_join(lfcs) %>% 
    rename(logFC = logFC.FALSE)
  
  logfc <- lfcs[[ 'logFC' ]]
  
  names(logfc) <- lfcs$entrez
  
  gsea_hallmark <- fgsea(pathways=Hs.H, stats=logfc, nperm=10000)
  gsea_go <- fgsea(pathways=Hs.c5, stats=logfc, nperm=10000)
  gsea_oncogenic <- fgsea(pathways=Hs.c6, stats=logfc, nperm=10000)
  
  list(
    clone = clone,
    interesting_clones = interesting_clones,
    de_results = lfcs,
    gsea_hallmark = gsea_hallmark,
    gsea_go = gsea_go,
    gsea_oncogenic = gsea_oncogenic,
    cnv_data = cnv_data
  )

}

results <- lapply(interesting_clones, get_de_pathway_results)


saveRDS(results, snakemake@output[['rds']])
