
suppressPackageStartupMessages({
  library(tidyverse)
  library(scater)
  library(data.table)
  library(argparse)
  library(broom)
  library(clonealign)
})

set.seed(1345)


parser <- ArgumentParser(description = "Run CloneAlign")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="SingleCellExperiment")
parser$add_argument('--cnv', metavar='FILE', type='character',
                    help="CNV gene file")
parser$add_argument('--clonealign_fit', metavar='FILE', type='character',
                    help="Clonealign fit RDS")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for fit.")
parser$add_argument('--sample', type = 'character', metavar = 'FILE',
                    help="Sample name.")
args <- parser$parse_args()

sce <- readRDS(args$sce)
ca <- readRDS(args$clonealign_fit)
sample <- args$sample

sce <- sce %>%
  scater::filter(id %in% unlist(args$sample))

save.image(paste0("/cellassign/fitness-scrna/results/dumps/eval-", args$sample, ".rds"))

raw_cnvs <- fread(args$cnv)

rownames(sce) <- rowData(sce)$ID
colnames(sce) <- sce$Barcode

ca <- clonealign:::recompute_clone_assignment(ca, 0.5)

## Get CNV matrix with genes _not_ used
cnv <- dplyr::filter(raw_cnvs, !use_gene) %>%
  dplyr::rename(clone = cluster,
                copy_number=median_cnmode) %>% 
  dplyr::select(ensembl_gene_id, clone, copy_number) %>% 
  spread(clone, copy_number)

cnv_mat <- cnv %>%
  as.data.frame %>%
  column_to_rownames("ensembl_gene_id") %>%
  as.matrix


## Inferred clones
inferred_clones <- unique(ca$clone)
inferred_clones <- setdiff(inferred_clones, "unassigned")

## Need to break out cnv_mat if clones have been collapsed
collapsed_clones <- grepl("_", inferred_clones)
if(any(collapsed_clones)) {
  for(i in which(collapsed_clones)) {
    cclone <- inferred_clones[i]
    uclones <- unlist(strsplit(cclone, "_"))
    new_clone <- rowMedians(cnv_mat[, uclones])
    cnv_mat <- cbind(cnv_mat, new_clone)
    colnames(cnv_mat)[ncol(cnv_mat)] <- cclone
  }
}

## 

cnv_mat <- cnv_mat[, inferred_clones]
cnv_mat <- cnv_mat[matrixStats::rowVars(cnv_mat) > 0,]
cnv_mat <- cnv_mat[matrixStats::rowMaxs(cnv_mat) < 6,]

sce <- sce[,ca$clone_fit$Barcode]

sce <- sce[rowSums(as.matrix(counts(sce))) > 100, ]

common_genes <- intersect(rownames(sce), rownames(cnv_mat))

sce <- sce[common_genes,]
cnv_mat <- cnv_mat[common_genes,]

clones <- ca$clone

assigned_cells <- clones != "unassigned"

sce <- sce[, assigned_cells]
clones <- clones[assigned_cells]

logcs <- logcounts(sce)
cnv_mat_full <- cnv_mat[, clones]

test_estimates <- lapply(seq_len(nrow(sce)), function(i) {
  lc <- logcs[i,]
  cnv_dist <- cnv_mat_full[i,]
  tidy(lm(lc ~ cnv_dist))[2,]
}) %>% 
  bind_rows()


cnv_mat_full <- cnv_mat[, sample(clones)]
null_estimates <- lapply(seq_len(nrow(sce)), function(i) {
  lc <- logcs[i,]
  cnv_dist <- cnv_mat_full[i,]
  tidy(lm(lc ~ cnv_dist))[2,]
}) %>% 
  bind_rows()

df <- bind_rows(
  dplyr::mutate(test_estimates, dist = "observed"),
  dplyr::mutate(null_estimates, dist = "null")
)

tt <- t.test(test_estimates$estimate, null_estimates$estimate)

round2 <- function(x) format(round(x, 2), nsmall = 2)

ggplot(df, aes(x = dist, y = estimate)) +
  geom_boxplot(outlier.shape = NA, size = .4) +
  labs(x = "Distribution", y = "Coefficient expression ~ copy number",
       title = sample,
       subtitle = paste0("Genes not used by clonealign, p = ", round2(tt$p.value))) +
  theme_bw() +
  ylim(-.2, .2)

ggsave(args$outfname, width = 4, height = 4)




