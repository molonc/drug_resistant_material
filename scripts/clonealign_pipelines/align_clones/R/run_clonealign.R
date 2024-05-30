#' Run CloneAlign
library(argparse)
library(annotables)

parser <- ArgumentParser(description = "Run CloneAlign")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="SingleCellExperiment")
parser$add_argument('--cnv', metavar='FILE', type='character',
                   help="CNV gene file")
parser$add_argument('--clone_prevs', metavar='FILE', type='character',
                    help="Clonal prevalences file")
# parser$add_argument('--gex_var_quantile', type='double',
#                     help="Gene expression quantile", default = 0.5)
parser$add_argument('--n_extra_genes', type='integer',
                    help="Number of extra copy-same genes to add", default = 50)
parser$add_argument('--samples', type = 'character', nargs = '+',
                    help="Sample IDs to subset")
parser$add_argument('--initial_shrink', type = 'double', default = 0,
                    help="Initial shrink")
parser$add_argument('--conda_env', type = 'character',
                    help="Name of conda environment with tensorflow", default = "r-tensorflow")
parser$add_argument('--conda_path', type = 'character',
                    help="Conda path", default = "auto")
parser$add_argument('--clones_to_ignore', type = 'character',
                    help="Clones to ignore (comma separated)", default = "auto")
parser$add_argument('--max_copy_number', type='integer',
                    help="Maximum copy number for input to clonealign", default = 11)
parser$add_argument('--mean_counts', type='double',
                    help="Min mean counts per gene to be included after filtering", default = 100)
## MA: Feb 7, 2024: adding min_counts_per_cell
parser$add_argument('--min_counts_per_cell', type='double',
                    help="Min counts per cell. If sum of counts over all selected genes is <= this threshold, the cell is removed", default = 100)
## MA: 1 Dec 2020: adding data_init_mu
parser$add_argument('--data_init_mu', type='logical',
                    help="Data init mu for clonealign", default = TRUE)
parser$add_argument('--outfname_ca', type = 'character', metavar = 'FILE',
                    help="Output path for fit.")
parser$add_argument('--outfname_fit', type = 'character', metavar = 'FILE',
                    help="Output path for the csv text file fit.")
args <- parser$parse_args()

# Reset snakemake default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, 
                         conda = args$conda_path,
                         required = TRUE)

print(reticulate::py_config())


suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(data.table)
  library(methods)
  library(scran)
  library(yaml)
  
  library(scrna.utils)
  library(scrna.sceutils)
  library(cellassign.utils)
  
  # library(clonealign)
})

set.seed(1345)

devtools::load_all("/cellassign/clonealign/")


sce <- readRDS(args$sce)

# This is removing the rows with NA IDs, also done in identify-murine.R
sce <- sce[which(!is.na(rowData(sce)$ID)),]
rownames(sce) <- rowData(sce)$ID

raw_cnvs <- fread(args$cnv)
clone_prevalences <- fread(args$clone_prevs)
#gex_var_quantile <- args$gex_var_quantile  # not used



## Remove clones to ignore
if(args$clones_to_ignore != "None") {
  clones_to_ignore <- strsplit(args$clones_to_ignore, ",")[[1]]
  raw_cnvs <- raw_cnvs[!(raw_cnvs$cluster %in% clones_to_ignore),]
}

## Get list of timepoints to evaluate

sce <- sce %>%
  scater::filter(id %in% unlist(args$samples))
timepoints <- unique(sce$timepoint)

colnames(sce) <- paste0(sce$id, "_", sce$Barcode)


## Important -- select autosomes only

autosomal_genes <- annotables::grch37 %>% 
  dplyr::select(ensembl_gene_id = ensgene, chr) %>% 
  dplyr::filter(chr %in% as.character(1:22)) %>% 
  .$ensembl_gene_id


#save.image(paste0("/cellassign/fitness-scrna/results/dumps/_clonealign-dump-", sce$id[1], ".Rdata"))
# clone_prevalences <- clone_prevalences %>%
#   dplyr::filter(time %in% timepoints)

# if(nrow(clone_prevalences) == 0) {
#   stop("No matching clone prevalences!")
# }

# present_clones <- unique(clone_prevalences$cluster)

## Filter CNV data
cnv <- dplyr::filter(raw_cnvs, use_gene) %>%
  dplyr::rename(clone = cluster,
                copy_number=median_cnmode) %>% 
  dplyr::select(ensembl_gene_id, clone, copy_number) %>% 
  # dplyr::filter(clone %in% present_clones) %>%
  spread(clone, copy_number)

cnv_mat <- cnv %>%
  as.data.frame %>%
  column_to_rownames("ensembl_gene_id") %>%
  as.matrix

print("Dimensions of cnv_mat")
print(dim(cnv_mat))
numclones <- dim(cnv_mat)[2]
myclone <- paste(colnames(cnv_mat), collapse='_')

fixsce <- function (sce) {
  # some of them don;t have these, removing them for now
  # rowData(sce)$chromosome_name <- rowData(sce)$chr
  # rowData(sce)$start_position <- rowData(sce)$start
  # rowData(sce)$end_position <- rowData(sce)$end
  drops <- c("chr","start","end","strand","symbol","ensgene","entrez","biotype","description")
  rdcopy <- rowData(sce)
  rowData(sce) <- rdcopy[ , !(names(rdcopy) %in% drops)]
  return(sce)
}



common_genes <- intersect(intersect(rownames(cnv_mat), rownames(sce)), autosomal_genes)
sce <- fixsce(sce)
sce <- sce[common_genes,]
cnv_mat <- cnv_mat[common_genes,]
# MA: this type of subsetting doesn't work if there is only one column
#cnv_mat <- subset(cnv_mat, rownames(cnv_mat) %in% common_genes)
#cnv_mat <- as.matrix(cnv_mat[common_genes,])

# MA: If there is only one clone, write the result without running clonealign
# MA: or if all CNs are the same for all genes for all clones, do the same
if (numclones == 1) {
  print("WARNING: Only one clone, bypassing clonealign!!!")
  sce_filtered <- sce
  
  clone_fit <- data.frame(
    Sample=sce_filtered$Sample,
    Barcode=sce_filtered$Barcode,
    cell_id=colnames(sce_filtered),
    id=sce_filtered$id,
    clone=myclone,
    stringsAsFactors = FALSE
  )
  # Not writing ca file, just the fit.csv file
  # Actually I should also write ca file, otherwise the pipeline gives a MissingOutputException
  ca <- NULL
  saveRDS(ca, args$outfname_ca)
  
  write_csv(clone_fit, args$outfname_fit)
  print("Completed")
  
} else if (sum(rowSds(cnv_mat)) <= 1) {   ## for SA609X7XB03573, there is only 1 difference
  print("WARNING: the cnv matrix has 0 or 1 CN differences, bypassing clonealign!!!")
  sce_filtered <- sce
  
  clone_fit <- data.frame(
    Sample=sce_filtered$Sample,
    Barcode=sce_filtered$Barcode,
    cell_id=colnames(sce_filtered),
    id=sce_filtered$id,
    clone=myclone,
    stringsAsFactors = FALSE
  )
  # Not writing ca file, just the fit.csv file
  # Actually I should also write ca file, otherwise the pipeline gives a MissingOutputException
  ca <- NULL
  saveRDS(ca, args$outfname_ca)
  
  write_csv(clone_fit, args$outfname_fit)
  print("Completed")  
} else {
  processed_data <- preprocess_for_clonealign(sce, 
                                              cnv_mat,
                                              mean_counts_per_gene = args$mean_counts,
                                              min_counts_per_cell = args$min_counts_per_cell) 
  ## For SA604X6XB01979 we need fewer min_counts_per_cell. Normally we use 100
  
  
  cnv_mat_filtered <- processed_data$copy_number_data
  
  dist_matrix <- as.matrix(dist(t(cnv_mat_filtered), method = "manhattan"))/nrow(cnv_mat_filtered)
    
  clone_remaps <- reshape2::melt(dist_matrix) %>% 
    dplyr::filter(value < 0.005) %>% 
    dplyr::group_by(Var1) %>% 
    dplyr::summarise(merged_id = paste(sort(unique(as.character(Var2))), collapse = "_"))
  
  cnv_mat_remapped <- cnv_mat_filtered[,clone_remaps$Var1[!duplicated(clone_remaps$merged_id)],drop=FALSE]
  
  colnames(cnv_mat_remapped) <- colnames(cnv_mat_remapped) %>%
    plyr::mapvalues(from = clone_remaps$Var1, clone_remaps$merged_id)
  
  
  processed_data <- preprocess_for_clonealign((processed_data$gene_expression_data), 
                                              cnv_mat_remapped,
                                              mean_counts_per_gene = args$mean_counts,
                                              min_counts_per_cell = args$min_counts_per_cell)
  
  
  # save.image(paste0("/cellassign/fitness-scrna/results/dumps/clonealign-dump-", sce$id[1], ".Rdata"))
  
  print(paste("Final dimensions expression data", dim(processed_data$gene_expression_data)))
  
  ca <- run_clonealign(processed_data$gene_expression_data, 
                       processed_data$copy_number_data,
                       #initial_shrinks = c(2,5,10),     #### On Feb 8, 2024 it was like this. some of the clonealign files maye have run c(2,5,10) instead of args$initial_shrink
                       initial_shrinks = 0,
                       n_repeats = 3,
                       mc_samples = 1,
                       learning_rate = 0.07,
                       max_iter = 5e2,
                       ## for SA609X5XB03231, we need data_init_mu = FALSE. Needs to be added to the pipeline, but didn't manage to.
                       data_init_mu = TRUE,
                       saturation_threshold = args$max_copy_number,
                       clone_call_probability = 0.9,
                       seed = 12345)
  
  sce_filtered <- sce[, rownames(processed_data$gene_expression_data)]
  
  clone_fit <- data.frame(
    Sample=sce_filtered$Sample,
    Barcode=sce_filtered$Barcode,
    cell_id=rownames(processed_data$gene_expression_data),
    id=sce_filtered$id,
    clone=ca$clone,
    stringsAsFactors = FALSE
  )
  
  ca$clone_fit <- clone_fit
  ca$cnv_mat <- cnv_mat_remapped
  
  # Save results
  saveRDS(ca, args$outfname_ca)
  
  # Now saving the csv file
  cat("Saving the csv file\n")
  ca <- clonealign:::recompute_clone_assignment(ca, 0.5)
  
  df <- ca$clone_fit
  df$clone <- ca$clone
  
  write_csv(df, args$outfname_fit)  
  
  cat("Completed.\n")  
  
}  

  
