source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation_v1/')
output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/')

datatag <- 'SA609'
sce_scMerge <- readRDS(paste0(output_dir,'SA609_scMerge_correction_without_cosine_v2.rds'))  #SA609_scMerge_correction.rds dim(sce_transform)
exprs <- 'scMerge_fast'
sce_scMerge$treatmentSt <- get_treatment_status(sce_scMerge$series)



metadata <- as.data.frame(colData(sce_scMerge))
metadata$cell_id <- colnames(sce_scMerge)
metadata <- metadata %>%
  dplyr::select(cell_id, treatmentSt, clone, batch, series)

# ext_cells <- downsampling_data(metadata, col_use='treatmentSt', downsample_ratio=0.3, thres_small_clone=600)
# sce_scMerge <- sce_scMerge[,ext_cells]
# print(dim(sce_scMerge))
# metadata <- metadata %>%
#   dplyr::filter(cell_id %in% ext_cells)
# print(dim(metadata))
assay(sce_scMerge, exprs) = as.matrix(assay(sce_scMerge, exprs))
# cond <- sce_scMerge$treatmentSt %in% c("UTTT", "UUUU") & sce_scMerge$clone %in% c('R', 'H')
plot_umap(sce_scMerge, exprs, metadata, output_dir, datatag, dims=1:15)
# plot_tsne(sce_scMerge, exprs, metadata, output_dir, npcs=20, tsne_perplex=90)


# Get SCTransform output for 3 series, using filtered data with less zeros
# Genes clustering, group genes in the same cluster, get in-cis, in-trans infos.  
# Line plots, keep values 
# Heatmap plots 
#

