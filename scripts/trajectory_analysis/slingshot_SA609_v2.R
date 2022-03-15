script_dir <- '/home/htran/Projects/farhia_project/rnaseq/trajectory_analysis/'
source(paste0(script_dir, "slingshot_utils.R"))
source(paste0(script_dir, "tradeseq_utils.R"))
# Slingshot steps: 
# Read normalized data
# Extract PCA, UMAP, prepare_data_Seurat
# get lineage
# embedding to a feature reduction map
# plot lineage

datatag <- 'SA609'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')

# save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
# norm_sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_normalized.rds'))

# # sce <- readRDS(paste0(save_dir,'SA609_3000_sce.rds'))
# sce <- readRDS(paste0(save_dir,'sce_raw_3000hvg.rds'))
# dir.create(save_dir)
# convert_corrected_mtx_to_sce(norm_sce, save_dir, 
#                              corrected_mtx_fn=NULL, datatag, 
#                              return_sce=F, save_data=T)
# sce <- readRDS(paste0(input_dir,'slingshot_trajectory/BE_mtx_v2/SA609_norm_BE_sce.rds'))
# sce <- readRDS(paste0(input_dir,'slingshot_trajectory/withBE_SA609_v2/SA609_3000_rd_sce.rds'))
print(dim(sce))
assayNames(sce)
df <- as.matrix(counts(sce))

library(SingleCellExperiment)
output_dir <- save_dir
prepare_data_Seurat(sce, output_dir, datatag, save_srt=FALSE)
sce <- readRDS(paste0(save_dir,'SA609_3000_rd_sce.rds'))
# sce$clone[1]
# summary(as.factor(sce$clone))
metacells <- data.table::fread(paste0(output_dir,'SA609_3000_meta_cells.csv')) %>% as.data.frame()
dim(metacells)

meta_info <- data.table::fread(paste0(output_dir,'BE_mtx_v2/SA609_meta_info.csv')) %>% as.data.frame()
dim(meta_info)

umap_df <- data.table::fread(paste0(output_dir, "withBE_SA609_v2/SA609_3000_norm_umap.csv")) %>% as.data.frame()
dim(umap_df)
meta_info <- meta_info %>%
  dplyr::filter(cell_id %in% umap_df$cell_id)

data.table::fwrite(meta_info, paste0(output_dir,'BE_mtx_v2/SA609_meta_info.csv'))

sce <- sce[,!sce$clone %in% c('None','unassigned','R1','B_R_R1')]
sce <- sce[,metacells$cell_id]
print(dim(sce))
excluded_samples <- c('SA609X5XB03235','SA609X4XB03084')

sce <- sce[,!sce$sample %in% excluded_samples]
print(dim(sce))
# readRDS(paste0(save_dir,'sce_raw_3000hvg.rds'))
# sce <- sce[,!sce$clone %in% c('R1','B_R_R1')]

# remove cls 12
# remove replicates
# change R to SR1, SR2, SR3

nfeatures_use <- 3000
meta_data <- data.table::fread(paste0(output_dir, datatag,'_',nfeatures_use, "_meta_cells.csv")) %>% as.data.frame()
summary(as.factor(meta_data$seurat_clusters))
meta_data <- meta_data %>%
  dplyr::filter(seurat_clusters!=12)
dim(meta_data)
meta_data$clone
# meta_data <- meta_data %>%
#   dplyr::filter(clone=='R')

meta_data <- meta_data %>%
  dplyr::filter(seurat_clusters %in% c(0, 1, 4, 10))
meta_data <- meta_data %>%
  dplyr::select(cell_id, seurat_clusters)
meta_info2 <- meta_info2 %>% inner_join(meta_data, by = c('cell_id'))
meta_info2$clone <- paste0('S',meta_info2$clone,meta_info2$seurat_clusters)
summary(as.factor(meta_info2$clone))
dim(sce)

sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))

sce1 <- sce1[,colnames(sce1) %in% meta_info$cell_id]
dim(sce1)
cells_use <- colnames(sce1)
sce <- sce[,colnames(sce) %in% cells_use]
dim(sce)
print(unique(sce$clone))
output_dir <- save_dir
save_srt=TRUE
prepare_data_Seurat(sce, output_dir, datatag, save_srt)



# Load data for trajectory plotting
output_dir <- save_dir
nfeatures_use <- 3000
sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
dim(sce)
# sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
# pca_df <- data.table::fread(paste0(output_dir, datatag,'_',nfeatures_use, "_norm_pca.csv")) %>% as.data.frame()
# umap_df <- data.table::fread(paste0(output_dir, datatag,'_',nfeatures_use, "_norm_umap.csv")) %>% as.data.frame()
# dim(umap_df)
start_cls <- '10'
rd_use <- 'PCA'
nfeatures_use <- 3000
crv_umap_embed <- readRDS(paste0(save_dir, "slingshot_",datatag,'_',paste(start_cls, collapse='_'),"_UMAP_embed_crv.rds"))
# crv1 <- readRDS(paste0(save_dir, "slingshot_pseudotime_SA609_10_PCA_crv.rds"))

res <- plot_all_lingeages(sce, crv_umap_embed, output_dir, datatag)

segment_CNV_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/SA609_cnv_mat.csv'
sign_genes_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz'
get_cis_trans_genes(sce, crv_umap_embed, output_dir, datatag, 
                    segment_CNV_fn, sign_genes_fn, start_cls)

# Plotting
# Look at the script Figure7_SA609.R



