# Run slingshot SA535


script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/'
source(paste0(script_dir, "trajectory_analysis/slingshot_utils.R"))

datatag <- 'SA535'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
# sce <- readRDS(paste0(save_dir,'SA535_sctransform_normalized.rds'))
# dim(sce)
# assayNames(sce)
output_dir <- save_dir
# "/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/"

# prepare_data_Seurat(sce, save_dir, datatag)



# Load pca vectors 
# Load umap vectors 
nfeatures_use <- 3000
# pca_df <- data.table::fread(paste0(output_dir, datatag,'_',nfeatures_use, "_norm_pca.csv")) %>% as.data.frame()
# umap_df <- data.table::fread(paste0(output_dir, datatag,'_',nfeatures_use, "_norm_umap.csv")) %>% as.data.frame()
# var_genes_df <- data.table::fread(paste0(output_dir,'SA535_3000_hvg_genes.csv')) %>% as.data.frame()
# meta_info <- data.table::fread(paste0(output_dir, datatag,'_',nfeatures_use, "_meta_cells.csv")) %>% as.data.frame()

var_genes <- var_genes_df$bc[1:300]
sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_sce.rds'))
# sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_sce_withoutX9_cls18.rds'))
# saveRDS(sce, paste0(output_dir, datatag,'_',nfeatures_use,'_sce_withoutX9_cls18.rds'))
unique(sce$timepoint)
sce <- sce[,sce$timepoint!="X9"]
sce <- sce[,sce$cluster_label!=18]
dim(sce)
sce
unique(sce$cluster_label)
meta_info <- as.data.frame(colData(sce))
prepare_data_Seurat(sce, output_dir, datatag, save_srt=TRUE)
  
  
  

start_cls <- '7'
# start_cls <- '12'
rd_use <- 'PCA'
nfeatures_use <- 3000
crv_umap_embed <- readRDS(paste0(save_dir, "slingshot_",datatag,start_cls,"_UMAP_embed_crv.rds"))
# crv1 <- readRDS(paste0(save_dir, "slingshot_pseudotime_",datatag,"_",start_cls,'_',rd_use,"_crv.rds"))
# sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce_clones.rds'))
dim(sce)
# colnames(sce)[1]
colnames(sce) <- sce$cell_id
sce$treatmentSt <- sce$treat
# rownames(sce)[1]
# sce$treat[1]
# saveRDS(sce, paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce_clones.rds'))
# saveRDS(sce, paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
plot_all_lingeages(sce, crv_umap_embed, output_dir, datatag)

segment_CNV_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/SA535_cnv_mat.csv'
sign_genes_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/SA535_total_genes_modules_act_repr_trans_08_Dec.csv.gz'
get_cis_trans_genes(sce, crv_umap_embed, output_dir, datatag, 
                    segment_CNV_fn, sign_genes_fn, start_cls)


