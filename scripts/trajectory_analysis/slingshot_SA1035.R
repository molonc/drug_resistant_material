# script_dir <- '/home/htran/Projects/farhia_project/rnaseq/trajectory_analysis/'
script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/trajectory_analysis/'
source(paste0(script_dir, "slingshot_utils.R"))
source(paste0(script_dir, "tradeseq_utils.R"))


# Slingshot steps: 
# Read normalized data
# Extract PCA, UMAP, prepare_data_Seurat
# get lineage
# embedding to a feature reduction map
# plot lineage
datatag <- 'SA1035'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
if(!file.exists(save_dir)){
  dir.create(save_dir)
} 
output_dir <- save_dir
nfeatures_use <- 3000


# Prepare input data
# sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA1035-v6/SA1035_total_sce.rds'))
# dim(sce)
# print(rowData(sce)$Symbol[1])
# mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
# print(sum(mito_genes==TRUE))
# ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
# print(sum(ribo_genes==TRUE))
# print(dim(sce))
# assayNames(sce)
# norm_sce <- normalize_SCTransform(sce, output_dir, datatag, return_data=F, output_fn=NULL)
# norm_sce <- sce2
# norm_sce <- prepare_data_Seurat(sce, output_dir, datatag, save_srt=FALSE)

# sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
# saveRDS(sce, paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce_modified_clone.rds'))
# sce <- readRDS(paste0(output_dir, 'tradeseq/',datatag,'_',nfeatures_use,'_rd_sce_modified_clone.rds'))
# saveRDS(sce, paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce_v2.rds')) # the most updated
sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce_v2.rds'))

# summary(as.factor(sce$clone))
# summary(as.factor(sce$cluster_label))
# print(dim(sce))
# table(sce$cluster_label, sce$clone)
start_cls <- '0'
# get_slingshot_pseudotime_v2(sce, save_dir, datatag, start_cls=NULL, cl=NULL, rd_use='PCA')
crv_umap_embed <- readRDS(paste0(save_dir, "slingshot_output/slingshot_SA1035_0_UMAP_embed_crv.rds"))
output_figs_dir <- paste0(output_dir,'figs_v3/')
if(!dir.exists(output_figs_dir)){
  dir.create(output_figs_dir) 
}
## For pseudotime clonal labels, revision manuscript
# data.table::fwrite(umap_df, '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/clone_labels_unique_SA1035.csv.gz')

plot_all_lingeages(sce, crv_umap_embed, output_figs_dir, datatag)


segment_CNV_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/SA1035_cnv_mat.csv'
sign_genes_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/SA1035_total_genes_modules_act_repr_trans_08_Dec.csv.gz'
get_cis_trans_genes(sce, crv_umap_embed, output_dir, datatag, 
                    segment_CNV_fn, sign_genes_fn, start_cls)



