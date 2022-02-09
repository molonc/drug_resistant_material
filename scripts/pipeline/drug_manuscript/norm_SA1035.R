# Rscript compute_scSEG_stable_genes.R 
# -i /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA1035-v6/total_sce_v3.rds 
# -o /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/SA1035/segIndx_total.csv

source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA1035'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
# dir.create(save_dir)
output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
dir.create(output_dir)

sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/',datatag, '_total_sce_v3.rds')
sce_raw <- readRDS(sce_fn)
dim(sce_raw)
rowData(sce_raw)$ID[1]
rowData(sce_raw)$Symbol[1]
res <- plot_raw_variance(sce_raw, save_dir, datatag)
# options("scipen"=-100, "digits"=4)
# p1035_raw <- plot_rawdata_qc(sce_raw, save_dir, datatag)
# p1035_raw



print("Twice scran normalization")
# source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/twice_scran_normalization_corrected.R"))
# sce_scran <- twice_scran_normalize_v2(sce_raw, input_dir, output_dir, datatag, return_data=T)
sce_scran <- readRDS(paste0(output_dir, 'SA1035__twice_scran_normalized_v3.rds'))  #twice_scran_normalized_v3.rds
dim(sce_scran)
print("Seurat normalization")
# sce_seurat <- normalize_Seurat(sce_raw, input_dir, output_dir, return_data=T)
sce_seurat <- readRDS(paste0(output_dir, 'SA1035__seurat_normalized_v3.rds')) #seurat_normalized_v3.rds
dim(sce_seurat)

print("SCTransform normalization")
# normalize_SCTransform(sce_raw, output_dir, datatag, return_data=F)
sce_transform <- readRDS(paste0(output_dir, datatag,'_sctransform_normalized.rds'))
dim(sce_transform)


# treatmentSts <- c('UTTTT', 'UUUUU')
# clones <- c('H', 'E')
# cond <- sce_transform$treatmentSt %in% treatmentSts & sce_transform$clone %in% clones
# sce <- sce_transform[,cond]
# dim(sce)
# de_desc <- 'SA1035: UTTTT-H vs UUUUU-E'
# desc <- '(X8-Rx H - X8-UnRx E)'
# base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# de_genes <- read.csv(paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/SA1035_UTTTT_H_UUUUU_E/signif_genes.csv'), check.names = F, stringsAsFactors = F)
# dim(de_genes)
# # c_1035 <- eval_edgeR_v2(sce, de_genes, treatmentSts, clones, de_desc, desc, datatag, save_dir)
# rc_1035 <- eval_edgeR_counts(sce_raw, de_genes, treatmentSts, clones, paste0('raw data ',de_desc), desc, datatag, save_dir)
## eval_edgeR(sce, vol_1035, de_genes, de_desc, datatag, save_dir)


# 
# p1035 <- eval_edgeR(sce, de_genes, de_desc, datatag, save_dir)
# eval_edgeR_v2(sce, markers_ls, treatmentSts, clones, de_desc, datatag, output_dir)



# sce_scMerge <- readRDS(paste0(output_dir, datatag,'_corrected_scMerge.rds'))
# dim(sce_scMerge)


print("SCTransform normalization + Ridge Regression")
corrected_fn <- paste0(save_dir, datatag, "_norm_batch_corrected_sce.rds")
sce_corrected <- readRDS(corrected_fn)




metacells_fn <-  paste0(input_dir,'SA1035_rna/snakemake_10x_SA1035/metasample_SA1035.csv')
# metacells_fn <-  paste0(input_dir,'SA609_rna/snakemake_10x/SA609_10x.csv')
sample_df <- read.csv(metacells_fn, stringsAsFactors = F, check.names = F)

print(dim(sample_df))
sample_df <- sample_df %>%
  dplyr::select(mouse_id, batch_info) %>%
  dplyr::rename(batch=batch_info)

rownames(sample_df) <- sample_df$mouse_id
sce_raw$batch[1]
sce_transform$batch[1]
sce_scran$batch[1]
sce_seurat$batch[1]
# sce_raw$batch <- 'None'
# sce_raw$batch <- as.character(sample_df[sce_raw$sample,]$batch)
sce_raw$series <- paste0(sce_raw$series,'_',gsub("CHIP0","C",sce_raw$batch))
summary(as.factor(sce_raw$series))

# sce_transform$batch <- as.character(sample_df[sce_transform$sample,]$batch)
sce_transform$series <- paste0(sce_transform$series,'_',gsub("CHIP0","C",sce_transform$batch))

# sce_scran$batch <- as.character(sample_df[sce_scran$sample,]$batch)
sce_scran$series <- paste0(sce_scran$series,'_',gsub("CHIP0","C",sce_scran$batch))

# sce_seurat$batch <- as.character(sample_df[sce_seurat$sample,]$batch)
sce_seurat$series <- paste0(sce_seurat$series,'_',gsub("CHIP0","C",sce_seurat$batch))

hk_ref <- res$filtered_gene_attr_HK
hk_ref <- hk_ref %>%
  dplyr::filter(abs(log_var)<0.2)
stable_genes <- hk_ref$gene_id
length(stable_genes)
# length(intersect(stable_genes,rownames(sce_scMerge)))
summary(res$filtered_gene_attr_HK$log_var)




metacells_SA1035 <- plot_batch_estimation(sce_raw, sce_scran, sce_seurat, sce_transform, 
                                  stable_genes, save_dir, datatag, 'HK', T, c(0,1.4))
dim(metacells_SA1035)

p1035_hk <- plot_batch_effects(metacells_SA1035, xstring="series", ystring="mean_exp", 
                   plottype="batch_label", plottitle=paste0(datatag,"\nNormalization Evaluation"),
                   xlabel=NULL, ylabel=paste0(length(stable_genes),' housekeeping genes average expression'),
                   lg_pos="bottom", c(0,1.4), save_dir)

summary(as.factor(p1035_hk$tag))


# Get batch effect output here
p_ridgeBER <- plot_stable_genes_exp_v2(sce_corrected, stable_genes, use_raw=F, exprs='logcounts', 
                                          plottitle=paste0(datatag,': Batch Effect Correction '),
                                          subtitle = paste0(length(stable_genes),' ',genes_type,' genes avg exp'),
                                          xlabel='', ylabel="Avg. genes exp", yl=NULL, legend_visible=T)
bls <- unique(p_ridgeBER$batch)
length(bls)
dim(p_ridgeBER)
sum(bls %in% unique(meta_SA1035$batch))
length(unique(meta_SA1035$cell_id))
p_ridgeBER <- p_ridgeBER %>%
  dplyr::filter(cell_id %in% unique(meta_SA1035$cell_id))

col_use <- c("tag","sample","batch","cell_id","clone","library_id","treatmentSt","mean_exp","timepoint","series")
p_ridgeBER$tag <- 'Batch Corrected'
p_ridgeBER <- p_ridgeBER %>%
  dplyr::select(all_of(col_use))
dim(p_ridgeBER)

meta_SA1035 <- dplyr::bind_rows(meta_SA1035, p_ridgeBER) 
unique(meta_SA1035$tag)
saveRDS(meta_SA1035, paste0(save_dir,datatag,"_HK_eval_metadata_ridgeReg.rds"))

meta_SA1035 <- meta_SA1035 %>%
  dplyr::filter(tag != 'Seurat')
meta_SA1035 <- get_metadata(meta_SA1035)

p1035_hk <- plot_batch_effects(meta_SA1035, xstring="series", ystring="mean_exp", 
                              plottype="batch_label", plottitle=paste0(datatag),
                              xlabel=NULL, ylabel=paste0(length(stable_genes),' HK genes avg exp'),
                              lg_pos="bottom", save_dir)
saveRDS(p1035_hk, paste0(save_dir,"batch_effect_HK_plots_BE.rds"))

