source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA609'
# metacells_fn <-  paste0(input_dir,'SA609_rna/snakemake_10x/SA609_10x.csv')
# sce <- load_sce(metacells_fn, datatag, input_dir, return_data=T)

output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')

save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
dir.create(save_dir)


sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/','total_sce_v3.rds')
sce_raw <- readRDS(sce_fn)
dim(sce_raw)




  
# rowData(sce_raw)$ID[1]
# rowData(sce_raw)$Symbol[1]


res <- plot_raw_variance(sce_raw, save_dir)
# plot_rawdata_qc(sce_raw, save_dir, datatag)

output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
dir.create(output_dir)


print("Twice scran normalization")
source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/twice_scran_normalization_corrected.R"))
# sce_scran <- twice_scran_normalize_v2(sce_raw, input_dir, output_dir, datatag, return_data=F)
sce_scran <- readRDS(paste0(output_dir, datatag, '_twice_scran_normalized_v3.rds'))
dim(sce_scran)

print("Seurat normalization")
# sce_seurat <- normalize_Seurat(sce_raw, input_dir, output_dir, datatag, return_data=F)
sce_seurat <- readRDS(paste0(output_dir, datatag, '_seurat_normalized_v3.rds'))
dim(sce_seurat)

print("SCTransform normalization")
# normalize_SCTransform(sce, output_dir, datatag, return_data=F)
sce_transform <- readRDS(paste0(output_dir, datatag,'_sctransform_normalized.rds'))
dim(sce_transform)
# unique(sce_transform$timepoint)
# sce <- sce_transform[,sce_transform$timepoint %in% c('X3','X7')]
# unique(sce$sample)
# dim(logcounts(sce))
# saveRDS(sce, file=paste0(output_dir, datatag,'_sctransform_normalized_X3_X7_Sierra.rds'))
# cond <- sce_transform$treatmentSt %in% c('UTTTT', 'UUUUU') & sce_transform$clone %in% c('R', 'H')
# sce <- sce_transform[,cond]
# de_desc <- 'SA609: UTTTT-A vs UUUUU-H'
# # eval_edgeR(sce, vol_609, de_genes, de_desc, datatag, save_dir)
# base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# de_genes <- read.csv(paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTTT_R_UUUUU_H/signif_genes.csv'), check.names = F, stringsAsFactors = F)
# # p609 <- eval_edgeR(sce, de_genes, de_desc, datatag, save_dir)
# 
# treatmentSts <- c('UTTTT', 'UUUUU')
# clones <- c('R', 'H')
# desc <- '(X7-Rx A - X7-UnRx H)'
# # c_609 <- eval_edgeR_v2(sce, de_genes, treatmentSts, clones, de_desc, desc, datatag, save_dir)
# rc_609 <- eval_edgeR_counts(sce_raw, de_genes, treatmentSts, clones, paste0('raw data ',de_desc), 
#                             desc, datatag, save_dir)
# sce_scMerge <- readRDS(paste0(output_dir, datatag,'_corrected_scMerge.rds'))
# dim(sce_scMerge)


print("SCTransform normalization + Ridge Regression")
corrected_fn <- paste0(save_dir, datatag, "_norm_batch_corrected_sce.rds")
corrected_sce <- readRDS(corrected_fn)


metacells_fn <-  paste0(input_dir,'SA609_rna/snakemake_10x/SA609_10x.csv')
sample_df <- read.csv(metacells_fn, stringsAsFactors = F, check.names = F)
print(dim(sample_df))
sample_df <- sample_df %>%
  dplyr::select(mouse_id, batch_info) %>%
  dplyr::rename(batch=batch_info)

rownames(sample_df) <- sample_df$mouse_id
sce_raw$batch <- 'None'
sce_raw$batch <- sample_df[sce_raw$sample,'batch']
# sce_raw$series <- paste0(sce_raw$series,'_',gsub("CHIP0","C",sce_raw$batch))
summary(as.factor(sce_raw$batch))
sce_scran$batch <- sample_df[sce_scran$sample,'batch']
# sce_scran$series <- paste0(sce_scran$series,'_',gsub("CHIP0","C",sce_scran$batch))

sce_seurat$batch <- sample_df[sce_seurat$sample,'batch']
# sce_seurat$series <- paste0(sce_seurat$series,'_',gsub("CHIP0","C",sce_seurat$batch))

sce_transform$batch <- sample_df[sce_transform$sample,'batch']
# sce_transform$series <- paste0(sce_transform$series,'_',gsub("CHIP0","C",sce_transform$batch))


unique(sce_transform$series)
hk_ref <- res$filtered_gene_attr_HK
hk_ref <- hk_ref %>%
  dplyr::filter(abs(log_var)<0.2)
stable_genes <- hk_ref$gene_id
length(stable_genes)
# length(intersect(stable_genes,rownames(sce_scMerge)))
summary(res$filtered_gene_attr_HK$log_var)
metacells <- plot_batch_estimation(sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, 
                              datatag, 'HK', T, c(0.2,1.7))

dim(metacells)
p_hk <- plot_batch_effects(metacells, xstring="series", ystring="mean_exp", 
                               plottype="batch_label", plottitle=paste0(datatag,"\nNormalization Evaluation"),
                               xlabel=NULL, ylabel=paste0(length(stable_genes),' housekeeping genes average expression'),
                               lg_pos="bottom", save_dir)


# Get batch effect output here
meta_ridgeBER <- plot_stable_genes_exp_v2(corrected_sce, stable_genes, use_raw=F, exprs='logcounts', 
                                  plottitle=paste0(datatag,': Batch Effect Correction '),
                                  subtitle = paste0(length(stable_genes),' ',genes_type,' genes avg exp'),
                                  xlabel='', ylabel="Avg. genes exp", yl=NULL, legend_visible=T)
bls <- unique(meta_ridgeBER$batch)
length(bls)x
dim(meta_ridgeBER)
sum(bls %in% unique(meta_SA609$batch))
length(unique(meta_SA609$cell_id))
meta_ridgeBER <- meta_ridgeBER %>%
  dplyr::filter(cell_id %in% unique(meta_SA609$cell_id))

col_use <- c("tag","sample","batch","cell_id","clone","library_id","treatmentSt","mean_exp","timepoint","series")
meta_ridgeBER$tag <- 'Batch Corrected'
meta_ridgeBER <- meta_ridgeBER %>%
  dplyr::select(all_of(col_use))
dim(meta_ridgeBER)
unique(meta_SA609_v2$tag)
meta_SA609_v2 <- dplyr::bind_rows(meta_SA609, meta_ridgeBER) 
saveRDS(meta_SA609_v2, paste0(save_dir,"batch_effect_HK_plots_ridgeReg.rds"))
saveRDS(meta_SA609_v2, paste0(save_dir,datatag,"_HK_eval_metadata_ridgeReg.rds"))

meta_SA609_v2 <- meta_SA609_v2 %>%
  dplyr::filter(tag != 'Seurat')
meta_SA609_v2 <- get_metadata(meta_SA609_v2)

p609_hk <- plot_batch_effects(meta_SA609_v2, xstring="series", ystring="mean_exp", 
                              plottype="batch_label", plottitle=paste0(datatag),
                              xlabel=NULL, ylabel=paste0('512 HK genes avg exp'),
                              lg_pos="bottom", save_dir)
saveRDS(p609_hk, paste0(save_dir,"batch_effect_HK_plots_BE.rds"))

# p_hk <- plot_batch_estimation_v2(sce_scMerge, sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'HK', T)


## Plot scSEG genes list 
# ref <- res$filtered_gene_attr_scSEG
# ref <- ref %>%
#   dplyr::filter(abs(log_var)<0.2)
# stable_genes <- ref$gene_id
# length(stable_genes)
# p_ref_scSEG <- plot_batch_estimation_v2(sce_scMerge, sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'ref scSEG', T)
# 
# ref <- res$filtered_our_gene_attr_scSEG
# ref <- ref %>%
#   dplyr::filter(abs(log_var)<0.35)
# stable_genes <- ref$gene_id
# length(stable_genes)
# 
# p_custom_scSEG <- plot_batch_estimation_v2(sce_scMerge, sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'custom scSEG', T)
# 
# p_total <- cowplot::plot_grid(p_hk, p_ref_scSEG, p_custom_scSEG, ncol = 3, align='hv', rel_widths = c(1.45,1.2,1.2))
# png(paste0(save_dir,datatag,"_total_eval.png"), height = 2*970, width=2*1200,res = 2*72)
# print(p_total)
# dev.off()
# 



