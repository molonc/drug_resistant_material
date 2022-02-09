source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/others_dataset/SA535_cisplatin/'
datatag <- 'SA535'
tag <- 'CX5461'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'_',tag,'/')
dir.create(save_dir)
# if(tag!=''){
#   tag <- paste0(tag,'_')
# }
output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')


sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/',datatag, '_CIS_CX_untreated_total_sce_v3.rds') # _total_total_sce_v3.rds
sce_raw <- readRDS(sce_fn)

dim(sce_raw)

metacells_fn <-  paste0(input_dir,'SA535_total_rna_v2/snakemake/metasample_SA535.csv')
# metacells_fn <-  paste0(input_dir,'SA609_rna/snakemake_10x/SA609_10x.csv')
sample_df <- read.csv(metacells_fn, stringsAsFactors = F, check.names = F)
# sample_df$batch_info <- get_batch_infos(as.character(sample_df$library_id))
# write.csv(sample_df, file=metacells_fn, quote = F, row.names = F)
# head(sample_df)
if(tag=='CX5461' & datatag=='SA535'){
  print(tag)
  sample_df <- sample_df %>%
    dplyr::filter(PDX %in% c('SA535_untreated','SA535_CX'))
}
print(dim(sample_df))
# sample_df$mouse_id
# View(sample_df)
sample_df <- sample_df %>%
  dplyr::select(library_id, mouse_id, batch_info) %>%
  dplyr::rename(batch=batch_info)

sce_raw$library_id[1]
summary(as.factor(sce_raw$library_id))
summary(as.factor(sce_raw$batch))
length(sample_df$library_id)
dim(sce_raw)
sce_raw <- sce_raw[ ,sce_raw$library_id %in% sample_df$library_id]
dim(sce_raw)
res <- plot_raw_variance(sce_raw, save_dir, paste0(datatag,' ',tag))
# plot_rawdata_qc(sce_raw, save_dir, paste0(datatag,' ',tag))


# dir.create(output_dir)


# print("Twice scran normalization")
# source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/twice_scran_normalization_corrected.R"))
# sce_scran <- twice_scran_normalize_v2(sce_raw, input_dir, output_dir, paste0(datatag,'_',tag), return_data=T)
sce_scran <- readRDS(paste0(output_dir, datatag,'_total_twice_scran_normalized_v3.rds'))
dim(sce_scran)
sce_scran <- sce_scran[ ,sce_scran$library_id %in% sample_df$library_id]
dim(sce_scran)
print("Seurat normalization")
# sce_seurat <- normalize_Seurat(sce_raw, input_dir, output_dir, paste0(datatag,'_',tag), return_data=T)
sce_seurat <- readRDS(paste0(output_dir, datatag,'_total_seurat_normalized_v3.rds'))
dim(sce_seurat)
sce_seurat <- sce_seurat[ ,sce_seurat$library_id %in% sample_df$library_id]
dim(sce_seurat)
print("SCTransform normalization")
# sce_transform <- normalize_SCTransform(sce_raw, output_dir, paste0(datatag,'_',tag), return_data=F)
sce_transform <- readRDS(paste0(output_dir, datatag,'_total_sctransform_normalized.rds'))
dim(sce_transform)
sce_transform <- sce_transform[ ,sce_transform$library_id %in% sample_df$library_id]
dim(sce_transform)


print("SCTransform normalization + Ridge Regression")
corrected_fn <- paste0(save_dir, datatag, "_norm_batch_corrected_sce.rds")
sce_corrected <- readRDS(corrected_fn)
dim(sce_corrected)
sce_corrected$library_id[1]
dim(sample_df)
sce_corrected_cx <- sce_corrected[ ,sce_corrected$library_id %in% sample_df$library_id]
dim(sce_corrected_cx)

# ts <- unique(sce_transform$treatmentSt)
# ts <- ts[grepl('X$',ts)]
# cond <- sce_transform$treatmentSt %in% ts & sce_transform$clone %in% c('U','R')
# sce <- sce_transform[,cond]
# dim(sce)
# meta_cells <- as.data.frame(colData(sce))
# table(meta_cells$clone, meta_cells$treatmentSt)
# ts <- unique(sce$treatmentSt)
# genes_df <- data.frame(gene_id=rownames(sce))
# save_dir <- '/home/htran/storage/python_workspace/velocity/result_SA535/cisplatin/'
# 
# write.csv(genes_df, paste0(save_dir, datatag, "_genes_use.csv"), row.names = F, quote = F)
# 
# srt <- Seurat::as.Seurat(sce, counts = "counts", data="logcounts")
# srt[["RNA"]]@scale.data <- as.matrix(normcounts(sce))



# cond <- sce_transform$treatmentSt %in% c('UXXXX', 'UUUUU') & sce_transform$clone %in% c('U', 'J')
# sce <- sce_transform[,cond]
# dim(sce)
# de_desc <- 'SA535 CX5461: UXXXX-U vs UUUUU-J'
# # unique(sce$series)
# # sce$series <- ifelse(sce$series=='X8-UUUUU','2_X8-UUUUU','1_X8-UXXXX')
# 
# treatmentSts <- c('UXXXX', 'UUUUU')
# clones <- c('U', 'J')
# desc <- '(X8-Rx U - X9-UnRx J)'
# base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# de_genes <- read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UXXXX_U_UUUUU_J/signif_genes.csv'), check.names = F, stringsAsFactors = F)
# dim(de_genes)
# # c_535_cx <- eval_edgeR_v2(sce, de_genes, treatmentSts, clones, de_desc, desc, paste0(datatag,' ',tag), save_dir)
# 
# rc_535_cx <-eval_edgeR_counts(sce_raw, de_genes, treatmentSts, clones, paste0('raw data ',de_desc), desc, paste0(datatag,' ',tag), save_dir)


# # eval_edgeR(sce, vol_535_cx, de_genes, de_desc, datatag, save_dir)
# 

# p535_cx <- eval_edgeR(sce, de_genes, de_desc, datatag, save_dir)



hk_ref <- res$filtered_gene_attr_HK
hk_ref <- hk_ref %>%
  dplyr::filter(abs(log_var)<0.2)
stable_genes <- hk_ref$gene_id
length(stable_genes)
# length(intersect(stable_genes,rownames(sce_scMerge)))
# sce_raw$batch <- fix_batch(sce_raw)
sce_transform$batch <- fix_batch(sce_transform)
sce_scran$batch <- fix_batch(sce_scran)
sce_seurat$batch <- fix_batch(sce_seurat)

unique(sce_raw$batch)
unique(sce_transform$batch)
unique(sce_scran$batch)
unique(sce_seurat$batch)
# p_hk <- plot_batch_estimation(sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'HK', T)
# sce_raw$series <- paste0(sce_raw$series,'_',gsub("CHIP0","C",sce_raw$batch))
# sce_transform$series <- paste0(sce_transform$series,'_',gsub("CHIP0","C",sce_transform$batch))
# sce_seurat$series <- paste0(sce_seurat$series,'_',gsub("CHIP0","C",sce_seurat$batch))
# sce_scran$series <- paste0(sce_scran$series,'_',gsub("CHIP0","C",sce_scran$batch))


metacells <- plot_batch_estimation(sce_raw, sce_scran, sce_seurat, sce_transform, 
                              stable_genes, save_dir, paste0(datatag,':CX5461'), 'HK', T, c(0.1,1.4))

dim(metacells)
p_hk <- plot_batch_effects(metacells, xstring="series", ystring="mean_exp", 
                           plottype="batch_label", plottitle=paste0(datatag,":CX5461\nNormalization Evaluation"),
                           xlabel=NULL, ylabel=paste0(length(stable_genes),' housekeeping genes average expression'),
                           lg_pos="bottom", save_dir)

summary(as.factor(p_hk$tag))

# rm(p_ridgeBER)
p_ridgeBER <- plot_stable_genes_exp_v2(sce_corrected_cx, stable_genes, use_raw=F, exprs='logcounts', 
                                       plottitle=paste0(datatag,'_',tag,': Batch Effect Correction '),
                                       subtitle = paste0(length(stable_genes),' ',genes_type,' genes avg exp'),
                                       xlabel='', ylabel="Avg. genes exp", yl=NULL, legend_visible=T)
bls <- unique(p_ridgeBER$batch)
length(bls)
dim(p_ridgeBER)
sum(bls %in% unique(meta_SA535_cx$batch))
length(unique(meta_SA535_cx$cell_id))
p_ridgeBER <- p_ridgeBER %>%
  dplyr::filter(cell_id %in% unique(meta_SA535_cx$cell_id))

col_use <- c("tag","sample","batch","cell_id","clone","library_id","treatmentSt","mean_exp","timepoint","series")
p_ridgeBER$tag <- 'Batch Corrected'
p_ridgeBER <- p_ridgeBER %>%
  dplyr::select(all_of(col_use))
dim(p_ridgeBER)

meta_SA535_cx <- dplyr::bind_rows(meta_SA535_cx, p_ridgeBER) 
unique(meta_SA535_cx$tag)
batch_ls <- as.character(meta_SA535_cx$batch)
meta_SA535_cx$batch <- ifelse(meta_SA535_cx$batch=='CHIP001','CHIP0047',
                   ifelse(meta_SA535_cx$batch=='CHIP003' | meta_SA535_cx$batch=='CHIP004','CHIP0192',meta_SA535_cx$batch))
saveRDS(meta_SA535_cx, paste0(save_dir,datatag,"_HK_eval_metadata_ridgeReg.rds"))

meta_SA535_cx <- meta_SA535_cx %>%
  dplyr::filter(tag != 'Seurat')
meta_SA535_cx <- get_metadata(meta_SA535_cx)

p535_cx_hk <- plot_batch_effects(meta_SA535_cx, xstring="series", ystring="mean_exp", 
                                 plottype="batch_label", plottitle=paste0(datatag,":CX5461"),
                                 xlabel=NULL, ylabel=paste0(length(stable_genes),' HK genes avg exp'),
                                 lg_pos="bottom", save_dir)

saveRDS(p535_cx_hk, paste0(save_dir,"batch_effect_HK_plots_BE.rds"))


