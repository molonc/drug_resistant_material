
# Rscript compute_scSEG_stable_genes.R 
# -i /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA535-v6/SA535_cisplatin_total_sce_v3.rds 
# -o /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/SA535_cisplatin/segIndx_total.csv
# source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/preprocess_data.R"))

# metacells_fn <-  paste0(input_dir,'SA535_total_rna_v2/snakemake/metasample_SA535.csv')
# sce <- load_sce(metacells_fn, datatag, input_dir, return_data=T, tag)

source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/others_dataset/SA535_cisplatin/'
datatag <- 'SA535'
tag <- 'cisplatin'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/SA535_cisplatin/')
# dir.create(save_dir)
# if(tag!=''){
#   tag <- paste0(tag,'_')
# }
output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
# sce_fn <- paste0(paste0(input_dir,'rnaseq_v6/',datatag,'-v6/',datatag, '_',tag, '_total_sce_v3.rds'))
sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/',datatag, '_CIS_CX_untreated_total_sce_v3.rds') # _total_total_sce_v3.rds
sce_raw <- readRDS(sce_fn)
dim(sce_raw)

metacells_fn <-  paste0(input_dir,'SA535_total_rna_v2/snakemake/metasample_SA535.csv')
# metacells_fn <-  paste0(input_dir,'SA609_rna/snakemake_10x/SA609_10x.csv')
sample_df <- read.csv(metacells_fn, stringsAsFactors = F, check.names = F)
# sample_df$batch_info <- get_batch_infos(as.character(sample_df$library_id))
# write.csv(sample_df, file=metacells_fn, quote = F, row.names = F)
# head(sample_df)
if(tag=='cisplatin' & datatag=='SA535'){
  print(tag)
  sample_df <- sample_df %>%
    dplyr::filter(PDX %in% c('SA535_untreated','SA535_CY'))
}
print(dim(sample_df))
sample_df <- sample_df %>%
  dplyr::select(library_id, mouse_id, batch_info) %>%
  dplyr::rename(batch=batch_info)

# sce_raw$library_id[1]
summary(as.factor(sce_raw$library_id))
summary(as.factor(sce_raw$batch))
length(sample_df$library_id)
dim(sce_raw)
# sce_raw_origin <- sce_raw
sce_raw <- sce_raw[ ,sce_raw$library_id %in% sample_df$library_id]
dim(sce_raw)
res <- plot_raw_variance(sce_raw, save_dir, paste0(datatag,' ',tag))
# p535_cis_raw <- plot_rawdata_qc(sce_raw, save_dir, paste0(datatag,' ',tag))


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

# ts <- unique(sce_transform$treatmentSt)
# ts <- ts[grepl('T$',ts)]
# cond <- sce_transform$treatmentSt %in% ts & sce_transform$clone %in% c('J','Q')
# sce <- sce_transform[,cond]
# dim(sce)
# 
# 
# ts <- ts[grepl('U$',ts) & !grepl('TU$',ts)]
# cond <- sce_transform$treatmentSt %in% ts & sce_transform$clone %in% c('T','S','R')
# sce <- sce_transform[,cond]
# dim(sce)
# table(sce$treatmentSt, sce$clone)
# 
# genes_df <- data.frame(gene_id=rownames(sce))
# write.csv(genes_df, paste0(save_dir, datatag, "_genes_use.csv"), row.names = F, quote = F)
# corrected_data <- data.table::fread(paste0(save_dir, 'batch_effect_test/',datatag, "_untreated_batch_corrected_mtx.csv"), 
#                                     check.names=F, stringsAsFactors=F, header=T)
# dim(corrected_data)
# 
# rownames(corrected_data) <- corrected_data$V1
# corrected_data$V1 <- NULL
# sum(colnames(corrected_data)==colnames(sce))
# sum(rownames(corrected_data)==rownames(sce))
# 
# dim(sce)
# colnames(sce)[1]
# colnames(corrected_data)[1:10]
# sce1 <- sce
# logcounts(sce1) <- as.matrix(corrected_data)
# max(logcounts(sce1))
# max(logcounts(sce))
# dim(sce1)
# srt <- Seurat::as.Seurat(sce, counts = "counts", data="logcounts")
# srt <- Seurat::as.Seurat(sce1, counts = "counts", data="logcounts")
# 
# srt[["RNA"]]@scale.data <- as.matrix(normcounts(sce))
# # srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
# srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 3000)
# top_genes <- VariableFeatures(srt)
# top_genes[1]
# length(top_genes)
# head(top_genes)
# sce1 <- sce[top_genes,]
# dim(sce1)
# genes_df <- data.frame(ens_gene_id=rownames(sce1), gene_symbol=rowData(sce1)$Symbol)
# head(genes_df)
# save_dir <- '/home/htran/storage/python_workspace/velocity/result_SA535/cisplatin/'
# write.csv(genes_df, paste0(save_dir, datatag, "_top_3000_filtered_genes.csv"), row.names = F, quote = F)
# 
# # all.genes <- rownames(pbmc)
# # pbmc <- ScaleData(pbmc, features = all.genes)
# # srt
# srt <- Seurat::ScaleData(srt, features = rownames(srt))
# srt <- RunPCA(srt, verbose = F)
# srt <- RunUMAP(srt, dims = 1:30, verbose = FALSE)
# srt <- FindNeighbors(srt, dims = 1:30, verbose = FALSE)
# srt <- FindClusters(srt, verbose = FALSE, resolution = 0.3)
# # p_umap <- DimPlot(srt, label = TRUE, group.by = 'treatmentSt') #+ NoLegend()
# # p_umap1 <- DimPlot(srt, label = TRUE, group.by = 'batch') #+ NoLegend()
# # p_umap2 <- DimPlot(srt, label = TRUE, group.by = 'seurat_clusters') #+ NoLegend()
# # p_umap3 <- DimPlot(srt, label = TRUE, group.by = 'clone')
# 
# c_umap <- DimPlot(srt, label = TRUE, group.by = 'treatmentSt') #+ NoLegend()
# c_umap1 <- DimPlot(srt, label = TRUE, group.by = 'batch') #+ NoLegend()
# 
# c_umap2 <- DimPlot(srt, label = TRUE, group.by = 'seurat_clusters') #+ NoLegend()
# c_umap3 <- DimPlot(srt, label = TRUE, group.by = 'clone')
# 
# pm <- p_umap1 + p_umap3 + p_umap2 + patchwork::plot_layout(ncol=3)
# cm <- c_umap1 + c_umap3 + c_umap2 + patchwork::plot_layout(ncol=3)
# png(paste0(save_dir, datatag, "_clusters_umap_sctransform_normalization_untreated.png"), height = 2*300, width=2*1000,res = 2*72)
# print(pm)
# dev.off()
# 
# png(paste0(save_dir, datatag, "_clusters_umap_corrected_normalization_untreated.png"), height = 2*300, width=2*1000,res = 2*72)
# print(cm)
# dev.off()
# 
# meta_info <- srt@meta.data
# class(meta_info)
# # meta_info$seurat_clusters
# class(meta_info)
# sum(meta_info$cell_id==pca_df$cell_id)
# dim(srt@reductions$umap)
# umap_df <- srt[["umap"]]
# pca_df <- as.data.frame(Embeddings(object = srt, reduction = "pca"))
# umap_df <- as.data.frame(Embeddings(object = srt, reduction = "umap"))
# class(umap_df)
# class(pca_df)
# head(pca_df)
# dim(umap_df)
# head(umap_df)
# p <- DimPlot(srt, reduction = "pca")
# p
# unique(meta_info$seurat_clusters)
# meta_info <- meta_info %>%
#   dplyr::select(batch, seurat_clusters, cell_id, Barcode, 
#                 library_id, sample, clone, treatmentSt, timepoint)
# 
# write.csv(meta_info, paste0(save_dir, datatag, "_meta_info_untreated.csv"), row.names = F, quote = F)
# 
# norm_data <- as.data.frame(logcounts(sce))
# data.table::fwrite(t(norm_data), file=paste0(save_dir, datatag, "_norm_untreated.csv.gz"), row.names=F, quote=F)
# dim(norm_data)
# 
# umap_df$cell_id <- rownames(umap_df)
# umap_df <- umap_df %>% inner_join(meta_info, by=c("cell_id"))
# umap_df$cell_id <- paste0(umap_df$library_id,'_',umap_df$Barcode)
# 
# pca_df$cell_id <- rownames(pca_df)
# pca_df <- pca_df %>% inner_join(meta_info, by=c("cell_id"))
# pca_df$cell_id <- paste0(pca_df$library_id,'_',pca_df$Barcode)
# dim(umap_df)
# dim(pca_df)
# # test_sce <- readRDS(paste0(input_dir,'rnaseq_v6/SA535-v6/SA535X10XB03696.rdata'))
# # t <- colnames(colData(sce))
# # write.table(t, file = paste0(save_dir,'colnames_sce.txt'), sep='\n')
# 
# write.csv(pca_df, paste0(save_dir, datatag, "_norm_pca.csv"), row.names = F, quote = F)
# write.csv(umap_df, paste0(save_dir, datatag, "_norm_umap.csv"), row.names = F, quote = F)

# cond <- sce_transform$treatmentSt %in% c('UUTTTT', 'UUUUUU') & sce_transform$clone %in% c('T', 'J')
# sce <- sce_transform[,cond]
# dim(sce)
# de_desc <- 'SA535 cisplatin: UUTTTT-T vs UUUUUU-J'
# # plttitle <-'SA535 cisplatin: UUTTTT-T vs UUUUUU-J'
# treatmentSts <- c('UUTTTT', 'UUUUUU')
# clones <- c('T', 'J')
# desc <- '(X10-Rx T - X9-UnRx J)'
# base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# de_genes <- read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UUTTTT_T_UUUUUU_J/signif_genes.csv'), check.names = F, stringsAsFactors = F)
# dim(de_genes)
# # c_535_cis <- eval_edgeR_v2(sce, de_genes, treatmentSts, clones, de_desc, desc, paste0(datatag,' ',tag), save_dir)
# rc_535_cis <- eval_edgeR_counts(sce_raw, de_genes, treatmentSts, clones, paste0('raw data ',de_desc), 
#                                 desc, paste0(datatag,' ',tag), save_dir)

# rm(de_genes)

# 
# p535_cis <- eval_edgeR(sce, de_genes, de_desc, datatag, save_dir)
# eval_edgeR(sce, vol_535_cis, de_genes, de_desc, datatag, save_dir)



# sce_scMerge <- readRDS(paste0(output_dir, datatag,'_',tag,'corrected_scMerge.rds'))
# dim(sce_scMerge)
# sce_scMerge$treatmentSt <- get_treatment_status(sce_scMerge$series)


print("SCTransform normalization + Ridge Regression")
corrected_fn <- paste0(save_dir, datatag, "_norm_batch_corrected_sce.rds")
sce_corrected <- readRDS(corrected_fn)
dim(sce_corrected)
sce_corrected$library_id[1]
sce_corrected_cis <- sce_corrected[ ,sce_corrected$library_id %in% sample_df$library_id]
dim(sce_corrected_cis)







# s <- unique(sce_raw$library_id)
# sce_raw_1 <- sce_raw[,is.na(sce_raw$library_id)]
# dim(sce_raw_1)
# colnames(sce_raw_1)[1:3]
# unique(sce_raw_1$sample)

# rownames(sample_df) <- sample_df$library_id
# sce_raw$batch <- 'None'
# metacells <- as.data.frame(colData(sce))
# metacells <- metacells %>% left_join(sample_df, by=c("sample"="mouse_id"))
# colData(sce) <- as.matrix(metacells)
# sapply(sample_df$mouse_id, function(s) {
#   sce[,sce$sample==s]$batch <- sample_df[s,'batch']
# })
# sce_raw$library_id <- get_lib_info(sce_raw)
# sce_raw$batch <- sample_df[sce_raw$library_id,'batch']
# # sce_raw$batch <- ifelse(is.na(sce_raw$batch),'CHIP0000',sce_raw$batch)
# summary(as.factor(sce_raw$batch))
# 
# table(sce_raw$batch, sce_raw$library_id)
# 
# unique(sce_transform$library_id)
# sce_transform$library_id <- get_lib_info(sce_transform)
# sce_transform$batch <- sample_df[sce_transform$library_id,'batch']
# # sce_transform$batch <- ifelse(is.na(sce_transform$batch),'CHIP0000',sce_transform$batch)
# sce_scran$Sample[1]
# sce_scran$library_id <- get_lib_info(sce_scran)
# sce_scran$batch <- sample_df[sce_scran$library_id,'batch']
# # sce_scran$batch <- ifelse(is.na(sce_scran$batch),'CHIP0000',sce_scran$batch)
# 
# sce_seurat$Sample[1]
# sce_seurat$library_id <- get_lib_info(sce_seurat)
# sce_seurat$batch <- sample_df[sce_seurat$library_id,'batch']
# # sce_seurat$batch <- ifelse(is.na(sce_seurat$batch),'CHIP0000',sce_seurat$batch)


hk_ref <- res$filtered_gene_attr_HK
hk_ref <- hk_ref %>%
  dplyr::filter(abs(log_var)<0.2)
stable_genes <- hk_ref$gene_id
length(stable_genes)
# length(intersect(stable_genes,rownames(sce_scMerge)))

# p_hk <- plot_batch_estimation(sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'HK', T)
# sce_raw$series <- paste0(sce_raw$series,'_',gsub("CHIP0","C",sce_raw$batch))
# sce_transform$series <- paste0(sce_transform$series,'_',gsub("CHIP0","C",sce_transform$batch))
# sce_seurat$series <- paste0(sce_seurat$series,'_',gsub("CHIP0","C",sce_seurat$batch))
# sce_scran$series <- paste0(sce_scran$series,'_',gsub("CHIP0","C",sce_scran$batch))
metacells <- plot_batch_estimation(sce_raw, sce_scran, sce_seurat, sce_transform, 
                              stable_genes, save_dir, paste0(datatag,':Cisplatin'), 'HK', T, c(0,1.5))

dim(metacells)
p_hk <- plot_batch_effects(metacells, xstring="series", ystring="mean_exp", 
                           plottype="batch_label", plottitle=paste0(datatag,":Cisplatin\nNormalization Evaluation"),
                           xlabel=NULL, ylabel=paste0(length(stable_genes),' housekeeping genes average expression'),
                           lg_pos="bottom", save_dir)

summary(as.factor(p_hk$tag))


p_ridgeBER <- plot_stable_genes_exp_v2(sce_corrected_cis, stable_genes, use_raw=F, exprs='logcounts', 
                                       plottitle=paste0(datatag,': Batch Effect Correction '),
                                       subtitle = paste0(length(stable_genes),' ',genes_type,' genes avg exp'),
                                       xlabel='', ylabel="Avg. genes exp", yl=NULL, legend_visible=T)
bls <- unique(p_ridgeBER$batch)
length(bls)
dim(p_ridgeBER)
sum(bls %in% unique(meta_SA535_cis$batch))
length(unique(meta_SA535_cis$cell_id))
p_ridgeBER <- p_ridgeBER %>%
  dplyr::filter(cell_id %in% unique(meta_SA535_cis$cell_id))

col_use <- c("tag","sample","batch","cell_id","clone","library_id","treatmentSt","mean_exp","timepoint","series")
p_ridgeBER$tag <- 'Batch Corrected'
p_ridgeBER <- p_ridgeBER %>%
  dplyr::select(all_of(col_use))
dim(p_ridgeBER)

meta_SA535_cis <- dplyr::bind_rows(meta_SA535_cis, p_ridgeBER) 
unique(meta_SA535_cis$tag)
saveRDS(meta_SA535_cis, paste0(save_dir,datatag,"_HK_eval_metadata_ridgeReg.rds"))

meta_SA535_cis <- meta_SA535_cis %>%
  dplyr::filter(tag != 'Seurat')
meta_SA535_cis <- get_metadata(meta_SA535_cis)

p535_cis_hk <- plot_batch_effects(meta_SA535_cis, xstring="series", ystring="mean_exp", 
                               plottype="batch_label", paste0(datatag,":Cisplatin"),
                               xlabel=NULL, ylabel=paste0(length(stable_genes),' HK genes avg exp'),
                               lg_pos="bottom", save_dir)
saveRDS(p535_cis_hk, paste0(save_dir,"batch_effect_HK_plots_BE.rds"))




## Plot scSEG
# p_hk <- plot_batch_estimation_v2(sce_scMerge, sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'HK', T)
# 
# ref <- res$filtered_gene_attr_scSEG
# ref <- ref %>%
#   dplyr::filter(abs(log_var)<0.2)
# stable_genes <- ref$gene_id
# length(stable_genes)
# p_ref_scSEG <- plot_batch_estimation_v2(sce_scMerge, sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'ref scSEG', T)
# 
# ref <- res$filtered_our_gene_attr_scSEG
# ref <- ref %>%
#   dplyr::filter(abs(log_var)<0.2)
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
# 
# sce_scMerge = scater::runPCA(sce_scMerge,
#                              exprs_values = "scMerge_fast")
# p1 <- scater::plotPCA(
#   sce_scMerge, 
#   colour_by = "clone")
# p2 <- scater::plotPCA(
#   sce_scMerge, 
#   colour_by = "treatmentSt")
# p <- cowplot::plot_grid(p1, p2, ncol = 2, align='h')
# png(paste0(save_dir,datatag,"_corrected_PCA.png"), height = 2*450, width=2*1000,res = 2*72)
# print(p)
# dev.off()
# 
save_dir <- '/home/htran/storage/python_workspace/velocity/result_SA535/cisplatin/'
top_genes <- read.csv(paste0(save_dir,'figures_cis_dynamic/SA535_treated_topgenes.csv'), check.names=F, stringsAsFactors=F)
dim(top_genes)
sum(!is.na(top_genes$likelihood))
top_genes_ext <- top_genes %>%
  dplyr::filter(!is.na(likelihood))  # & likelihood>=0.3
dim(top_genes_ext)
mito_genes <- str_detect(top_genes_ext$gene_symbol, "^MT\\-")
print(sum(mito_genes==TRUE))
ribo_genes <- str_detect(top_genes_ecoxt$gene_symbol, "^RP(L|S)")  # or ^RP[L|S]?
print(sum(ribo_genes==TRUE))
top_genes_ext <- top_genes_ext %>%
  dplyr::filter(!str_detect(gene_symbol, "^RP(L|S)"))

lm_genes <- read.csv(paste0(save_dir,'SA535cisplatin_all_changing_genes2.csv'), check.names=F, stringsAsFactors=F)
dim(lm_genes)
head(lm_genes)
lm_genes <- lm_genes %>%
  dplyr::filter(abs(slopeX)>=0.01 & abs(slopeY)>=0.001)

lm_genes <- lm_genes %>%
  dplyr::filter(abs(slopeX)>=0.01 & abs(slopeY)>=0.01)

genes_signf <- intersect(lm_genes$gene_symbol, top_genes_ext$gene_symbol)
length(genes_signf)

# 272/4121 mono-increase/decrease genes
# 272/526 velocity driver genes

