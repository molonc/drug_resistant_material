
# remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
# remotes::install_github("jlmelville/uwot")
# devtools::install_github(repo = 'ChristophH/sctransform')
# BiocManager::install(c("RcppAnnoy"))

# Load sce data
# Combining into 1 object 
# Then remove MT, ribo genes first 
# Then test normalization methods 
# Test twice scran here 
# Test Seurat, with scale factor 10e6 
# Test SCTransform 

source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation_v1/')
output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/')

datatag <- 'SA609'
sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce.rds')

# sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce_treated.rds')
sce <- readRDS(sce_fn)
dim(sce)
# print("Twice scran normalization")
# source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/twice_scran_normalization_corrected.R"))
# twice_scran_normalize_v2(sce, input_dir, output_dir, datatag, return_data=F)
# 
# print("Seurat normalization")
# normalize_Seurat(sce, input_dir, output_dir, return_data=F)
# 
# 
# print("SCTransform normalization")
# normalize_SCTransform(sce, input_dir, output_dir, return_data=F)


# sce_scran <- readRDS(paste0(output_dir,'twice_scran_normalized.rds'))
# sce_seurat <- readRDS(paste0(output_dir,'seurat_normalized.rds'))
# sce_transform <- readRDS(paste0(output_dir,'sctransform_normalized.rds'))
sce_scMerge <- readRDS(paste0(output_dir,'SA609_scMerge_correction_without_cosine.rds'))  #SA609_scMerge_correction_without_cosine_v2.rds  SA609_scMerge_correction.rds dim(sce_transform)
dim(sce_scMerge)

sce$batch <- as.character(sce_scMerge[,colnames(sce)]$batch)
# sce_transform$batch <- as.character(sce_scMerge[,colnames(sce_transform)]$batch)
# sce_scran$batch <- as.character(sce_scMerge[,colnames(sce_scran)]$batch)
# sce_seurat$batch <- as.character(sce_scMerge[,colnames(sce_seurat)]$batch)
# sce$batch <- as.character(sce_scMerge[,colnames(sce)]$batch)

# dim(assay(sce_scMerge, "scMerge_fast"))
# t <- assay(sce_scMerge, "scMerge_fast")
# max(t)
# min(t)
# assayNames(sce_scMerge) 


# Load scSEG genes list
scSEG_df <- read.csv(paste0(output_dir,'segIndx_filtered_80.csv'),check.names = F, stringsAsFactors = F)
dim(scSEG_df)
colnames(scSEG_df)
View(head(scSEG_df))
summary(scSEG_df$segIdx)
scSEG_df <- scSEG_df %>%
            dplyr::filter(segIdx >= median(segIdx))
scSEG_df <- scSEG_df[order(scSEG_df$segIdx,decreasing = T),]  
ntop <- 300
if(ntop>nrow(scSEG_df)){
  ntop <- nrow(scSEG_df)
}
scSEG_df <- scSEG_df[1:ntop,]
dim(scSEG_df)
# View(head(scSEG_df))
stable_genes <- intersect(scSEG_df$gene_ens, rownames(sce))
# stable_genes <- intersect(stable_genes, rownames(sce_scran))
# stable_genes <- intersect(stable_genes, rownames(sce_seurat))
# stable_genes <- intersect(stable_genes, rownames(sce_transform))
stable_genes <- intersect(stable_genes, rownames(sce_scMerge))
length(stable_genes)

p_raw <- plot_stable_genes_exp(sce, stable_genes, use_raw=T, exprs='logcounts', plottitle=paste0('Raw data - ',length(stable_genes),' stable genes mean exp'),
                      xlabel='', ylabel="Mean exp of stably expressed genes", yl=NULL)

p_scran <- plot_stable_genes_exp(sce_scran, stable_genes, use_raw=F, exprs='logcounts', plottitle=paste0('Twice scran normalize - ',length(stable_genes),' stable genes mean exp'),
                                  xlabel='', ylabel="Mean exp of stably expressed genes", yl=NULL)

# p_seurat <- plot_stable_genes_exp(sce_seurat, stable_genes, use_raw=F, exprs='normcounts', plottitle=paste0('Seurat + Scale - ',length(stable_genes),' stable genes mean exp'),
#                                xlabel='', ylabel="Mean exp of stably expressed genes", yl=NULL)

# p_sctransform <- plot_stable_genes_exp(sce_transform, stable_genes, use_raw=F, exprs='normcounts', plottitle=paste0('SCTransform + Scale - ',length(stable_genes),' stable genes mean exp'),
#                                   xlabel='', ylabel="Mean exp of stably expressed genes", yl=NULL)

p_scmerge <- plot_stable_genes_exp(sce_scMerge, stable_genes, use_raw=F, exprs='scMerge_fast', plottitle=paste0('scMerge - batch corrected - ',length(stable_genes),' stable genes mean exp'),
                                       xlabel='', ylabel="Mean exp of stably expressed genes", yl=NULL)


p <- cowplot::plot_grid(p_raw, p_scmerge, ncol = 1, align='v')
# p <- cowplot::plot_grid(p_raw, p_sctransform, p_scmerge, ncol = 1, align='v')
png(paste0(output_dir,"scSEG_corrected_v2.png"), height = 2*600, width=2*510,res = 2*72)
print(p)
dev.off()

p_raw <- plot_stable_genes_exp_v2(sce, stable_genes, use_raw=T, exprs='logcounts', plottitle=paste0('Raw data - ',length(stable_genes),' stable genes mean exp'),
                               xlabel='', ylabel="Mean exp of stably expressed genes", yl=NULL)
p_scran <- plot_stable_genes_exp_v2(sce_scran, stable_genes, use_raw=F, exprs='logcounts', plottitle=paste0('Twice scran normalize - ',length(stable_genes),' stable genes mean exp'),
                                 xlabel='', ylabel="Mean exp of stably expressed genes", yl=NULL)

p_scmerge <- plot_stable_genes_exp_v2(sce_scMerge, stable_genes, use_raw=F, exprs='scMerge_fast', plottitle=paste0('scMerge - batch corrected - ',length(stable_genes),' stable genes mean exp'),
                                   xlabel='', ylabel="Mean exp of stably expressed genes", yl=NULL)
p <- cowplot::plot_grid(p_raw, p_scran, p_scmerge, ncol = 1, align='v')
# p <- cowplot::plot_grid(p_raw, p_sctransform, p_scmerge, ncol = 1, align='v')
png(paste0(output_dir,"scSEG_evaluation_batch_effect_wholedata_SA609_batch_infos.png"), height = 2*900, width=2*510,res = 2*72)
print(p)
dev.off()

# png(paste0(output_dir,"scSEG_evaluation.png"), height = 2*1200, width=2*510,res = 2*72)
# print(p)
# dev.off()
# 
# 
# Load housekeeping genes list
scHK_df <- readxl::read_excel(paste0(output_dir,'HKGenes.xlsx'))
dim(scHK_df)
head(scHK_df)
scHK_df$stability_index <- scHK_df$`Stability index`
scHK_df <- scHK_df %>%
      dplyr::filter(stability_index > quantile(x=stability_index, probs = 0.75))

dim(scHK_df)
summary(scHK_df$stability_index)
stable_genes <- intersect(scHK_df$GeneSymbol, rowData(sce)$Symbol)
length(stable_genes)
# 
meta_genes <- data.frame(ens_gene=rowData(sce)$ID, gene_symb=rowData(sce)$Symbol,
                         row.names = rowData(sce)$ID, stringsAsFactors = F)
meta_genes <- meta_genes[meta_genes$gene_symb %in% stable_genes,]
stable_genes <- meta_genes$ens_gene
stable_genes <- intersect(stable_genes, rownames(sce_scran))
# stable_genes <- intersect(stable_genes, rownames(sce_seurat))
stable_genes <- intersect(stable_genes, rownames(sce_transform))
stable_genes <- intersect(stable_genes, rownames(sce_scMerge))
length(stable_genes)
p_raw <- plot_stable_genes_exp(sce, stable_genes, use_raw=T, exprs='logcounts', plottitle=paste0('Raw data - ',length(stable_genes),' HK genes mean exp'),
                               xlabel='', ylabel="Mean exp of HK genes", yl=NULL)

p_scran <- plot_stable_genes_exp(sce_scran, stable_genes, use_raw=F, exprs='logcounts', plottitle=paste0('Twice scran - ',length(stable_genes),' HK genes mean exp'),
                                 xlabel='', ylabel="Mean exp of HK genes", yl=NULL)

p_seurat <- plot_stable_genes_exp(sce_seurat, stable_genes, use_raw=F, exprs='normcounts', plottitle=paste0('Seurat + Scale - ',length(stable_genes),' HK genes mean exp'),
                                  xlabel='', ylabel="Mean exp of HK genes", yl=NULL)

p_sctransform <- plot_stable_genes_exp(sce_transform, stable_genes, use_raw=F, exprs='normcounts', plottitle=paste0('SCTransform + Scale - ',length(stable_genes),' HK genes mean exp'),
                                       xlabel='', ylabel="Mean exp of HK genes", yl=NULL)

p_scmerge <- plot_stable_genes_exp(sce_scMerge, stable_genes, use_raw=F, exprs='scMerge_fast', plottitle=paste0('scMerge - batch corrected - ',length(stable_genes),' HK genes mean exp'),
                                   xlabel='', ylabel="Mean exp of HK genes", yl=NULL)


p <- cowplot::plot_grid(p_raw, p_scran, p_sctransform, p_scmerge, ncol = 1, align='v')
p <- cowplot::plot_grid(p_raw, p_scran, p_seurat, p_sctransform, ncol = 1, align='v')

png(paste0(output_dir,"HK_genes_batch_correction_eval_treated_cells_SA609.png"), height = 2*1200, width=2*510,res = 2*72)
print(p)
dev.off()


p_raw <- plot_stable_genes_exp_v2(sce, stable_genes, use_raw=T, exprs='logcounts', plottitle=paste0('Raw data - ',length(stable_genes),' HK genes mean exp'),
                                  xlabel='', ylabel="Mean exp of HK genes", yl=NULL)
p_scran <- plot_stable_genes_exp_v2(sce_scran, stable_genes, use_raw=F, exprs='logcounts', plottitle=paste0('Twice scran normalize - ',length(stable_genes),' HK genes mean exp'),
                                    xlabel='', ylabel="Mean exp of HK genes", yl=NULL)

p_scmerge <- plot_stable_genes_exp_v2(sce_scMerge, stable_genes, use_raw=F, exprs='scMerge_fast', plottitle=paste0('scMerge - batch corrected - ',length(stable_genes),' HK genes mean exp'),
                                      xlabel='', ylabel="Mean exp of HK genes", yl=NULL)
p <- cowplot::plot_grid(p_raw, p_scran, p_scmerge, ncol = 1, align='v')
# p <- cowplot::plot_grid(p_raw, p_sctransform, p_scmerge, ncol = 1, align='v')
png(paste0(output_dir,"HK_genes_evaluation_batch_effect_wholedata_SA609_batch_infos.png"), height = 2*900, width=2*510,res = 2*72)
print(p)
dev.off()

# summary(as.factor(sce_scran$clone))
# summary(as.factor(sce_seurat$clone))
# summary(as.factor(sce_transform$clone))



# View(head(results$gene_cluster))
# dim(results$gene_cluster)
# gene_cluster <- results$gene_cluster
# gene_cluster <- gene_cluster %>%
#   dplyr::filter(ens_gene %in% rownames(sce_raw))
# dim(gene_cluster)
# meta_cells_df$cell_id[1:3]
# meta_cells_df <- meta_cells_df %>%
#                     dplyr::filter(cell_id %in% colnames(sce_raw))
# p_old_scran <- plot_genes_exp_v3(gene_cluster, meta_genes, norm_data,
#                   meta_cells_df, obs_clones, 
#                   output_dir, datatag, clone_aware=TRUE)
# 
# observed_cells <- meta_cells_df$cell_id
# observed_genes <- gene_cluster$ens_gene
# norm_scran <- get_norm_data(sce_scran, observed_genes, observed_cells, exprs='logcounts')
# norm_seurat <- get_norm_data(sce_seurat, observed_genes, observed_cells, exprs='normcounts')
# norm_sctransform <- get_norm_data(sce_transform, observed_genes, observed_cells, exprs='normcounts')
# 
# write.csv(gene_cluster, paste0(output_dir, 'gene_cluster.csv'), quote = F, row.names = F)
# write.csv(meta_genes, paste0(output_dir, 'meta_genes.csv'), quote = F, row.names = F)
# write.csv(meta_cells_df, paste0(output_dir, 'meta_cells_df.csv'), quote = F, row.names = F)
# p_twice_scran <- plot_genes_exp_v3(gene_cluster, meta_genes, norm_scran,
#                                  meta_cells_df, obs_clones, 
#                                  output_dir, datatag, clone_aware=TRUE)
# p_seurat <- plot_genes_exp_v3(gene_cluster, meta_genes, norm_seurat,
#                                  meta_cells_df, obs_clones, 
#                                  output_dir, datatag, clone_aware=TRUE)
# p_sctransform <- plot_genes_exp_v3(gene_cluster, meta_genes, norm_sctransform,
#                                  meta_cells_df, obs_clones, 
#                                  output_dir, datatag, clone_aware=TRUE)
# 
# 
# p <- cowplot::plot_grid(p_old_scran,p_twice_scran, p_seurat, p_sctransform, 
#                         labels=c('Previous Scran','Twice Scran','Seurat_Scaled','SCTransform_Scaled'),
#                         ncol=1, align = 'v')
# 
# png(paste0(output_dir,"genes_clustering_mean_exp_evaluation.png"), height = 2*1200, width=2*600,res = 2*72)
# print(p)
# dev.off()
# 
# 
# 
meta_cells <- as.data.frame(colData(sce))
# meta_cells$total_features_by_counts
# meta_cells$total_counts

meta_cells <- meta_cells %>%
          dplyr::group_by(series) %>%
          dplyr::mutate(nb_cells = n())

meta_cells$series <- paste0(meta_cells$series,'_',meta_cells$nb_cells,' cells')

p <- plot_variation_function(meta_cells, xstring="series", ystring="total_features_by_counts",
                             plottype="series", "Raw data total_features_by_counts",
                             '', 'total_features_by_counts')

p_c <- plot_variation_function(meta_cells, xstring="series", ystring="total_counts",
                             plottype="series", "Raw data total_counts",
                             '', 'total_counts')
p <- cowplot::plot_grid(p, p_c,
                        ncol=2, align = 'h')

png(paste0(output_dir,"raw_data_features.png"), height = 2*500, width=2*1300,res = 2*72)
print(p)
dev.off()
# 
# get_norm_data <- function(observed_sce, observed_genes, observed_cells, exprs='logcounts'){
#   observed_sce <- observed_sce[observed_genes, observed_cells]
#   # norm_data <- logcounts(observed_sce)
#   norm_data <- assay(observed_sce, exprs)
#   norm_data <- as.data.frame(as.matrix(norm_data))
#   norm_data[norm_data==0] <- NA
#   norm_data$gene <- rownames(observed_sce)
#   print(dim(norm_data))
#   return(norm_data)
# }
# 
# rownames(sce_scran)[1:3]
# dim(rowData(sce_seurat))
# head(rowData(sce_seurat))



datatag <- 'SA609'
obs_clones <- c('R')
obs_clones_untreated <- c('H')
obs_treatment_st <- c('UTTTT')
obs_untreated_st <- c('UUUUU')
de_x3 = read.csv(paste0(input_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTT_R_UUUU_H/signif_genes.csv'))
de_x4 = read.csv(paste0(input_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTTT_R_UUUUU_H/signif_genes.csv'))
exprs <- 'scMerge_fast'
sce_scMerge$treatmentSt <- get_treatment_status(sce_scMerge$series)
print(summary(as.factor(sce_scMerge$treatmentSt)))
print(summary(as.factor(sce_scMerge$clone)))
class(assay(sce_scMerge, exprs))
assay(sce_scMerge, exprs) = as.matrix(assay(sce_scMerge, exprs))

# up_regulated_genes <- de_x3 %>%
#                     dplyr::filter(logFC>0) %>%
#                     dplyr::select(ensembl_gene_id)
# sce_scMerge_untreated <- sce_scMerge[,sce_scMerge$treatmentSt %in% obs_untreated_st]
# sce_scMerge_treated <- sce_scMerge[,sce_scMerge$treatmentSt %in% obs_treatment_st]
# sce_scMerge <- assay(sce_scMerge, exprs)
# meta_data <- as.data.frame(colData(sce_scMerge))
# meta_data$mean_exp_untreated <- DelayedArray::colMeans(assay(sce_scMerge_untreated, exprs))
# meta_data$mean_exp_treated <- DelayedArray::colMeans(assay(sce_scMerge_treated, exprs))

treated_cond <- sce_scMerge$treatmentSt %in% obs_treatment_st & sce_scMerge$clone %in% obs_clones
gexp_treated <- get_mean_value_by_treatment_condition(sce_scMerge, treated_cond, exprs)
gexp_treated <- gexp_treated %>%
                      dplyr::filter(gene %in% de_x3$ensembl_gene_id)
dim(gexp_treated)
treated_cond <- sce_scMerge$treatmentSt %in% obs_untreated_st & sce_scMerge$clone %in% obs_clones_untreated
gexp_untreated <- get_mean_value_by_treatment_condition(sce_scMerge, treated_cond, exprs)
gexp_untreated <- gexp_untreated %>%
  dplyr::filter(gene %in% de_x3$ensembl_gene_id)
dim(gexp_untreated)

gexp_total <- rbind(gexp_treated, gexp_untreated)
dim(gexp_total)
head(gexp_total)

gexp_total <- gexp_total %>% inner_join(de_x4, by=c("gene"="ensembl_gene_id"))
gexp_total$gene_type <- ifelse(gexp_total$logFC>0,'up-regulated','down-regulated')

gexp_total_up <- gexp_total %>%
  dplyr::filter(gene_type == 'up-regulated')
summary(gexp_total_up$logFC)
gexp_total_down <- gexp_total %>%
  dplyr::filter(gene_type == 'down-regulated')
summary(gexp_total_down$logFC)

p <- plot_variation_function(gexp_total_up, xstring="treatmentSt", ystring="mean_exps", plottype="treatmentSt", plottitle="Up-regulated in UTTTT vs UUUUU",
                        xlabel='Treatment status - batch corrected data', ylabel="Mean gene exp")

p2 <- plot_variation_function(gexp_total_down, xstring="treatmentSt", ystring="mean_exps", plottype="treatmentSt", plottitle="Down-regulated in UTTTT vs UUUUU",
                             xlabel='Treatment status - batch corrected data', ylabel="Mean gene exp")
p2




dim(sce)
length(intersect(ref_scgenes, stable_genes))
data("segList_ensemblGeneID", package = "scMerge") 
ref_scgenes <- segList_ensemblGeneID$human$human_scSEG
ref_scgenes <- intersect(ref_scgenes, rownames(sce))
length(ref_scgenes)
length(stable_genes)

var_genes <- read.csv(paste0(output_dir,'var_genes.csv'), stringsAsFactors = F, check.names = F)
obs_genes <- rownames(markers_ls[markers_ls$avg_diff<0,])
p_raw <- plot_stable_genes_exp(sce, var_genes$var_gene, use_raw=T, exprs='logcounts', plottitle=paste0('Raw data - ',length(var_genes$var_gene),' stable genes mean exp'),
                                  xlabel='', ylabel="Mean exp of top variable genes", yl=NULL)
p_raw_ref <- plot_stable_genes_exp(sce_scMerge, obs_genes, use_raw=F, exprs=exprs, plottitle=paste0('scMerge - ',length(obs_genes),' var genes mean exp'),
                                  xlabel='', ylabel="Mean exp of top variable genes", yl=NULL)
p <- cowplot::plot_grid(p_raw, p_raw_ref, ncol = 1, align='v')
# p <- cowplot::plot_grid(p_raw, p_sctransform, p_scmerge, ncol = 1, align='v')
png(paste0(output_dir,"scSEG_ours_vs_ref_v3.png"), height = 2*600, width=2*510,res = 2*72)
print(p)
dev.off()
dim(cancer_ref_genes_df)
length(intersect(cancer_ref_genes_df$ensemble_id, rownames(sce_scMerge)))
length(intersect(cis_20$gene, rowData(sce_scMerge)$Symbol))
length(intersect(cosmic_genes$Gene_Symbol, rowData(sce_scMerge)$Symbol))
length(intersect(ref_cis_genes$gene_symbol, rowData(sce_scMerge)$Symbol))
dim(ref_cis_genes)



# sce_scMerge
# exprs <- 'scMerge_fast'
# norm_data <- as.data.frame(assay(sce_scMerge, exprs))
# sum(norm_data<0)
# dim(norm_data)
# count_data <- as.data.frame(assay(sce_scMerge, "counts"))
# sum(count_data<0)
# summary(as.numeric(norm_data[2,]))
# norm_data1 <- log2(count_data+1)
# sum(norm_data1<0)
dim(srt[["RNA"]]@scale.data)
sum(srt[["RNA"]]@scale.data>0)

summary(as.data.frame(assay(sce_scMerge, exprs)))
median(assay(sce_scMerge, exprs))
min(assay(sce_scMerge, exprs))
max(assay(sce_scMerge, exprs))

sce_scMerge
cond <- sce_scMerge$treatmentSt %in% c(obs_treatment_st, obs_untreated_st) & sce_scMerge$clone %in% c(obs_clones, obs_clones_untreated)
cond <- sce_scMerge$treatmentSt %in% c("UTTT", "UUUU") & sce_scMerge$clone %in% c(obs_clones, obs_clones_untreated)
sum(cond==TRUE)
assay(sce_scMerge, exprs) = as.matrix(assay(sce_scMerge, exprs))
assay(sce_scMerge, "counts") = as.matrix(assay(sce_scMerge, "counts"))
dim(sce_scMerge[,cond])
srt <- Seurat::as.Seurat(sce_scMerge[,cond], counts = "counts", data = exprs) 
Idents(object = srt) <- 'clone'
# print(Idents(object = srt2))
# rownames(genes_map) <- genes_map$gene_ens
# Differential Expression Analysis 
groups_use <- c('R','H')
minLogFC <- 0.25
srt <- Seurat::ScaleData(object = srt, verbose = FALSE)
markers_ls <- FindMarkers(srt, slot="scale.data",  # using normalized data   
                          ident.1=groups_use[1], ident.2=groups_use[2],    #Vector of cell names belonging to group 1, group 2
                          logfc.threshold=minLogFC, test.use='wilcox')
print(paste0("Total DE marker genes: ",nrow(markers_ls)))
summary(as.numeric(markers_ls$avg_diff))
colnames(markers_ls)
markers_ls <- markers_ls %>%
              dplyr::filter(p_val_adj<0.05)
# write.csv(markers_ls, file=paste0(output_dir,'DE_signif_genes_R-UTTTT_vs_H-UUUUU.csv'), quote=F, row.names = F)
sum(markers_ls$avg_diff<0)
markers_ls$p_val_adj
pcm <- markers_ls %>% 
  ggplot(aes(avg_diff, -log10(p_val_adj))) +
  geom_point(size=2) 
pcm
obs_genes <- rownames(markers_ls[markers_ls$avg_log2FC<-0.3399,])
length(obs_genes)
p_down <- plot_stable_genes_exp(sce_scMerge[,cond], obs_genes, use_raw=F, exprs=exprs, plottitle=paste0('scMerge - ',length(obs_genes),' genes R-UTTT vs H-UUUU'),
                                   xlabel='', ylabel="Mean exp of down-regulated genes", yl=NULL)

obs_genes <- rownames(markers_ls[markers_ls$avg_log2FC>-0.3399,])
p_up <- plot_stable_genes_exp(sce_scMerge[,cond], obs_genes, use_raw=F, exprs=exprs, plottitle=paste0('scMerge - ',length(obs_genes),' genes R-UTTT vs H-UUUU'),
                                xlabel='', ylabel="Mean exp of up-regulated genes", yl=NULL)
p <- cowplot::plot_grid(p_up, p_down, nrow=1)
p

png(paste0(output_dir,"DE_analysis_corrected_R-UTTT_vs_H-UUUU.png"), height = 2*400, width=2*600,res = 2*72)
print(p)
dev.off()
write.csv(markers_ls, file=paste0(output_dir,'DE_signif_genes_R-UTTT_vs_H-UUUU.csv'), quote=F, row.names = F)

var_genes <- read.csv(paste0(output_dir,'var_genes.csv'), stringsAsFactors = F, check.names = F)
dim(var_genes)
obs_genes <- var_genes$var_gene
p <- plot_stable_genes_exp(sce_scMerge, obs_genes, use_raw=F, exprs=exprs, plottitle=paste0('scMerge - ',length(obs_genes),' top variable genes'),
                              xlabel='', ylabel="Mean exp of variable genes", yl=NULL, scale=F)

p1 <- plot_stable_genes_exp(sce_scMerge, obs_genes, use_raw=F, exprs=exprs, plottitle=paste0('scMerge - scaled genes ',length(obs_genes),' top variable genes'),
                           xlabel='', ylabel="Mean exp of scaled variable genes", yl=NULL, scale=T)

p2 <- cowplot::plot_grid(p, p1, nrow=1)
p2

png(paste0(output_dir,"top_200_variable_genes_batch_corrected_data.png"), height = 2*400, width=2*510,res = 2*72)
print(p)
dev.off()



