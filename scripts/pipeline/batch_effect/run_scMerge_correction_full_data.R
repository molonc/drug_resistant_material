source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))
# BiocManager::install("batchelor")
library(batchelor)
library(BiocSingular)
library(dplyr)
library(Matrix)
library(DelayedArray)

## SEG list in ensemblGene ID
data("segList_ensemblGeneID", package = "scMerge") 

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/')
dir.create(output_dir)
datatag <- 'SA609'
# sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce_treated.rds')
ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
cancer_ref_genes_fn <- paste0(ref_dif,'Behan_CFgenes.csv')
cancer_ref_genes_df <- read.csv(cancer_ref_genes_fn, stringsAsFactors=F, check.names = F)

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
cosmic_genes <- read.table(paste0(base_dir, 'biodatabase/oncogene_cosmic.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
ref_cis_genes <- read.table(paste0(base_dir, 'biodatabase/cisplatin_resistance_genes.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
dim(cosmic_genes)
dim(ref_cis_genes)
dim(cancer_ref_genes_df)





sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce.rds')
sce <- readRDS(sce_fn)
dim(sce)


meta_genes <- data.frame(ens_gene=rowData(sce)$ID, gene_symb=rowData(sce)$Symbol,
                         row.names = rowData(sce)$ID, stringsAsFactors = F)
cosmic_genes_ls <- meta_genes[meta_genes$gene_symb %in% cosmic_genes$Gene_Symbol,'ens_gene']
ref_cis_genes_ls <- meta_genes[meta_genes$gene_symb %in% ref_cis_genes$gene_symbol,'ens_gene']



scSEG_df <- read.csv(paste0(output_dir,'segIndx_filtered_80.csv'),check.names = F, stringsAsFactors = F)
dim(scSEG_df)
colnames(scSEG_df)
# View(head(scSEG_df))
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
print(length(stable_genes))
# Get number of highly variable features
sce_seurat <- readRDS(paste0(output_dir,'sctransform_normalized.rds')) #seurat_normalized.rds
zero_cbind1 <- DelayedArray::rowMeans(assay(sce_seurat, 'counts') == 0)
cut_off_overall <- 0.025
print(max(zero_cbind1))
print(min(zero_cbind1))
sce_seurat <- sce_seurat[names(zero_cbind1[zero_cbind1 <= (1 - cut_off_overall)]), ]
print(dim(sce_seurat))

srt <- Seurat::as.Seurat(sce_seurat, counts = "counts", data = "logcounts")
srt <- FindVariableFeatures(object = srt, verbose = FALSE, selection.method = "vst", nfeatures = 6500)
var_genes <- VariableFeatures(object = srt)
print(paste0('Nb var genes: ',length(var_genes)))

# var_genes <- union(segList_ensemblGeneID$human$human_scSEG, var_genes)
var_genes <- union(stable_genes, var_genes)
var_genes <- union(cosmic_genes_ls, var_genes)
var_genes <- union(ref_cis_genes_ls, var_genes)
var_genes <- union(cancer_ref_genes_df$ensemble_id, var_genes)
var_genes <- intersect(var_genes, rownames(sce))
print(length(var_genes))
rm(sce_seurat)
rm(srt)
sce <- sce[var_genes,]
# print("Cosine normalization as logcounts input")
# cosdata <- cosineNorm(as.matrix(log2(counts(sce)+1)), BPPARAM = SerialParam())
# colnames(cosdata) <- colnames(sce)
# rownames(cosdata) <- rownames(sce)
# logcounts(sce) <- as.matrix(cosdata)

print("Simple logcounts input")
logcounts(sce) <- as.matrix(log2(counts(sce)+1))
print(dim(sce))

metacells_fn <-  paste0(input_dir,'SA609_rna/snakemake_10x/SA609_10x.csv')
sample_df <- read.csv(metacells_fn)
# head(sample_df)
rownames(sample_df) <- sample_df$mouse_id
print(dim(sample_df))
sample_df <- sample_df %>%
  dplyr::select(mouse_id, batch_info) %>%
  dplyr::rename(batch=batch_info)
sce$batch <- 'None'
# metacells <- as.data.frame(colData(sce))
# metacells <- metacells %>% left_join(sample_df, by=c("sample"="mouse_id"))
# colData(sce) <- as.matrix(metacells)
# sapply(sample_df$mouse_id, function(s) {
#   sce[,sce$sample==s]$batch <- sample_df[s,'batch']
# })
sce$batch <- sample_df[sce$sample,'batch']
summary(as.factor(sce$batch))
kmeansK_params <- c()
for(b in unique(sce$batch)){
  # sce_tmp <- sce[,sce$batch==b]
  # length(unique(sce_tmp$clone))
  kmeansK_params <- c(kmeansK_params, length(unique(sce[,sce$batch==b]$clone)))
}

print(kmeansK_params)

assay(sce, "counts") = as(counts(sce), "dgeMatrix")
assay(sce, "logcounts") = as(logcounts(sce), "dgeMatrix")

t1 = Sys.time()
#segList_ensemblGeneID$human$human_scSEG,
scMerge_res <- scMerge::scMerge(
  sce_combine = sce, 
  ctl = stable_genes,  
  assay_name = "scMerge_fast",
  replicate_prop = 0.4,
  cell_type = NULL, # unsupervised
  kmeansK = kmeansK_params,
  verbose=T,
  BSPARAM = IrlbaParam(), 
  svd_k = 20)

t2 = Sys.time()
print(t2-t1)


saveRDS(scMerge_res, file=paste0(output_dir,datatag,"_scMerge_correction_without_cosine_v3.rds"))


t <- assay(scMerge_res, "scMerge_fast")
print(max(t))
print(min(t))




