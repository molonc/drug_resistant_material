# Processing data for CPA 
library(stringr)
library(dplyr)

library(Seurat)
generate_input_data_CPA <- function(norm_sce, save_dir, datatag='SA'){
  # require("SingleCellExperiment")
  mito_genes <- str_detect(rowData(norm_sce)$Symbol, "^MT\\-")
  print(sum(mito_genes==TRUE))
  srt <- Seurat::as.Seurat(norm_sce, counts = "counts", data="logcounts", assay = "RNA", project = "SingleCellExperiment") #set to NULL if only normalized data are present
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 3000)
  
  length(VariableFeatures(srt))
  hvg <- VariableFeatures(srt)
  hvg[1]
  norm_data <- as.data.frame(logcounts(norm_sce[hvg,]))
  dim(norm_data)
  mtx_fn <- paste0(save_dir, datatag, "_norm.csv.gz")
  data.table::fwrite(norm_data, file=mtx_fn, row.names=T, quote=F)
  raw_data <- as.data.frame(counts(norm_sce[hvg,]))
  dim(raw_data)
  raw_mtx_fn <- paste0(save_dir, datatag, "_raw.csv.gz")
  data.table::fwrite(raw_data, file=raw_mtx_fn, row.names=T, quote=F)
  print(dim(norm_data))
  print(paste0('Save normalized mtx output as: ',mtx_fn))
  
}




get_metadata_CPA <- function(norm_sce, meta_fn){
  meta_info <- as.data.frame(colData(norm_sce))
  # t <- unique(meta_info$treatmentSt)
  meta_info <- meta_info %>% 
    dplyr::select(treatmentSt, clone, sample, cell_id, clone, library_id, batch)
  
  meta_info$drug <- ifelse(grepl('T',meta_info$treatmentSt),'CIS',
                           ifelse(grepl('X',meta_info$treatmentSt),'CX','Vehicle'))
  
  meta_info$dose <- ifelse(meta_info$drug=='CIS',str_count(meta_info$treatmentSt,'T'),
                           ifelse(meta_info$drug=='CX',str_count(meta_info$treatmentSt,'X'),1.0))
  meta_info$dose_val <- meta_info$dose/max(meta_info$dose)
  meta_info$dose_val <- ifelse(meta_info$drug=='Vehicle',1.0,meta_info$dose_val)
  print(table(meta_info$treatmentSt, meta_info$dose))
  print(summary(meta_info$dose_val))
  # unique(meta_info$dose_val)
  meta_info$cell_type <- 'TNBC'
  meta_info$drug_dose_name <- paste0(meta_info$drug,'_',meta_info$dose_val)
  unique(meta_info$drug_dose_name)
  meta_info$drug_dose_name <- ifelse(meta_info$drug_dose_name=='Vehicle_1', 'Vehicle_1.0', meta_info$drug_dose_name)
  meta_info$drug_dose_name <- ifelse(meta_info$drug_dose_name=='CIS_1', 'CIS_1.0', meta_info$drug_dose_name)
  
  
  meta_info$cov_drug_dose_name <- paste0(meta_info$cell_type,'_',meta_info$drug_dose_name)
  meta_info$condition <- meta_info$drug
  meta_info$control <- ifelse(meta_info$drug=='Vehicle', 1, 0)
  # meta_fn <- paste0(save_dir,'CPA/', datatag, "_meta_info_CPA.csv")
  write.csv(meta_info, meta_fn, row.names = F, quote = F)
  unique(meta_info$drug)
  return(meta_info)
  
}


library(dplyr)
datatag <- 'SA609'
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
norm_sce <- readRDS(paste0(output_dir, datatag,'_sctransform_normalized.rds'))
dim(norm_sce)
metacells_fn <-  paste0(input_dir,'SA609_rna/snakemake_10x/SA609_10x.csv')
sample_df <- read.csv(metacells_fn, stringsAsFactors = F, check.names = F)
print(dim(sample_df))
length(unique(sample_df$mouse_id))
sample_df <- sample_df %>%
  dplyr::select(mouse_id, batch_info) %>%
  dplyr::rename(batch=batch_info)
rownames(sample_df) <- sample_df$mouse_id
norm_sce$batch <- sample_df[norm_sce$sample,'batch']
unique(norm_sce$treatmentSt)
exclude_ts <- c('UU','UUU','UUUU','UUUUU')
unique(norm_sce1$treatmentSt)
norm_sce1 <- norm_sce[,!norm_sce$treatmentSt %in% exclude_ts]
sample_df$batch
View(sample_df)
# generate_input_data_ridge_regression(norm_sce, save_dir, datatag)
save_dir <- paste0(save_dir,'bbknn/')
convert_corrected_mtx_to_sce(norm_sce, save_dir, corrected_mtx_fn=NULL, datatag,
                             return_sce=F, save_data=T)
save_dir <- paste0(save_dir,'CPA/')
dir.create(save_dir)
meta_fn <- paste0(save_dir,datatag,'_meta_info_CPA.csv')
meta_info <- get_metadata_CPA(norm_sce, meta_fn)
meta_data <- data.table::fread(meta_fn) %>% as.data.frame()
dim(meta_data)
summary(as.factor(meta_data$treatmentSt))

summary(as.factor(meta_data$drug))
summary(as.factor(meta_data$drug_dose_name))
meta_info <- meta_data
library(stringr)
meta_info$dose <- ifelse(meta_info$drug=='CIS',str_count(meta_info$treatmentSt,'T'),
                         ifelse(meta_info$drug=='Vehicle',str_count(meta_info$treatmentSt,'U')-1,0))


meta_info$dose_val <- meta_info$dose/max(meta_info$dose)
# meta_info$dose_val <- ifelse(meta_info$drug=='Vehicle',1.0,meta_info$dose_val)
print(table(meta_info$treatmentSt, meta_info$dose))
print(summary(meta_info$dose_val))
# unique(meta_info$dose_val)
meta_info$cell_type <- 'TNBC'
meta_info$drug_dose_name <- paste0(meta_info$drug,'_',meta_info$dose_val)
unique(meta_info$drug_dose_name)
meta_info$drug_dose_name <- ifelse(meta_info$drug_dose_name=='Vehicle_1', 'Vehicle_1.0', meta_info$drug_dose_name)
meta_info$drug_dose_name <- ifelse(meta_info$drug_dose_name=='CIS_1', 'CIS_1.0', meta_info$drug_dose_name)


meta_info$cov_drug_dose_name <- paste0(meta_info$cell_type,'_',meta_info$drug_dose_name)
meta_info$condition <- meta_info$drug
meta_info$control <- ifelse(meta_info$drug=='Vehicle', 1, 0)
summary(as.factor(meta_info$control))
summary(as.factor(meta_info$condition))
summary(as.factor(meta_info$drug_dose_name))
summary(as.factor(meta_info$dose_val))
# meta_fn <- paste0(save_dir,'CPA/', datatag, "_meta_info_CPA.csv")
write.csv(meta_info, meta_fn, row.names = F, quote = F)
unique(meta_info$drug)
meta_info$cell_id
rownames(meta_info) <- meta_info$cell_id
baseline_cells <- meta_info %>%
  dplyr::filter(treatmentSt=='U')%>%
  dplyr::pull(cell_id)

length(baseline_cells)
cis_cells <- sample(baseline_cells, 200, replace=F)
untreated_cells <- baseline_cells[!baseline_cells %in% cis_cells]
length(cis_cells)
length(untreated_cells)

meta_info[cis_cells,'drug'] <- 'CIS'
meta_info$drug_dose_name <- ifelse(meta_info$drug_dose_name=='CIS_0','CIS_0.01',meta_info$drug_dose_name)
meta_info$drug_dose_name <- ifelse(meta_info$drug_dose_name=='Vehicle_0','Vehicle_0.01',meta_info$drug_dose_name)

meta_info$dose_val <- ifelse(meta_info$dose_val==0.00,0.01,meta_info$dose_val)
unique(meta_info$dose_val)
