

# TO DO: get results of batch effect models for all series of data 
# check SA535 cisplatin output
# check cx output 
# do 1 correction for 2 drugs only 
# sa1035 output 
# SA609 output 


source(paste0("/home/htran/Projects/farhia_project/rnaseq/pipeline/utils/normalize_utils.R"))


datatag <- 'SA535'
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/SA535_cisplatin/')
# save_dir <- paste0(save_dir,'bbknn/')
output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
norm_sce <- readRDS(paste0(output_dir, datatag,'_total_sctransform_normalized.rds'))
# dim(norm_sce)
# meta_info <- read.csv(paste0(save_dir,datatag, "_meta_info.csv"))
# dim(meta_info)
# View(head(meta_info))
# generate_input_data_ridge_regression(norm_sce, save_dir, datatag)
# convert_corrected_mtx_to_sce(norm_sce, save_dir, corrected_mtx_fn=NULL, datatag)
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/correction_output/')
corrected_fn <- paste0(save_dir, datatag, "_norm_batch_corrected_sce.rds")
corrected_sce <- readRDS(corrected_fn)
dim(corrected_sce)
saveRDS(corrected_sce, file=corrected_fn)



# adata.obs['dose_val'] = adata.obs.dose.astype(float) / np.max(adata.obs.dose.astype(float))
# adata.obs['dose_val'][adata.obs['drug'].str.contains('Vehicle')] = 1.0
# 
# adata.obs['cell_type'] = 'A549'
# adata.obs['drug_dose_name'] = adata.obs.drug.astype(str) + '_' + adata.obs.dose_val.astype(str)
# adata.obs['cov_drug_dose_name'] = adata.obs.cell_type.astype(str) + '_' + adata.obs.drug_dose_name.astype(str)
# adata.obs['condition'] = adata.obs.drug.copy()
# adata.obs['control'] = [1 if x == 'Vehicle_1.0' else 0 for x in adata.obs.drug_dose_name.values]
# rm(norm_sce)
# rm(corrected_data)
datatag <- 'SA1035'
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
# output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
# norm_sce <- readRDS(paste0(output_dir, datatag,'_sctransform_normalized.rds'))
# dim(norm_sce)
# generate_input_data_ridge_regression(norm_sce, save_dir, datatag)
# convert_corrected_mtx_to_sce(norm_sce, save_dir, corrected_mtx_fn=NULL, datatag,
#                              return_sce=F, save_data=T)
corrected_fn <- paste0(save_dir, datatag, "_norm_batch_corrected_sce.rds")
corrected_sce <- readRDS(corrected_fn)



# rm(norm_sce)
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
# unique(norm_sce$treatmentSt)
# exclude_ts <- c('UU','UUU','UUUU','UUUUU')
# unique(norm_sce1$treatmentSt)
# norm_sce1 <- norm_sce[,!norm_sce$treatmentSt %in% exclude_ts]
# sample_df$batch
# View(sample_df)
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


  
meta_info <- data.table::fread(meta_fn)
dim(meta_info)
meta_info$cell_id[1]
paste0(save_dir, datatag, "_norm_batch_corrected_sce.rds")
meta_info$sample_id <- get_sample_id(meta_info$cell_id)
meta_info$treatmentSt <- sample_df[meta_info$sample_id,'treatmentSt']
meta_info$sample_id <- paste0(meta_info$sample_id,'_',meta_info$treatmentSt)
library(dplyr)
meta_info <- meta_info %>%
  dplyr::filter(clone=='R')
t <- table(meta_info$clone, meta_info$sample_id)
View(t)
get_sample_id <- function(cell_id) {
  labels <- sapply(strsplit(cell_id, "_"), function(x) {
    return(x[1])
  })
  return(as.character(labels))
}
# meta_info <- read.csv("/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/SA609/SA609_meta_info.csv")
# unique(meta_info$clone)
  