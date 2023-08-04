tag <- 'SA609'
# output_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones_v2/'
output_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones_v3/'

fn <- paste0(output_dir,'mapped_clones_',tag,'.csv.gz')
mapped_clones <- data.table::fread(fn) %>% as.data.frame()
dim(mapped_clones)
mapped_clones <- mapped_clones %>%
    dplyr::filter(pct_pure>=0.5 & avg_mappability>=0.9)

mapped_clones <- mapped_clones %>%
  dplyr::select(-avg_mappability, -pct_pure) %>%
  tidyr::pivot_wider(names_from='clone_id', values_from='cn_median')

head(mapped_clones)
dim(mapped_clones)
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat_testing.csv.gz')
data.table::fwrite(mapped_clones, cnv_fn)

cnv_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/clonealign/SA535_cnv_mat.csv'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')

data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v4.csv')
save_dir <- paste0(data_dir,'signif_genes/')
# datatag <- 'SA609'
datatag <- 'SA535'
# cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
res <- load_data(pair_groups_fn, datatag, cnv_fn)

pair_groups <- res$pair_groups
# idx <- 'scrande_SA609_10'
idx <- 'scrande_SA535_4'
pair_groups <- pair_groups %>%
  dplyr::filter(file_header==idx)
pair_groups
head(res$cnv_mat)
cnv_mat <- res$cnv_mat
# rownames(cnv_mat) <- cnv_mat$V1
# cnv_mat$V1 <- NULL
# head(cnv_mat)
df <- get_gene_type_edgeR_v2(pair_groups, cnv_mat, 
                       datatag, input_dir, save_dir,
                       FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)

df <- df %>%
  dplyr::filter(abs(logFC)>0.5 & FDR<0.01 & PValue<0.05)
summary(as.factor(df$Gene_Type))
sum(grepl('In_cis',df$Gene_Type))
(86 +3)/131
(36+6)/131
## DEBUG

# idx <- 'scrande_SA609_10'
# head(res$cnv_mat)
# fn <- paste0(save_dir,idx,'/signif_genes.csv')
# 
# df <- data.table::fread(fn) %>% as.data.frame()
# dim(df)
# View(head(df))
# df <- df %>%
#   dplyr::filter(abs(logFC)>0.5 & FDR<0.01 & PValue<0.05)
# summary(as.factor(df$Gene_Type))


# Do DE analysis
# Comparing with output from edgeR logFC positive or negative? 
# How many DE genes
# CNV output: using bin regions, or segment regions? 

# datatag <- 'SA609'
# base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_normalized.rds'))
# dim(sce)
# sce1 <- sce[,sce$clone %in% c('H','C')]
# dim(sce1)
# summary(as.factor(sce$clone))
# 
# sids <- unique(sce1$sample)

# SA535: 
# old version:
#   In_cis_Decrease_DownRegulated 
# 57 
# In_cis_Decrease_UpRegulated 
# 288 
# In_cis_Increase_DownRegulated 
# 10 
# In_cis_Increase_UpRegulated 
# 92 
# In_trans_DownRegulated 
# 328 
# In_trans_UpRegulated 
# 1686 
# unmapped 
# 373 


# In_cis_Decrease_DownRegulated 
# 97 
# In_cis_Decrease_UpRegulated 
# 393 
# In_cis_Increase_DownRegulated 
# 15 
# In_cis_Increase_UpRegulated 
# 116 
# In_trans_DownRegulated 
# 354 
# In_trans_UpRegulated 
# 1792 
# unmapped 
# 67
