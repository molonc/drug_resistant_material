library(dplyr)
library(tidyr)
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
cnv_fn <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones_v2/SA609_mapped_genes/mapped_wholedata_SA609.csv.gz'
p1 <- 'scrande_SA609_3_SA609_UTTT_R_UUUU_H_logfc_results.csv'
fn <- paste0(input_dir, p1)

df <- data.table::fread(fn) %>% as.data.frame()
dim(df)
df <- df %>%
  dplyr::filter(abs(logFC)>0.5 & FDR < 0.01)

cnv <- data.table::fread(cnv_fn) %>% as.data.frame()
dim(cnv)
cnv <- cnv %>% drop_na()

sum(df$ensembl_gene_id %in% cnv$ensembl_gene_id)

df <- df %>%
  dplyr::filter(ensembl_gene_id %in% cnv$ensembl_gene_id)
dim(df)
cnv <- cnv %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(avg_mappability_R=mean(avg_mappability_R),
                   pct_pure_R=mean(pct_pure_R),
                   avg_mappability_H=mean(avg_mappability_H),
                   pct_pure_H=mean(pct_pure_H),
                   cn_median_R=mean(cn_median_R),
                   cn_median_H=mean(cn_median_H)
                   ) %>%
  dplyr::filter(ensembl_gene_id %in% df$ensembl_gene_id)
dim(cnv)
length(unique(cnv$ensembl_gene_id))
View(cnv)
colnames(cnv)
map_thrs <- 0.9
pct_pure_thrs <- 0.6

dim(t)
cols_use <- colnames(cnv)
cols_use <- cols_use[!grepl('_H', cols_use) & !grepl('clone_id', cols_use)]
cnv <- cnv %>%
  dplyr::select(all_of(cols_use))

df <- df %>% left_join(cnv, by='ensembl_gene_id')
dim(df)
View(head(cnv))
length(unique(df$ensembl_gene_id))
t <- df %>%
  dplyr::filter(avg_mappability_R>=map_thrs &
                  # pct_pure_R >= pct_pure_thrs & 
                  avg_mappability_H>=map_thrs #&pct_pure_H >= pct_pure_thrs
                  )

t1 <- df %>%
  dplyr::filter(avg_mappability_R>=map_thrs &
                  pct_pure_R < pct_pure_thrs & 
                  avg_mappability_H>=map_thrs & pct_pure_H < pct_pure_thrs
  )
dim(t1)
t2 <- df %>%
  dplyr::filter(avg_mappability_R>=map_thrs &
                  pct_pure_R >= pct_pure_thrs & 
                  avg_mappability_H>=map_thrs & pct_pure_H < pct_pure_thrs
  )
t3 <- df %>%
  dplyr::filter(avg_mappability_R>=map_thrs &
                  pct_pure_R < pct_pure_thrs & 
                  avg_mappability_H>=map_thrs & pct_pure_H >= pct_pure_thrs
  )
dim(t2)
dim(t3)
dim(df)
# sum(cnv$avg_mappability_R)