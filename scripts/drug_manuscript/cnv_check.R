

library(dplyr)
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/clonealign/sohrab_cnv/Mirela_cnv/'

# SA1035X8XB03425 T
# SA1035X8XB03631 U

sH <- read.table(paste0(base_dir,'SA1035X8XB03425.tsv'), sep='\t', header = T, check.names = F, stringsAsFactors = F)
sE <- read.table(paste0(base_dir,'SA1035X8XB03631.tsv'), sep='\t', header = T, check.names = F, stringsAsFactors = F)

sH <- sH %>%
  dplyr::filter(cluster=='H')
dim(sH)

sE <- sE %>%
  dplyr::filter(cluster=='E')
dim(sE)

colnames(sH)
genes_use <- intersect(sH$ensembl_gene_id, sE$ensembl_gene_id)
length(genes_use)

sH <- sH %>%
  dplyr::filter(ensembl_gene_id %in% genes_use) %>%
  dplyr::select(median_cnmode,ensembl_gene_id)
rownames(sH) <- sH$ensembl_gene_id
sH <- sH[genes_use,]
sE <- sE %>%
  dplyr::filter(ensembl_gene_id %in% genes_use) %>%
  dplyr::select(median_cnmode,ensembl_gene_id)

rownames(sE) <- sE$ensembl_gene_id
sE <- sE[genes_use,]
head(sE)
cnv <- data.frame(H=sH$median_cnmode,E=sE$median_cnmode, row.names = genes_use, stringsAsFactors = F)
dim(cnv)
head(cnv)

cnv <- cnv[(cnv$H-cnv$E)!=0,]
cnv_global <- readRDS('/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/clonealign/whole_data/SA1035_cnv_mat.rds')
dim(cnv_global)
head(cnv_global)
cnv_global <- cnv_global %>%
  dplyr::select(H,E)

cnv_global <- cnv_global[(cnv_global$H-cnv_global$E)!=0,]
g_global <- rownames(cnv_global)
g_sample <- rownames(cnv)
length(g_global)
length(g_sample)
head(g_sample)
length(g_global[g_global %in% g_sample])
length(g_sample[g_sample %in% g_global])


cnv_mat_wd <- read.csv(paste0(save_dir,base_name,'_cnv_mat.csv'), check.names = F, stringsAsFactors=F)

head(cnv_mat_wd)
cnv_mat_wd <- cnv_mat_wd[,colnames(cnv_mat_wd) %in% c(clones,"ensembl_gene_id")]
cnv_mat_wd <- cnv_mat_wd[cnv_mat_wd$E-cnv_mat_wd$H!=0,]
dim(cnv_mat_wd)
var_genes <- cnv_mat_wd$ensembl_gene_id

deg_df$ensembl_gene_id
cnv_mat$ensembl_gene_id <- rownames(cnv_mat)
deg_df <- deg_df %>% left_join(cnv_mat, by=c("ensembl_gene_id"))
head(deg_df)
# deg_df <- deg_df %>%
#          dplyr::select(-H.x,-H.y,-E.x,-)
deg_df$is_incis <- ifelse(deg_df$ensembl_gene_id %in% cnv_mat$ensembl_gene_id,T,F)
summary(as.factor(deg_df$is_incis))
var_genes <- intersect(deg_df$ensembl_gene_id, cnv_mat$ensembl_gene_id)
length(var_genes)

df_cnv$chr <- as.character(df_cnv$chr)
df_cnv <- df_cnv[df_cnv$chr==2,]
df_cnvH <- df_cnv %>%
  dplyr::filter(clone=='H')
dim(df_cnvH)
df_cnvE <- df_cnv %>%
  dplyr::filter(clone=='E')
dim(df_cnvE)
sum(is.na(df_cnvH$cnv))
df_cnv1 <- df_cnv %>%
  pivot_wider(names_from = clone, values_from = cnv)
df_cnvH <- as.data.frame(df_cnvH)
df_cnvE <- as.data.frame(df_cnvE)
rownames(df_cnvH) <- df_cnvH$ensembl_gene_id
rownames(df_cnvE) <- df_cnvE$ensembl_gene_id
df_cnvH <- df_cnvH[!duplicated(df_cnvH$ensembl_gene_id),]
df_cnvE <- df_cnvE[!duplicated(df_cnvE$ensembl_gene_id),]
genes <- intersect(df_cnvH$ensembl_gene_id, df_cnvE$ensembl_gene_id)
length(genes)
df_cnvE <- df_cnvE[genes,]
df_cnvH <- df_cnvH[genes,]
cnv_test <- data.frame(H=df_cnvH$cnv,E=df_cnvE$cnv, row.names = genes, stringsAsFactors = F)
dim(cnv_test)
head(cnv_test)
cnv_test <- cnv_test[as.numeric(cnv_test$H)-as.numeric(cnv_test$E)!=0,]

# H E
# ENSG00000243147 3 2

df_cnv_check <- df_cnv[df_cnv$ensembl_gene_id=='ENSG00000243147',]
df_cnv_check



write.csv(df_cnv, paste0(save_fig_dir,'SA1035_H_E_cnv_mat.csv'), row.names = F, quote=F)



