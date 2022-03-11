

results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_Tyler_v3/'
datatag <- 'SA535'
library_grouping <- read.csv(paste0(results_dir,'library_groupings.csv'), check.names=F, stringsAsFactors=F)
# View(head(library_grouping))
dim(library_grouping)
length(unique(library_grouping$sample))
sample_ids <- unique(library_grouping$sample)
length(sample_ids)

results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/'
output_file <- paste0(results_10x_dir,'deg_analysis/de_output.png')
input_dir <- paste0(results_10x_dir,'clonealign')
input_file <- paste0(input_dir,'/normalized_sce.rds')
sample_id <- ''
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results'

norm_data <- logcounts(sce)
class(norm_data)
norm_data <- as.data.frame(as.matrix(norm_data))
dim(norm_data)
z_score <- scale(t(norm_data))
dim(z_score)



summary(as.double(z_score[1,]))
View(z_score[1:10,1:10])
z_score <- t(z_score)
min(z_score)
z_score_backup <- z_score
rownames(z_score)[1:3]
dim(z_score)

med_scores <- apply(obs_score, 1, median)
length(med_scores)

dim(stat)
View(head(stat))
stat$cnvA <- stat$cnvA$cnvA

signif_genes_df <- read.csv(paste0(save_dir, 'increase_genes_strict.csv'), check.names=F, stringsAsFactors=F)
dim(signif_genes_df)
signif_genes_df
stat <- stat[,colnames(stat) %in% c("cnvA","cnvB","scoreA","scoreB")]

stat_v2 <- stat[rownames(stat) %in% signif_genes_df$gene_ens,]
stat_v2 <- stat
dim(stat_v2)
rownames(stat_v2)[1:3]
dim(marker_df)
View(marker_df)
marker_df <- data.frame(avg_logFC=stat_v2$scoreA,
                        gene_ens=rownames(stat_v2))
genes_map_symb
library(dplyr)
marker_df <- marker_df %>% inner_join(genes_map_symb, by='gene_ens')
method_use <- 'z_scores'
base_name <- 'resistance_cells'
gmt_dir <- paste0(input_dir,"biodatabase/")
groups_use <- 'CloneA'
for(i in rep(1:length(gmt_ls),1)){
  # gmt_fn <- paste0(gmt_dir, gmt_ls[i])
  print(gmt_ls[i])
  print(pathway_names[i])
  pathway_ls <- get_pathway_results(marker_df, method_use,
                                    base_name, paste0(gmt_dir, gmt_ls[i]), 
                                    pathway_names[i], 
                                    groups_use,
                                    save_dir, 30)
}


df_cnv <- stat_v2
df_cnv$ensembl_gene_id <- rownames(df_cnv)

annots <- annotables::grch37

annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 

df_cnv <- inner_join(df_cnv, annots)

df_cnv$ensembl_gene_id
df_cnv <- dplyr::distinct(df_cnv, ensembl_gene_id, .keep_all = TRUE)

dim(df_cnv)
View(head(df_cnv))
rownames(df_cnv) <- df_cnv$ensembl_gene_id
df_genes <- df_cnv[,colnames(df_cnv) %in% c("scoreA", "scoreB")]

df_copy <- df_cnv[,colnames(df_cnv) %in% c("cnvA", "cnvB")]
colnames(df_copy)
summary(df_genes)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
hm_genes <- ComplexHeatmap::Heatmap(as.matrix(t(df_genes)), na_col = "white",
                              show_column_names=F,
                              show_row_names = T,
                              col = col_fun,
                              cluster_rows=F,cluster_columns=F,
                              column_split = df_cnv$chr,
                              name = paste0('Median genes z-score'),
                              column_names_gp = grid::gpar(fontsize = 10),
                              row_names_gp = grid::gpar(fontsize = 10),
                              show_heatmap_legend = T)

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
summary(as.factor(df_copy))
df_copy$cnvB
df_copy$cnvB <- ifelse(df_copy$cnvB>5,5,df_copy$cnvB)
df_copy$cnvA$cnvB <- ifelse(df_copy$cnvA$cnvB>5,5,df_copy$cnvA$cnvB)
col_cn = colorRamp2(c(0, 0, 5), c("green", "white", "purple"))
col_cn = colorRamp2(seq(0, 5, length = 3), c("blue", "#EEEEEE", "red"), 
                space = "RGB")
hm_cn <- ComplexHeatmap::Heatmap(as.matrix(t(df_copy$cnvA)), na_col = "white",
                              show_column_names=F,
                              show_row_names = T,
                              col = col_cn,
                              cluster_rows=F,cluster_columns=F,
                              column_split = df_cnv$chr,
                              name = paste0('Copy Number'),
                              column_names_gp = grid::gpar(fontsize = 10),
                              row_names_gp = grid::gpar(fontsize = 10),
                              show_heatmap_legend = T)

library(ComplexHeatmap)
hm <- hm_genes %v% hm_cn


png(paste0(save_dir,paste0("genes_copynumber_correlation.png")), 
    height = 800, width=2*1300,res = 2*72)
print(hm)
dev.off()









input_cnv_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/clonealign/whole_data/SA535_whole_data.csv'
df_cnv <- read.csv(input_cnv_fn, check.names=F, stringsAsFactors=F, header=T)
dim(df_cnv)
colnames(df_cnv)
View(head(df_cnv))

dim(obs_scoreA)
rownames(obs_scoreA)[1:3]
summary(as.numeric(obs_scoreA))
summary(as.numeric(obs_cnv_mat$A))

View(obs_scoreA['ENSG00000000457',1:10])
# A: resistant, B: sensitive
dim(obs_cnv_mat)
median_genesA
obs_genes <- rownames(obs_cnv_mat)
median_genesA <- median_genesA[obs_genes]
median_genesB <- median_genesB[obs_genes]


dim(obs_cnv_mat)
View(head(obs_cnv_mat))
amp_genes <- rownames(obs_cnv_mat)[(obs_cnv_mat$A>obs_cnv_mat$B)]  
# FALSE 854 TRUE 382, 72-1164 in total of 1236
del_genes <- rownames(obs_cnv_mat)[(obs_cnv_mat$A<obs_cnv_mat$B)] 
# FALSE 885 TRUE 404, 67 - 1222 >= in total of 1289
# Summary: 94% use unstrict condition >=, 31% use unstrict condition
length(amp_genes)
length(del_genes)
dim(obs_cnv_mat_2)
View(obs_cnv_mat[amp_genes[1:3],])

length(obs_scoreA_med)
length(obs_scoreB_med)
dim(obs_score)

lcA <- as.data.frame(as.matrix(logcounts(sceA)))
cells_use <- sample(colnames(sceB), ncol(sceA), replace=FALSE)
length(cells_use)
lcB <- as.data.frame(as.matrix(logcounts(sceB[,cells_use])))
dim(lcA)
dim(lcB)
View(lcB[1:3,1:3])


de_genes <- c()
obs_clone1 <- 'B'
obs_clone2 <- 'A'
# obs_genes <- amp_genes
obs_genes <- del_genes
for(g in obs_genes){
  
  # supA <- sum(as.numeric(obs_scoreA[g,])>0)
  # sdownA <- sum(as.numeric(obs_scoreA[g,])<0)
  # supB <- sum(as.numeric(obs_scoreB[g,])>0)
  # sdownB <- sum(as.numeric(obs_scoreB[g,])<0)
  # print(paste0('A_',supA,'_A_',sdownA,'  B_',supB,'_B_',sdownB))
  # if(supA > sdownA & (supA/n_cells > 0.5)){
  #   desc <- 'Up-reg'
  # } else if(sdownA > supA & (sdownA/n_cells > 0.5)){
  #   desc <- 'Down-reg'
  # } else{
  #   desc <- 'Unkn'
  # }
  # desc <- median_genesA[g] > median_genesB[g]
  t_test <- t.test(lcB[g,],lcA[g,])
  desc <- genes_stat[(genes_stat$ens_gene==g & genes_stat$clone==obs_clone1),'z_score'] > genes_stat[(genes_stat$ens_gene==g & genes_stat$clone==obs_clone2),'z_score']
  if(desc==TRUE & t_test$p.value < 0.05){
    desc==TRUE
  } else{
    desc==FALSE
  }
  de_genes <- c(de_genes, desc)
  
}

combine_stat_df <- data.frame(obs_gene=obs_genes,DE_gene=de_genes,
           row.names = obs_genes)
summary(as.factor(combine_stat_df$DE_gene))
combine_stat_df <- combine_stat_df[combine_stat_df$DE_gene==T,]

combine_stat_df_AU <- combine_stat_df
dim(combine_stat_df_AU)
View(head(combine_stat_df_AU))

combine_stat_df_DD <- combine_stat_df

summary_df <- data.frame(gene_ens=c(as.character(combine_stat_df_AU$obs_gene),as.character(combine_stat_df_DD$obs_gene)),
                         classify=c(rep('AU_A',rep(nrow(combine_stat_df_AU))),
                                    rep('AU_B',rep(nrow(combine_stat_df_DD)))))


View(head(summary_df))
dim(summary_df)
g <- 'ENSG00000004142'
obs_cnv_mat[g,]
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/'
genes_map_symb <- read.csv(paste0(input_dir, "biodatabase/meta_genes.csv"), header=T, stringsAsFactors=F)
rownames(genes_map_symb) <- genes_map_symb$gene_ens

summary_df <- inner_join(summary_df,genes_map_symb, by=('gene_ens'))

write.csv(summary_df,paste0(save_dir,'monotonic_increase_genes_strict.csv'), quote=F, row.names = F)





