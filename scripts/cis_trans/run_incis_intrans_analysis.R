input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# input_sdir <- ''
datatag <- 'SA1035'
# 
sample_df <- read.csv(paste0(input_dir,'SA1035_rna/snakemake_10x_SA1035/metasample_SA1035.csv'))
head(sample_df)



ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
cancer_ref_genes_df <- read.csv(paste0(ref_dif,'Behan_CFgenes.csv'), 
                                stringsAsFactors=F, check.names = F)
dim(cancer_ref_genes_df)
colnames(cancer_ref_genes_df)[which(colnames(cancer_ref_genes_df) == "ADAM PanCancer Core-Fitness genes")] <- "PanCancer_Fitness_genes"


output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/Mirela_output/gene_regression/SA609/'
datatag <- 'SA609'
signif_genes_df <- read.table(paste0(output_dir,datatag,'_significant'), header = T, sep='\t')
dim(signif_genes_df)
get_gene_type_stat_testing(signif_genes_df, datatag, output_dir, outlier_FC_thrs=3)
signif_genes_df$de_analysis <- paste0(signif_genes_df$sample1,'_vs_',signif_genes_df$sample2)
cancer_genes_ls <- list()
for(de in unique(signif_genes_df$de_analysis)){
  signif_genes <- read.csv(paste0(output_dir, datatag,'_result/',datatag,'_',de,'_signif_genes.csv'), check.names = F, stringsAsFactors=F)
  signif_genes <- signif_genes[signif_genes$is_fitness_gene,]
  print(dim(signif_genes))
  # fitness_genes <- summary(as.factor(signif_genes$classified_gene_dlp))
  nb_incis <- sum(!is.na(signif_genes$classified_gene_dlp)==TRUE)
  nb_intrans <- sum(is.na(signif_genes$classified_gene_dlp)==TRUE)
  # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
  cancer_genes_ls[[paste0(datatag,'_',de,'_incis')]] <- round(nb_incis*100/nrow(cancer_ref_genes_df),2)
  cancer_genes_ls[[paste0(datatag,'_',de,'_intrans')]] <- round(nb_intrans*100/nrow(cancer_ref_genes_df),2)
}

output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/Mirela_output/gene_regression/SA1035/'
datatag <- 'SA1035'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/clonealign/whole_data/'
signif_genes_df <- read.table(paste0(output_dir,datatag,'_significant'), header = T, sep='\t')
cnv_mat <- readRDS(paste0(dlp_dir,datatag,'_cnv_mat.rds'))
dim(signif_genes_df)
get_gene_type_stat(signif_genes_df, cnv_mat, datatag, output_dir, outlier_FC_thrs=3)

signif_genes_df$de_analysis <- paste0(signif_genes_df$sample1,'_vs_',signif_genes_df$sample2)

for(de in unique(signif_genes_df$de_analysis)){
  signif_genes <- read.csv(paste0(output_dir, datatag,'_result/',datatag,'_',de,'_signif_genes.csv'), check.names = F, stringsAsFactors=F)
  signif_genes <- signif_genes[signif_genes$is_fitness_gene,]
  print(dim(signif_genes))
  # fitness_genes <- summary(as.factor(signif_genes$classified_gene_dlp))
  nb_incis <- sum(!is.na(signif_genes$classified_gene_dlp)==TRUE)
  nb_intrans <- sum(is.na(signif_genes$classified_gene_dlp)==TRUE)
  # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
  cancer_genes_ls[[paste0(datatag,'_',de,'_incis')]] <- round(nb_incis*100/nrow(cancer_ref_genes_df),2)
  cancer_genes_ls[[paste0(datatag,'_',de,'_intrans')]] <- round(nb_intrans*100/nrow(cancer_ref_genes_df),2)
}


output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/Mirela_output/gene_regression/SA535/'
datatag <- 'SA535_cis'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/clonealign/whole_data/'
signif_genes_df <- read.table(paste0(output_dir,datatag,'_significant'), header = T, sep='\t')
cnv_mat <- readRDS(paste0(dlp_dir,'SA535_total_cnv_mat.rds'))
get_gene_type_stat(signif_genes_df, cnv_mat, datatag, output_dir, outlier_FC_thrs=3)

signif_genes_df$de_analysis <- paste0(signif_genes_df$sample1,'_vs_',signif_genes_df$sample2)
for(de in unique(signif_genes_df$de_analysis)){
  signif_genes <- read.csv(paste0(output_dir, datatag,'_result/',datatag,'_',de,'_signif_genes.csv'), check.names = F, stringsAsFactors=F)
  signif_genes <- signif_genes[signif_genes$is_fitness_gene,]
  print(dim(signif_genes))
  nb_incis <- sum(!is.na(signif_genes$classified_gene_dlp)==TRUE)
  nb_intrans <- sum(is.na(signif_genes$classified_gene_dlp)==TRUE)
  # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
  cancer_genes_ls[[paste0(datatag,'_',de,'_incis')]] <- round(nb_incis*100/nrow(cancer_ref_genes_df),2)
  cancer_genes_ls[[paste0(datatag,'_',de,'_intrans')]] <- round(nb_intrans*100/nrow(cancer_ref_genes_df),2)
}



output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/Mirela_output/gene_regression/SA535/'
datatag <- 'SA535_cx'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/clonealign/whole_data/'
signif_genes_df <- read.table(paste0(output_dir,datatag,'_significant'), header = T, sep='\t')
cnv_mat <- readRDS(paste0(dlp_dir,'SA535_total_cnv_mat.rds'))
get_gene_type_stat(signif_genes_df, cnv_mat, datatag, output_dir, outlier_FC_thrs=3)


signif_genes_df$de_analysis <- paste0(signif_genes_df$sample1,'_vs_',signif_genes_df$sample2)
for(de in unique(signif_genes_df$de_analysis)){
  signif_genes <- read.csv(paste0(output_dir, datatag,'_result/',datatag,'_',de,'_signif_genes.csv'), check.names = F, stringsAsFactors=F)
  signif_genes <- signif_genes[signif_genes$is_fitness_gene,]
  print(dim(signif_genes))
  nb_incis <- sum(!is.na(signif_genes$classified_gene_dlp)==TRUE)
  nb_intrans <- sum(is.na(signif_genes$classified_gene_dlp)==TRUE)
  # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
  cancer_genes_ls[[paste0(datatag,'_',de,'_incis')]] <- round(nb_incis*100/nrow(cancer_ref_genes_df),2)
  cancer_genes_ls[[paste0(datatag,'_',de,'_intrans')]] <- round(nb_intrans*100/nrow(cancer_ref_genes_df),2)
}

cancer_pct_df <- do.call(rbind, cancer_genes_ls)
cancer_pct_df <- as.data.frame(cancer_pct_df)
cancer_pct_df$series <- c(rep(c('SA609','SA1035','SA535_cis','SA535_cx'),c(8, 4, 12, 8)))
cancer_pct_df$series <- factor(cancer_pct_df$series, levels=c('SA609','SA1035','SA535_cis','SA535_cx'))

colnames(cancer_pct_df)[which(colnames(cancer_pct_df)=='V1')] <- 'Fitness_genes_pct'
rn <- as.character(rownames(cancer_pct_df))
cancer_pct_df$PDX <- factor(rn, levels = rn)
write.csv(cancer_pct_df, paste0(output_dir,'total_cancer_pct_df.csv'), quote=F, row.names = F)
cancer_pct_df <- read.csv(paste0(output_dir,'total_cancer_pct_df.csv'), check.names = F, stringsAsFactors = F)
View(cancer_pct_df)
cancer_pct_df$gene_type <- 'intrans'
cancer_pct_df$gene_type <- ifelse(grepl("*incis",cancer_pct_df$PDX),"incis",cancer_pct_df$gene_type)
# viz_pancancer_genes(cancer_pct_df, output_dir, tag="totaldata")

desc <- mclapply(strsplit(cancer_pct_df$PDX, "_"), function(x) {
  if(length(x)==6){
    return(paste0(x[3],'_',x[4],'_',x[5]))
  }else{
    return(paste0(x[2],'_',x[3],'_',x[4]))
  }
  
}, mc.cores = 2)
cancer_pct_df$de_pair <- as.character(desc)
# cancer_pct_df$PDX <- gsub('_incis', '', grep("*incis",cancer_pct_df$PDX, value=T))
# cancer_pct_df$PDX <- gsub('_intrans', '', grep("*intrans",cancer_pct_df$PDX, value=T))
colorcode <- c('#FF3232','#40A0E0')
ylabel <- 'Pct fitness genes (%)'
xlabel <- 'PDX - DE analysis'
plottitle <- 'Percentage genes in ADAM PanCancer Core-Fitness genes'

cancer_pct_df$de_pair <- factor(cancer_pct_df$de_pair, levels = unique(cancer_pct_df$de_pair))

cancer_pct_df$series <- factor(cancer_pct_df$series, levels = unique(cancer_pct_df$series))
plot_stack_barplot(cancer_pct_df, colorcode, xlabel, ylabel, plottitle, output_dir,'totaldata',
                   fa='gene_type', xa='de_pair', ya='Fitness_genes_pct')




input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/'
datatag <- 'SA609'
source("/home/htran/Projects/farhia_project/rnaseq/method_testing/in_cis_trans_utils.R")
sample_fn <- paste0(input_dir,'SA609_rna/snakemake_10x/library_10x.csv')
save_dir <- paste0(input_dir,'SA609_rna/snakemake_10x/')
save_dir <- '/home/htran/storage/python_workspace/velocity/result_SA609/'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/Mirela_output/gene_regression/SA609/'
load_sce_data(input_dir, save_dir, datatag, sample_fn, cutoff_zero_thrs=0.05)
load_clone_labels(sce_combine, input_dir, datatag, sample_fn, save_dir)

de_dir <- '/home/htran/Projects/farhia_project/drug_resistance/rnaseq/differential_expression/results/SA609-v6/comps/'
fns <- list.files(de_dir, pattern = '*_logfc_results.csv')
f <- fns[3]  #fns[3]
signif_genes_df <- data.table::fread(paste0(de_dir,f))

dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/added_segments/clonealign/whole_data/'
signif_genes_df <- read.table(paste0(output_dir,datatag,'_significant'), header = T, sep='\t')
cnv_mat <- readRDS(paste0(dlp_dir,datatag,'_cnv_mat.rds'))
dim(signif_genes_df)
View(head(signif_genes_df))
signif_genes_df$gene_id <- signif_genes_df$ensembl_gene_id
dim(cnv_mat)
View(head(cnv_mat))
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/testing/'
get_gene_type_stat(signif_genes_df, cnv_mat, datatag, output_dir, outlier_FC_thrs=3)

library(scater)
sample_df <- sample_df[,c("treatment_st","library_id")]
colData <- colData %>% left_join(sample_df, by=c('library_id'))
dim(colData)
View(sample_df)

write.csv(colData, paste0(output_dir,'colData.csv'), row.names = F, quote = F)

summary(as.factor(sce_rc$treatmentSt))

sample_df$library_id[sample_df$library_id=='TENX063'] <- 'SCRNA10X_SA_CHIP0063_000'
sum(unique(colData$library_id) %in% unique(sample_df$library_id))
t <- unique(colData$library_id) 
t[!t %in% unique(sample_df$library_id)]
length(unique(colData$library_id))
sce <- readRDS(paste0(input_dir,'total_sce_clones.rds'))
sce_r <- sce[,sce$clone=='R']
sce_rc <- sce[,sce$clone %in% c('R','C')]
dim(sce_rc)
sce_rc <- sce_rc[,!sce_rc$sample %in% c('SA609X4XB03084','SA609X5XB03235','SA609X7XB03573')]
unique(signif_genes_df$de_analysis)
de <- 'R_vs_C'
View(head(signif_genes_df))
signif_genes_df <- signif_genes_df[signif_genes_df$de_analysis==de,]
signif_genes_df <- signif_genes_df[abs(signif_genes_df$log2FoldChange)>0.25,]
dim(signif_genes_df)
View(head(signif_genes_df))

cut_off_overall <- 0.1
dim(sce_rc)
zero_cbind_r <- DelayedArray::rowMeans(assay(sce_rc[,sce_rc$clone=='R'], "counts") == 0)
obs_genes <- names(zero_cbind_r[zero_cbind_r <= (1 - cut_off_overall)])
length(obs_genes)
zero_cbind_c <- DelayedArray::rowMeans(assay(sce_rc[,sce_rc$clone=='C'], "counts") == 0)
obs_genes_c <- names(zero_cbind_c[zero_cbind_c <= (1 - cut_off_overall)])
length(obs_genes_c)
# genes_use <- intersect(obs_genes, signif_genes_df$gene_id)
genes_use <- intersect(obs_genes_c, genes_use)
length(genes_use)

signif_genes_df <- signif_genes_df[signif_genes_df$gene_id %in% genes_use,]
rownames(signif_genes_df) <- signif_genes_df$gene_id
tss <- c('U','UT','UUU','UUUU','UUUUU') 
mean_exp <- c()
for(t in tss){
  signif_genes_df$mean_C <- as.numeric(DelayedArray::rowMeans(assay(sce_rc[as.character(signif_genes_df$gene_id),
                                                                           sce_rc$clone=='C' & sce_rc$treatmentSt==t], "counts")))
  mean_exp <- c(mean_exp, sum(signif_genes_df$mean_C>=signif_genes_df$mean_R))
  
}
signif_genes_df$mean_C <- as.numeric(DelayedArray::rowMeans(assay(sce_combine[as.character(signif_genes_df$gene_id),
                                                                              sce_combine$clone=='C' & sce_combine$treatmentSt=='U'], "normcounts")))

signif_genes_df$mean_R <- as.numeric(DelayedArray::rowMeans(assay(sce_combine[as.character(signif_genes_df$gene_id),
                                                                              sce_combine$clone=='R' & sce_combine$treatmentSt %in% c('UTT','UTTT')], "normcounts")))

sum(signif_genes_df$mean_C>=signif_genes_df$mean_R)
nrow(signif_genes_df)
View(normcounts(sce_rc_test)[1:5,1:5])
test <- data.frame(cloneC_treatment_st=tss, 
                   cloneR_treatment_st=rep('UTT_UTTT',length(tss)), sensitive_genes=mean_exp, resistant_genes=nrow(signif_genes_df)-mean_exp)
View(test)
write.csv(test, paste0(output_dir,'R_vs_C.csv'), quote=F, row.names = F)
signif_genes_df$mean_R <- as.numeric(DelayedArray::rowMeans(assay(sce_rc[as.character(signif_genes_df$gene_id),
                                                                         sce_rc$clone=='R' & sce_rc$treatmentSt %in% c('UTT','UTTT')], "counts")))
View(signif_genes_df[1:50,c('log2FoldChange','mean_C','mean_R','sample1','sample2')])
colnames(signif_genes_df)

# normalized_mtx <- data.table::fread(paste0(output_dir,datatag,"_normalized_count"))

dim(normalized_mtx)
View(normalized_mtx[1:3,1:5])
sum(signif_genes_df$mean_C<signif_genes_df$mean_R)

sum(signif_genes_df$log2FoldChange>0)
sum(signif_genes_df$log2FoldChange<0)
sce_list <- list()
sce_rc_test <- sce_rc
sce_rc_test <- sce_rc_test[,sce_rc_test$treatmentSt %in% c('UTT','UTTT','U')]
dim(sce_rc_test)
for (l in unique(sce_rc_test$library_ids)){
  sce_normalized <- sce_normalize_size_factors(sce_rc_test[,sce_rc_test$library_ids==l], min_size=300, rlog=FALSE, exprs="counts")
  sce_list[[l]] <- sce_normalized
}
sce_combine <- sce_cbind_func(sce_list, rowData(sce_normalized), 
                              cut_off_overall = 0, 
                              exprs = c("counts", "normcounts"), 
                              colData_names = colnames(colData(sce_normalized)), meta_data=NULL)
sce_normalized_total <- sce_normalize_size_factors(sce_combine, 
                                                   min_size=300, rlog=F, exprs="normcounts")
dim(sce_normalized_total)

obs_genes <- signif_genes_df$gene_id[signif_genes_df$log2FoldChange >-0.55 & signif_genes_df$log2FoldChange < -0.5]
length(obs_genes)
dim(logcounts(sce_r))
signif_genes_df <- signif_genes_df[,signif_genes_df]
scater::plotExpression(sce_rc, obs_gene, x = "treatmentSt",colour_by="clone") #, 
rownames(sce_r)[1:4]

output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/Mirela_output/gene_regression/SA609/'
signif_genes_df <- read.table(paste0(output_dir,datatag,'_significant'), header = T, sep='\t')
sum(signif_genes_df$log2FoldChange>0)
View(head(signif_genes_df))
obs_gene <- 'ENSG00000187608'

meanC1 <- mean(logcounts(sce_rc[obs_gene,sce_rc$clone=='C' & sce_rc$treatmentSt=='U']))
meanC2 <- mean(logcounts(sce_rc[obs_gene,sce_rc$clone=='C' & sce_rc$treatmentSt=='UU']))
meanC2T <- mean(logcounts(sce_rc[obs_gene,sce_rc$clone=='C' & sce_rc$treatmentSt=='UT']))
meanC3 <- mean(logcounts(sce_rc[obs_gene,sce_rc$clone=='C' & sce_rc$treatmentSt=='UUU']))
meanC4 <- mean(logcounts(sce_rc[obs_gene,sce_rc$clone=='C' & sce_rc$treatmentSt=='UUUU']))
meanC5 <- mean(logcounts(sce_rc[obs_gene,sce_rc$clone=='C' & sce_rc$treatmentSt=='UUUUU']))


dim(meanC)
mean(meanC)
meanR4 <- mean(logcounts(sce_rc[obs_gene,(sce_rc$clone=='R') & (sce_rc$treatmentSt=='UTTTT')]))
meanR3 <- mean(logcounts(sce_rc[obs_gene,sce_rc$clone=='R' & sce_rc$treatmentSt=='UTTT']))
meanR2 <- mean(logcounts(sce_rc[obs_gene,sce_rc$clone=='R' 
                                & sce_rc$treatmentSt=='UTT']))

unique(sce_rc$library_ids)

# SA609_qc$library_id[SA609_qc$library_id=='TENX063'] <- 'SCRNA10X_SA_CHIP0063_000'
SA609_qc <- SA609_qc %>% inner_join(sample_df, by='library_id')
dim(SA609_qc)
View(SA609_qc)  
# ts <- grep('*T',SA609_qc$treatment_st, value=T)
# grep('*T',t1, value=T)
ts <- c('UT','UTT','UTTT','UTTTT','UTTTTT')
SA609_qc_treated <- SA609_qc[SA609_qc$treatment_st %in% ts,]
tss <- c('U','UU','UUU','UUUU','UUUUU')  
SA609_qc_untreated <- SA609_qc[SA609_qc$treatment_st %in% tss,]
View(SA609_qc_treated)
View(SA609_qc_untreated)

View(head(colData(sce)))
mean(sce$total_counts)
meta_cells <- colData(sce)

mean_reads <- c()
median_genes <- c()
lids <- intersect(unique(sce$library_ids), sample_df$library_id)
sample_df <- sample_df[sample_df$library_id %in% lids,]
sample_df_untreated <- sample_df[sample_df$treatment_st %in% tss,]
sample_df_treated <- sample_df[sample_df$treatment_st %in% ts,]
for(l in sample_df_untreated$library_id){
  median_genes <- c(median_genes, median(meta_cells[meta_cells$library_ids==l,'total_features_by_counts']))
  mean_reads <- c(mean_reads, mean(meta_cells[meta_cells$library_ids==l,'total_counts']))
}

for(l in sample_df_treated$library_id){
  median_genes <- c(median_genes, median(meta_cells[meta_cells$library_ids==l,'total_features_by_counts']))
  mean_reads <- c(mean_reads, mean(meta_cells[meta_cells$library_ids==l,'total_counts']))
}

qc_df <- data.frame(library_id=c(sample_df_untreated$library_id,sample_df_treated$library_id),
                    median_genes=round(median_genes,0),mean_reads=round(mean_reads,0))
View(qc_df)
SA609_qc$preservation_method_new
SA609_qc <- SA609_qc[,c('library_id','preservation_method_new')]
qc_df <- qc_df %>% left_join(sample_df, by=c('library_id'))
qc_df <- qc_df %>% left_join(SA609_qc, by=c('library_id'))

write.csv(qc_df, paste0(output_dir,'qc_SA609.csv'), quote = F, row.names = F)


signif_genes_df$de_analysis <- paste0(signif_genes_df$sample1,'_vs_',signif_genes_df$sample2)

de_pairs <- unique(signif_genes_df$de_analysis)
dim(marker_df)

View(head(signif_genes_df))
dim(signif_genes_df)
signif_genes_df <- signif_genes_df[signif_genes_df$qval<0.05 & abs(signif_genes_df$log2FoldChange)>0.25,]
# c('R','C')
# R vs C, 
# DEseq2: 6570 genes 5921 down-regulated 649 up-regulated
# Seurat: 2106 genes  884 down-regulated 1222 up-regulated
# Seurat: 1278 genes 856 down-regulated 422 up-regulated   threshold: pAdjustThrs: 0.05, logFC: 0.25
#   

# R-UTTT vs C-UUUU
# edgeR:  6839 genes, 2667 down-regulated, 4172 up-regulated threshold: FDR<0.01 & PValue<0.05, logFC: 0.25
# DEseq2: 6494 genes 5856 down-regulated, 638 up-regulated  pValue, pAdj<0.05, qVal<0.05, logFC > 0.25

# Seurat: 2757 genes, 706 down-regulated, 2051 up-regulated  all markers output
# Seurat: 1636 genes, 706 down-regulated, 930 up-regulated   threshold: pAdjustThrs: 0.05, logFC: 0.25

# intersect of DESeq2 and edgeR: 3637 genes, ~35%
# intersect of Seurat and DESeq2: 1825 genes in total of 2757 seurat genes
# intersect of Seurat and edgeR: 2701 genes in total of 2757 seurat genes

marker_df <- read.table(paste0(save_dir,'/deg_analysis/SA609_R_C/total_markers.txt'), sep='\t', header = T, row.names = 1)
View(head(marker_df))
dim(marker_df)
sum(marker_df$avg_logFC>0)
sum(marker_df$avg_logFC<=0)

sum(signif_genes_df$log2FoldChange>0)
sum(signif_genes_df$log2FoldChange<=0)

marker_df <- read.table(paste0(save_dir,'/deg_analysis/SA609_R_C/de_significant_genes.txt'), sep='\t', header = T, row.names = 1)

edgeR_df <- read.csv(paste0(save_dir,'/deg_analysis/SA609_R_C/SA609edgeR_significant_genes.csv'), check.names = F)
SA609_edgeR_significant_genes <- SA609_edgeR_significant_genes[abs(SA609_edgeR_significant_genes$logFC)>0.25,]
dim(SA609_edgeR_significant_genes)
View(head(SA609_edgeR_significant_genes))
sum(SA609_edgeR_significant_genes$logFC>0)
sum(SA609_edgeR_significant_genes$logFC<=0)
View(head(edgeR_df))
g <- intersect(edgeR_df$gene_symbol, signif_genes_df$gene_id)
length(g)
g <- intersect(marker_df$GENEID,signif_genes_df$gene_id)
edgeR_df$ensembl_gene_id[1:3]






# First create pairgroup metadata file

# data_dir <- '/home/htran/Projects/farhia_project/drug_resistance/cis_trans_landscape_treated_and_untreated/differential_expression/'
data_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/'

fns <- list.files(data_dir)

fns <- fns[grepl('logfc_results.csv', fns)]
pair_groups <- data.frame(result_fn=fns)
fns[1]
library(parallel)
get_info <- function(cell_ids, cores_use=2, idx=1) {
  labels <- mclapply2(strsplit(cell_ids, "_"), function(x) {
    return(x[idx])
  }, mc.cores = cores_use)
  return(as.character(labels))
}
library(snpEnrichment)
pair_groups$series <- get_info(pair_groups$result_fn, cores_use=2, idx=2)
pair_groups$title <- pair_groups$series
pair_groups$clone1 <- get_info(pair_groups$result_fn, cores_use=2, idx=6)
pair_groups$clone2 <- get_info(pair_groups$result_fn, cores_use=2, idx=8)
# View(pair_groups)
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna_cis/clonealign/'
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/clonealign/'
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/clonealign/'

data.table::fwrite(pair_groups, paste0(data_dir,'comparisons_drug_res_v2.csv'), quote=F)
# input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/cis_trans_landscape_treated_and_untreated/differential_expression/'
pair_groups <- data.table::fread(paste0(input_dir,'comparisons_drug_res.csv'))
pair_groups$datatag <- pair_groups$series  
pair_groups$desc <- gsub('_logfc_results.csv','', pair_groups$result_fn)
datatag <- 'SA609'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6'
# output_file <- paste0(input_dir,'/SA609_rna/deg_analysis/SA609_deg_pathway.rds')
# pair_groups_fn <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/deg_analysis/pair_groups_SA1035_SA609.csv')
# sce_fn <- paste0(results_10x_dir,'/total_sce_clones.rds')
# edgeR_DE_analysis_by_clones_samples(pair_groups_fn, sce_fn, output_file, datatag, 1500, 0.1)
# input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results','/',datatag,'_rna/deg_analysis/SA609-v6/')
# dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/added_segments/clonealign/whole_data/'
# cnv_mat <- readRDS(paste0(dlp_dir,datatag,'_cnv_mat.rds'))
# pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
source("/home/htran/Projects/farhia_project/rnaseq/cis_trans/in_cis_trans_utils.R")

# rm(cnv_mat)
datatag <- 'SA609'
dlp_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/clonealign/')
cnv_mat <- data.table::fread(paste0(dlp_dir,datatag,'_cisplatin_cnv_mat.csv')) %>% as.data.frame()
cnv_mat <- data.table::fread(paste0(dlp_dir,datatag,'_cnv_mat.csv')) %>% as.data.frame()

input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/cis_trans_landscape_treated_and_untreated/differential_expression/'
pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v2.csv')
save_dir <- paste0(dirname(dlp_dir),'/cis_trans/')
head(cnv_mat)
rownames(cnv_mat) <- cnv_mat$V1
cnv_mat$V1 <- NULL
# cnv_mat <- cnv_mat %>% 
#   tibble::column_to_rownames(var="ensembl_gene_id")
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
pair_groups$datatag <- pair_groups$series
pair_groups <- pair_groups[pair_groups$datatag==datatag,]
View(pair_groups)
# pair_groups <- pair_groups[4,]
# pair_groups$result_fn <- 'SA609_scran_SA609_UTTTT_R_UUUUU_H_logfc_results.csv'
# pair_groups$desc <- "scran_SA609_UTTTT_R_UUUUU_H"

# Scran + edgeR: abs(logFC)> 0.25, FDR < 0.01, pvalue < 0.05
# In_cis_Decrease_DownRegulated   In_cis_Decrease_UpRegulated 
# 296                            82 
# In_cis_Increase_DownRegulated   In_cis_Increase_UpRegulated 
# 73                           900 
# In_trans_DownRegulated          In_trans_UpRegulated 
# 1243                          2554 

# Scran + edgeR:: abs(logFC) > 0.5, FDR < 0.01, pvalue < 0.05
# In_cis_Decrease_DownRegulated   In_cis_Decrease_UpRegulated 
# 156                            27 
# In_cis_Increase_DownRegulated   In_cis_Increase_UpRegulated 
# 28                           550 
# In_trans_DownRegulated          In_trans_UpRegulated 
# 536                          1107 


# TMM + edgeR : abs(logFC)> 0.25, FDR < 0.01, pvalue < 0.05
# In_cis_Decrease_DownRegulated 467 
# In_cis_Decrease_UpRegulated 20 
# In_cis_Increase_DownRegulated 344 [PROBLEM]
# In_cis_Increase_UpRegulated 371 
# In_trans_DownRegulated 3725 [PROBLEM, TOO HIGH COMPARED TO UpRegulated]
# In_trans_UpRegulated 685

# TMM + edgeR : abs(logFC)> 0.5, FDR < 0.01, pvalue < 0.05
# In_cis_Decrease_DownRegulated 374 
# In_cis_Decrease_UpRegulated 9
# In_cis_Increase_DownRegulated 143  [PROBLEM]
# In_cis_Increase_UpRegulated  155  
# In_trans_DownRegulated 2109 [PROBLEM, TOO HIGH COMPARED TO UpRegulated]
# In_trans_UpRegulated 273

save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/cis_trans_scran/'
pair_groups$desc <- pair_groups$result_fn
input_dir <- data_dir
get_gene_type_edgeR_v2(pair_groups, cnv_mat, 
                       datatag, input_dir, save_dir,
                       cancer_ref_genes_fn=NULL, outlier_FC_thrs=4.25,
                       FDR_cutoff=0.01, minLogFC=0.5, pValueThrs=0.05)

obs_comps <- c('rstmm1_SA609_UT_R_UU_H','rstmm2_SA609_UTT_R_UUU_H','rstmm3_SA609_UTTT_R_UUUU_H','rstmm4_SA609_UTTTT_R_UUUUU_H','fitness31_SA609_UUUU_H_UUUU_C')
for(comp in obs_comps){
  de_genes <- data.table::fread(paste0(save_dir, comp,'/signif_genes_FDR0.01.csv')) %>% as.data.frame()
  # dim(de_genes)
  # summary(as.factor(de_genes$Gene_Type))
  de_genes <- de_genes %>%
    dplyr::filter(abs(logFC)>=0.5)
  cis_pos <- de_genes %>%
    dplyr::filter(Gene_Type %in% c('In_cis_Decrease_DownRegulated','In_cis_Increase_UpRegulated'))%>%
    dplyr::pull(ensembl_gene_id)
  cis_neg <- de_genes %>%
    dplyr::filter(Gene_Type %in% c('In_cis_Decrease_UpRegulated','In_cis_Increase_DownRegulated'))%>%
    dplyr::pull(ensembl_gene_id)
  nb_cis <- length(cis_pos) + length(cis_neg)
  print(paste0(comp,' :   pos(DD,IU): ',round(length(cis_pos)*100/nb_cis,1),'  neg(DU, ID): ',round(length(cis_neg)*100/nb_cis,1)))  
  
}
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/cis_trans/'
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/cis_trans/'

obs_comps <- list.dirs(save_dir, full.names = F)
obs_comps <- obs_comps[obs_comps!=""]
for(comp in obs_comps){
  de_genes <- data.table::fread(paste0(save_dir, comp,'/signif_genes_FDR0.01.csv')) %>% as.data.frame()
  # dim(de_genes)
  # summary(as.factor(de_genes$Gene_Type))
  de_genes <- de_genes %>%
    dplyr::filter(abs(logFC)>=0.5)
  cis_pos <- de_genes %>%
    dplyr::filter(Gene_Type %in% c('In_cis_Decrease_DownRegulated','In_cis_Increase_UpRegulated'))%>%
    dplyr::pull(ensembl_gene_id)
  cis_neg <- de_genes %>%
    dplyr::filter(Gene_Type %in% c('In_cis_Decrease_UpRegulated','In_cis_Increase_DownRegulated'))%>%
    dplyr::pull(ensembl_gene_id)
  nb_cis <- length(cis_pos) + length(cis_neg)
  print(paste0(comp,' :   pos(DD,IU): ',round(length(cis_pos)*100/nb_cis,1),'  neg(DU, ID): ',round(length(cis_neg)*100/nb_cis,1)))  
  
}


df <- cnv_mat[rowSums(is.na(cnv_mat))==0,]
dim(df)

# datatag <- 'SA1035'
# datatag <- 'SA535'
datatag <- 'SA609'
# dlp_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/clonealign/')
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
# cnv_mat <- data.table::fread(paste0(dlp_dir,datatag,'_cisplatin_cnv_mat.csv'))%>% as.data.frame()
cnv_mat <- data.table::fread(paste0(dlp_dir,datatag,'_cnv_mat.csv'))%>% as.data.frame()
# input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/cis_trans_landscape_treated_and_untreated/differential_expression/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'

pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')

# input_dir <- data_dir
# save_dir <- paste0(dirname(dlp_dir),'/cis_trans_v2/')
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/pathways_IU_DD/'
head(cnv_mat)
# cnv_mat <- cnv_mat %>% 
#   tibble::column_to_rownames(var="ensembl_gene_id")
rownames(cnv_mat) <- cnv_mat$V1
cnv_mat$V1 <- NULL
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
# pair_groups$datatag <- pair_groups$series
# pair_groups$desc <- pair_groups$result_fn
pair_groups <- pair_groups[pair_groups$datatag==datatag,]
print(pair_groups)
# View(pair_groups$clone1)


get_gene_type_edgeR_v2(pair_groups, cnv_mat, 
                       datatag, input_dir, save_dir,
                       cancer_ref_genes_fn=NULL,
                       FDR_cutoff=0.05, minLogFC=0.05, pValueThrs=0.05)


"/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/cis_trans/"
datatag <- 'SA535'
dlp_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/clonealign/')
cnv_mat <- data.table::fread(paste0(dlp_dir,datatag,'_cisplatin_cnv_mat.csv')) %>% as.data.frame()
input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/cis_trans_landscape_treated_and_untreated/differential_expression/'
pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
# save_dir <- dlp_dir
save_dir <- paste0(dirname(dlp_dir),'/cis_trans/')
rownames(cnv_mat) <- cnv_mat$V1
cnv_mat$V1 <- NULL
head(cnv_mat)
# cnv_mat <- cnv_mat %>% 
#   tibble::column_to_rownames(var="ensembl_gene_id")
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
pair_groups <- pair_groups[pair_groups$datatag==datatag,]

# pair_groups <- pair_groups %>%
#   dplyr::filter(!clone1 %in% c('A-G','A-E-J') & clone2!='B-C')
# Just keep it here to run program first
pair_groups$clone1 <- ifelse(pair_groups$clone1=='A-G','A',
                             ifelse(pair_groups$clone1=='A-E-J','A',pair_groups$clone1))

pair_groups$clone2 <- ifelse(pair_groups$clone2=='B-C','B',pair_groups$clone2)

View(pair_groups)
get_gene_type_edgeR_v2(pair_groups, cnv_mat, 
                       datatag, input_dir, save_dir,
                       cancer_ref_genes_fn=NULL, outlier_FC_thrs=4.25,
                       FDR_cutoff=0.01, minLogFC=0.25, pValueThrs=0.05)




signif_genes <- read.csv(paste0(input_dir,pair_groups$result_fn[1]), check.names=F, stringsAsFactors=F)

stat <- as.data.frame(stat_dlp, stringsAsFactors=F)
stat$PDX <- c(rep('SA609_CIS',nrow(dlp609)),rep('SA1035_CIS',nrow(dlp1035)),
              rep('SA535_CIS',nrow(dlpSA535_cis)),rep('SA535_CX5461',nrow(dlpSA535_cx)))
View(stat)







save_dir <- paste0(input_dir,'SA609_SA609X7XB03505-R_SA609X7XB03554-H/')
deg_df <- read.csv(paste0(save_dir,'signif_genes.csv'), check.names=F, stringsAsFactors = F)
dim(df)
deg_df <- deg_df[deg_df$ensembl_gene_id %in% df_cnv$ensembl_gene_id,]
colnames(df)
de <- 'SA609_SA609X7XB03505-R_SA609X7XB03554-H'
cnv_mat <- read.csv(paste0(input_dir,de,"/cnv_mat_total.csv"),check.names = F)
dim(cnv_mat)
View(head(cnv_mat))
cnv_mat <- cnv_mat[,colnames(cnv_mat) != 'classified_gene_dlp']
cnv_mat <- cnv_mat[,colnames(cnv_mat) %in% obs_clones]
cnv_mat$ensembl_gene_id <- rownames(cnv_mat)
cnv_mat$classified_gene_dlp
df_cnv <- cnv_mat %>%
  pivot_longer(!ensembl_gene_id, names_to = "cluster", values_to = "cnv")
dim(df_cnv)
View(head(df_cnv))
plot_DE_genes_edgeR(df, topGenes, capstr=" ", FDRcutoff=0.01, logFCcutoff=0.25, 
                    plttitle="A versus B", save_dir="",legendVisible=F)
obs_clones <- c('R','H')
cnv_mat <- cnv_mat[,colnames(cnv_mat) != 'classified_gene_dlp']
saveRDS(track_plot2,paste0(input_dir,de,"/track_plot_deg.rds"))
save_dir <- paste0(input_dir,de,"/")
cnv_mat_AB <- readRDS('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna_cys/deg_analysis/SA535_A_B/cnv_mat_AB.rds')
View(head(cnv_mat_AB))

df <- summary_incis_intrans_genes
dim(df)
df$nb_genes
length(unique(df$desc))
df <- df[,c("gene_type","nb_genes",'desc')]
df$pct_genes
library(dplyr)
library(tidyr)
df <- df %>%
  pivot_wider(names_from = gene_type, values_from = nb_genes)
View(head(df))
dim(df)

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/summary_results_DE_v6/')

write.csv(df, paste0(input_dir,'summary_incis_intrans_nbgenes.csv'), quote = F, row.names = F)




# [1] "rstmm1_SA609_UT_R_UU_H :   pos(DD,IU): 85.9  neg(DU, ID): 14.1"
# [1] "rstmm2_SA609_UTT_R_UUU_H :   pos(DD,IU): 88.9  neg(DU, ID): 11.1"
# [1] "rstmm3_SA609_UTTT_R_UUUU_H :   pos(DD,IU): 68.8  neg(DU, ID): 31.2"
# [1] "rstmm4_SA609_UTTTT_R_UUUUU_H :   pos(DD,IU): 77.7  neg(DU, ID): 22.3"
# [1] "fitness31_SA609_UUUU_H_UUUU_C :   pos(DD,IU): 59.5  neg(DU, ID): 40.5"

# [1] "fitness32_SA535_UUUUU_G_UUUUU_B_C :   pos(DD,IU): 93.6  neg(DU, ID): 6.4"
# [1] "rstmm10_SA535_UUTTTT_A_UUUUUU_G :   pos(DD,IU): 88.4  neg(DU, ID): 11.6"
# [1] "rstmm11_SA535_UUTTTTT_A_UUUUUU_G :   pos(DD,IU): 95.5  neg(DU, ID): 4.5"
# [1] "rstmm7_SA535_UUT_A_G_UUU_G :   pos(DD,IU): 27.2  neg(DU, ID): 72.8"  # not accurate cis, trans results due to combined clones
# [1] "rstmm8_SA535_UUTT_A_E_J_UUUU_G :   pos(DD,IU): 52  neg(DU, ID): 48" # not accurate cis, trans results due to combined clones
# [1] "rstmm9_SA535_UUTTT_A_UUUUU_G :   pos(DD,IU): 92.6  neg(DU, ID): 7.4"

# [1] "fitness33_SA1035_UUUU_E_UUUU_B :   pos(DD,IU): 99.1  neg(DU, ID): 0.9"
# [1] "rstmm5_SA1035_UTTT_H_UUUU_E :   pos(DD,IU): 88.4  neg(DU, ID): 11.6"
# [1] "rstmm6_SA1035_UTTTT_H_UUUUU_E :   pos(DD,IU): 86.4  neg(DU, ID): 13.6"


