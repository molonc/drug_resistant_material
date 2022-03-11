# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# BiocManager::install('Subread')
# library(Subread)

# devtools::install('/home/htran/storage/install_software/HGC')
# library(HGC)

source('/home/htran/Projects/farhia_project/rnaseq/bulk/tximport_utils.R')





bulk_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_bulk/output_v2/'
script_dir <- '/home/htran/Projects/farhia_project/rnaseq/bulk/'
input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_resistance/SA609_kallisto_Dec/'
save_dir <- bulk_dir
# source(paste0(script_dir,'plot_utils_bulk.R'))
source(paste0(script_dir,'tximport_utils.R'))

if(!dir.exists(bulk_dir)){
  dir.create(bulk_dir)
}
datatag <- 'SA609'
group_df <- data.table::fread(paste0(input_dir,'bulk_info_SA609.csv')) %>% as.data.frame()
group_df$sample_id
dim(group_df)
group_df <- group_df %>% 
  dplyr::filter(grepl('SA609',sample_id))
dim(group_df)
rownames(group_df) <- group_df$bulk_output_id
# group_df$bulk_alignment_kalliso <- 'DONE'
# group_df$sample_id[1]
# data.table::fwrite(group_df, paste0(input_dir,'bulk_info_SA609.csv'))

counts_df <- load_raw_counts_kallisto(sample_ids=group_df$bulk_output_id, 
                                      input_dir, datatag, save_dir)
dim(counts_df)
head(counts_df)
counts_df1 <- counts_df
counts_df1$ens_gene_id <- NULL
colnames(counts_df1) <- group_df[colnames(counts_df1),'sample_id']
colnames(counts_df1)
lcpm_dge <- load_dge_counts(as.matrix(counts_df1), save_dir, datatag, groups_use=NULL)
class(lcpm_dge)
head(lcpm_dge)
max(lcpm_dge)
min(lcpm_dge)

colSums(lcpm_dge)
lcpm_dge$ens_gene_id[1]
lcpm_dge$ens_gene_id <- get_ens_gene_ids(lcpm_dge$ens_gene_id)
annots <- annotables::grch38 %>%
  dplyr::select(ens_gene_id = ensgene, gene_symbol = symbol, chr)
annots <- annots[!duplicated(annots$ens_gene_id),]
lcpm_dge <- lcpm_dge %>% left_join(annots, by='ens_gene_id')
head(lcpm_dge)

lcpm_dge1 <- lcpm_dge %>%
  dplyr::filter(chr!='Y') #%>%
  # dplyr::select(-description)
dim(lcpm_dge1)

data.table::fwrite(lcpm_dge1, paste0(save_dir, datatag, '_lcpm_dge_to_Mirela.csv.gz'))

obs_genes <- c('TOP2A','COX6C','ID1','CCDC26','HIST1H4C')

lcpm_dge2 <- lcpm_dge %>%
  dplyr::filter(gene_symbol %in% obs_genes)
dim(lcpm_dge2)
sids <- colnames(lcpm_dge2)[grepl(datatag, colnames(lcpm_dge2))]
lcpm_dge2 <- lcpm_dge2 %>% 
  tidyr::pivot_longer(cols=all_of(sids), names_to = 'bulk_sample', values_to = 'expr')
lcpm_dge2 <- lcpm_dge2 %>% as.data.frame()
head(lcpm_dge2)
dim(lcpm_dge2)
class(lcpm_dge2)
meta_samples <- data.table::fread(paste0(save_dir, 'bulk_metadata_v2.csv')) %>% as.data.frame()
class(meta_samples)
sum(sids %in% meta_samples$sample)
lcpm_dge2$
lcpm_dge2 <- lcpm_dge2 %>% left_join(meta_samples, by=c('bulk_sample'='sample'))

plt_ls <- list()
for(g in unique(lcpm_dge2$gene_symbol)){
  df <- lcpm_dge2 %>%
    dplyr::filter(gene_symbol==g)
  p <- viz_lineplot_gene_expr(df, g)
  plt_ls[[g]] <- p
}
ptotal <- cowplot::plot_grid(plotlist = plt_ls, ncol=4, align = 'hv')
save_fig_dir <- paste0(save_dir,'signif_genes_bulk_expr/')
dir.create(save_fig_dir)
png(paste0(save_fig_dir,"bulk_genes_exp_",datatag,".png"), height = 2*150*2, width=2*250*5,res = 2*72)
print(ptotal)
dev.off() 

# Trajectory output
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/'
results_dir <- paste0(output_dir,'figs_v3/')

demo_genes <- obs_genes
plt_genes_demo <- list()
for(g in demo_genes){
  exg <- readRDS(paste0(results_dir,g,'_with_lg.rds'))  
  plt_genes_demo[[g]] <- exg
}
plt_combines <- list(plt_genes_demo$TOP2A,plt_ls$TOP2A,
                     plt_genes_demo$COX6C, plt_ls$COX6C,
                     plt_genes_demo$ID1, plt_ls$ID1,
                     plt_genes_demo$HIST1H4C, plt_ls$HIST1H4C)
ptotal <- cowplot::plot_grid(plotlist = plt_combines, ncol=2, align = 'v')
save_fig_dir <- paste0(save_dir,'signif_genes_bulk_expr/')
dir.create(save_fig_dir)
png(paste0(save_fig_dir,"bulk_genes_exp_",datatag,".png"), height = 2*180*4, 
    width=2*160*4,res = 2*72)
print(ptotal)
dev.off() 

## Older version
# results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_bulk/'
# input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_resistance/SA609_bulk/'
# base_name <- 'SA609'
# 
# results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_bulk/'
# input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_resistance/SA535_bulk/'
# base_name <- 'SA535'
# load_data(results_dir, input_dir, base_name)



# Select signif genes from trajectory

output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/tradeSeq_3000/'
patternRes <- readRDS(paste0(output_dir, "patternRes_out.rds"))

patternRes <- patternRes %>%
  dplyr::filter(pvalue<0.05)
dim(patternRes)

meta_genes <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/Symbol_ensembl.csv') %>% as.data.frame()
rownames(meta_genes) <- meta_genes$Ensembl
patternRes$Ensembl <- rownames(patternRes)
dim(patternRes)
patternRes <- patternRes %>% inner_join(meta_genes, by=c('Ensembl'))
oPat <- order(patternRes$waldStat, decreasing = TRUE)
print(head(patternRes[oPat,'Symbol']))
View(head(patternRes))
obs_genes <- patternRes[oPat[11:20],'Symbol']
lcpm_dge_extracted <- lcpm_dge %>%
  dplyr::filter(gene_symb %in% patternRes$Symbol) %>%
  as.data.frame()
dim(lcpm_dge_extracted)
plt_ls <- list()
for(g in obs_genes[1:10]){
  df <- lcpm_dge_extracted %>%
    dplyr::filter(gene_symb==g)
  p <- viz_lineplot_gene_expr(df, g)
  plt_ls[[g]] <- p
}
ptotal <- cowplot::plot_grid(plotlist = plt_ls, ncol=5, align = 'hv')
save_dir <- paste0(results_dir,'signif_genes_bulk_expr/')
dir.create(save_dir)
png(paste0(save_dir,"bulk_genes_exp_top11_20_",base_name,".png"), height = 2*150*2, width=2*250*5,res = 2*72)
print(ptotal)
dev.off() 


save_dir <- paste0(results_dir,'signif_genes_bulk_expr/')
dir.create(save_dir)
dim(lcpm_dge3)
lcpm_dge3 <- lcpm_dge2 %>%
  dplyr::filter(gene_symb %in% patternRes$Symbol) %>%
  as.data.frame()
data.table::fwrite(lcpm_dge3, paste0(results_dir,'bulk_CPM_significant_genes.csv.gz'))

lcpm_dge3 <- lcpm_dge2 %>%
  dplyr::filter(gene_symb %in% obs_genes) %>%
  as.data.frame()
for(g in obs_genes){
  df <- lcpm_dge3 %>%
    dplyr::filter(gene_symb==g)
  if(dim(df)[1]>0){
    p <- viz_lineplot_gene_expr(df, g)
    png(paste0(save_dir,g,".png"), height = 2*200, width=2*350,res = 2*72)
    print(p)
    dev.off() 
  }else{
    print(g)
  }
  
}



# Select signif genes from trajectory
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/'
patternRes <- readRDS(paste0(output_dir, "patternRes_out.rds"))
dim(patternRes)
patternRes <- patternRes %>%
  dplyr::filter(pvalue<0.05)
meta_genes <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/Symbol_ensembl.csv') %>% as.data.frame()
rownames(meta_genes) <- meta_genes$Ensembl
meta_genes$Symbol
patternRes$Ensembl <- rownames(patternRes)
dim(patternRes)
patternRes <- patternRes %>% inner_join(meta_genes, by=c('Ensembl'))
oPat <- order(patternRes$waldStat, decreasing = TRUE)
print(head(patternRes[oPat,'Symbol']))
View(head(patternRes))
obs_genes <- patternRes[oPat[1:50],'Symbol']

save_dir <- paste0(results_dir,'signif_genes_bulk_expr/')
dir.create(save_dir)
dim(lcpm_dge3)
lcpm_dge3 <- lcpm_dge2 %>%
  dplyr::filter(gene_symb %in% patternRes$Symbol) %>%
  as.data.frame()
data.table::fwrite(lcpm_dge3, paste0(results_dir,'bulk_CPM_significant_genes.csv.gz'))

lcpm_dge3 <- lcpm_dge2 %>%
  dplyr::filter(gene_symb %in% obs_genes) %>%
  as.data.frame()
for(g in obs_genes){
  df <- lcpm_dge3 %>%
    dplyr::filter(gene_symb==g)
  if(dim(df)[1]>0){
    p <- viz_lineplot_gene_expr(df, g)
    png(paste0(save_dir,g,".png"), height = 2*200, width=2*350,res = 2*72)
    print(p)
    dev.off() 
  }else{
    print(g)
  }
  
}

plt_ls <- list()
for(g in obs_genes[1:10]){
  df <- lcpm_dge3 %>%
    dplyr::filter(gene_symb==g)
  p <- viz_lineplot_gene_expr(df, g)
  plt_ls[[g]] <- p
}
ptotal <- cowplot::plot_grid(plotlist = plt_ls, ncol=5, align = 'hv')

png(paste0(save_dir,"bulk_genes_exp_top10_SA609.png"), height = 2*150*2, width=2*250*5,res = 2*72)
print(ptotal)
dev.off() 



# sample <- '3080'
samples <- c('3505','3510','3554')  #3080
# conditions <- c('Rx','Rx','UnRx')
groups_use <- c('Rx','UnRx')


samples <- c('3083','3080')
conditions <- c('Rx','UnRx')
datatag <- 'SA609'
counts_df <- load_counts_data(samples, input_dir, datatag)
dim(counts_df)

t <- colSums(counts_df)
meta_coverage <- data.frame(sample=paste0(datatag,'_',names(t)), coverage=as.numeric(t))
data.table::fwrite(meta_coverage, paste0(input_dir,datatag,'_coverage.csv'))
View(head(counts_df))
meta_df <- data.frame(sample=samples, condition=conditions)
rownames(meta_df) <- meta_df$sample
meta_df$condition <- ifelse(meta_df$condition==groups_use[1],paste0("2_",meta_df$condition),
                            paste0("1_",meta_df$condition))
groups_use <- unique(meta_df$condition)

meta_df
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# genes_symb_df <- read.csv(paste0(base_dir,'biodatabase/meta_genes.csv'), check.names = F, stringsAsFactors = F)
# dim(genes_symb_df)
# rownames(genes_symb_df) <- genes_symb_df$gene_ens
# meta_genes <- data.frame(gene_symb=rownames(counts_df),stringsAsFactors=F)
# meta_genes <- meta_genes %>% left_join(genes_symb_df,by=c('gene_symb'))
# dim(meta_genes)
# meta_genes <- meta_genes[!duplicated(meta_genes$gene_symb),]

dge <- load_dge_counts(counts_df, meta_df$condition)

# predefined_disp <- 0.01342844
save_dir_de <- paste0(base_dir,'SA609_bulk/',paste(samples,collapse = '_'),'/')
dir.create(save_dir_de, recursive = T)
de_genes <- edgeR_de_v2(dge, meta_df, save_dir_de, tag=NULL)
de_genes$logFC

dge_late <- readRDS(paste0(save_dir_de, 'edgeR_dge2_Rx_1_UnRx.rds'))
predefined_disp <- dge_late$common.dispersion
dim(dge_late$counts)
View(head(dge_late$counts))
counts_df <- dge_late$counts
dim(de.com)
sum(de.com$table$PValue<0.05)
View(head(de.com$table))
tag <- NULL
saveRDS(de.com,paste0(save_dir_de, 'edgeR_dge',tag,'_',predefined_disp,'.rds'))
de_genes <- de.com$table
write.csv(de_genes, file=paste0(save_dir_de, 'edgeR_significant_genes.csv'), row.names=FALSE, quote=FALSE)



input_dir <- '/home/htran/storage/datasets/drug_resistance_RNAseq/SA609_bulk/'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_bulk/'

input_dir <- '/home/htran/storage/datasets/drug_resistance_RNAseq/SA535_bulk/'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_bulk/'

# sample <- '3080'
# samples <- c('3505','3510','3554')  #3080
# conditions <- c('Rx','Rx','UnRx')
sample_df <- data.table::fread(paste0(base_dir, 'bulk_metadata.csv'), header=T) %>% as.data.frame()
View(sample_df)
# sample_df$condition <- ifelse(grepl(groups_use[2], sample_df$treatmentSt),paste0("1_",groups_use[2]),
#                               paste0("2_",groups_use[1]))

meta_df <- sample_df %>%
  dplyr::filter(series %in% c('SA535_untreated','SA535_cis') &
                  time %in% c('X8','X9'))

meta_df <- sample_df %>%
  dplyr::filter(series %in% c('SA535_untreated','SA535_cis') &
                  time %in% c('X6','X8'))

meta_df <- sample_df %>%
  dplyr::filter(series %in% c('SA535_untreated','SA535_cx') &
                  time %in% c('X6','X8'))

meta_df <- sample_df %>%
  dplyr::filter(series %in% c('SA535_untreated','SA535_cx') &
                  time %in% c('X6','X8'))


meta_df <- sample_df %>%
  dplyr::filter(series %in% c('SA535_untreated','SA535_cx') &
                  sample != 'SA535X8XB03548')

meta_df <- sample_df %>%
  dplyr::filter(series %in% c('SA609') &
                  time %in% c('X7'))

datatag <- 'SA609'
datatag <- 'SA535 Cisplatin'

datatag <- 'SA535 CX5461'
rm(res)
res <- run_edgeR(datatag, meta_df, input_dir, base_dir)
viz_heatmap_rawdata(res$counts_df, res$meta_df, res$save_dir, res$datatag)
de_genes_fn <- paste0(res$save_dir, 'edgeR_significant_genes_2_Rx_1_UnRx.csv')
viz_heatmap_DE_genes(res$counts_df, res$meta_df, res$save_dir, res$datatag, de_genes_fn)


counts_df <- res$counts_df
meta_df <- res$meta_df
save_dir <- res$save_dir
datatag <- 'SA535 CX5461'

# 1 sample vs. 1 sample
dge_late <- readRDS(paste0(save_dir_de, 'edgeR_dge2_Rx_1_UnRx.rds'))

predefined_disp <- dge_late$common.dispersion
predefined_disp
groups_use <- c('Rx','UnRx')
sample_df <- data.table::fread(paste0(base_dir, 'bulk_metadata.csv'), header=T) %>% as.data.frame()
View(sample_df)

meta_df <- sample_df %>%
  dplyr::filter(series %in% c('SA535_untreated','SA535_cis') &
                  time %in% c('X8','X6'))
meta_df$sample
meta_df$sample <- c('3099','3101','3664') # using substring to get 4 final characters
rownames(meta_df) <- meta_df$sample
groups_use <- unique(meta_df$condition)
datatag <- 'SA535_cisplatin'
counts_df <- load_counts_data(meta_df$sample, input_dir, datatag)
head(counts_df)
rm(dge)
dge <- load_dge_counts(counts_df, meta_df$condition)
save_dir_de <- paste0(base_dir,paste(meta_df$sample,collapse = '_'),'/')
dir.create(save_dir_de, recursive = T)
de_genes <- edgeR_de_v2(dge, meta_df, save_dir_de, tag=NULL)



# dim(de.com)
# sum(de.com$table$PValue<0.05)
# View(head(de.com$table))
# tag <- NULL
# saveRDS(de.com,paste0(save_dir_de, 'edgeR_dge',tag,'_',predefined_disp,'.rds'))
# de_genes <- de.com$table
# write.csv(de_genes, file=paste0(save_dir_de, 'edgeR_significant_genes.csv'), row.names=FALSE, quote=FALSE)



# dds <- DESeqDataSetFromMatrix(countData = counts_df, colData = meta_df, design = ~ condition)
# # #----- estimate dispersions -----#
# dds <- estimateSizeFactors(dds, type = "poscounts")
# cores_use <- 3
# nbcores <- detectCores()
# if(cores_use > nbcores){
#   cores_use <- nbcores
# }
# p <- MulticoreParam(cores_use)
# dds_rg <- DESeq(dds, fitType = "local", parallel = TRUE, BPPARAM = p)  # TO DO: parallel computation


# Run edgeR get output, do correlation 
