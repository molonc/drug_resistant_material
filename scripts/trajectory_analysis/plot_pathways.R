source('/home/htran/Projects/farhia_project/rnaseq/pipeline/utils/pathway_utils.R')


save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/pathway_results/'
desc='total_genes_modules'
datatag <- 'SA609'
gsea_609 <- data.table::fread(paste0(save_dir, "signf_pathways_",datatag,"_",desc,".csv")) %>% as.data.frame()
dim(gsea_609)
# View(gsea_609)

save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/tradeSeq_3000/pathway_results/'
datatag <- 'SA535'
gsea_535 <- data.table::fread(paste0(save_dir, "signf_pathways_",datatag,"_",desc,".csv")) %>% as.data.frame()
dim(gsea_535)


save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/tradeseq_v2/results/pathway_results/'
datatag <- 'SA1035'
gsea_1035 <- data.table::fread(paste0(save_dir, "signf_pathways_",datatag,"_",desc,".csv")) %>% as.data.frame()
dim(gsea_1035)

# TO DO: check genes fall into the signif list first
pathway_stat <- dplyr::bind_rows(gsea_609, gsea_535, gsea_1035)



pw_ls <- list()  
for(obs_pw in unique(pathway_stat$pathway)){
  pw_tmp <- pathway_stat %>%
    dplyr::filter(pathway==obs_pw)
  genes_set <- list()  
  for(dt in pw_tmp$datatag){
    tmp <- pw_tmp %>%
      dplyr::filter(datatag==dt)
    genes_set[[dt]] <- unlist(strsplit(tmp$signf_genes,','))
  }
  
  pw <- viz_vennDiagram_ggplot(genes_set, save_fig_dir, 'total', 
                               plottitle=obs_pw, genes_set_names=NULL,
                               ht=600, wd=600)
  pw_ls[[obs_pw]] <- pw
  
}
p1 <- cowplot::plot_grid(p, NULL, rel_widths = c(2,1))
ptotal <- cowplot::plot_grid(p1, pw_ls$HALLMARK_E2F_TARGETS, 
                             pw_ls$HALLMARK_G2M_CHECKPOINT, pw_ls$HALLMARK_MITOTIC_SPINDLE, 
                             ncol=2)

png(paste0(save_fig_dir,"trajectory_pathways_summary.png"), 
    height = 2*1000, width=2*800, res = 2*72)
print(ptotal)
dev.off()


save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/'
genes_df <- data.table::fread(paste0(save_dir, 'SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(genes_df)
unique(genes_df$gene_type)
datatag <- 'SA609'
genes_df$gene_type <- ifelse(genes_df$gene_type_module %in% c('Module2','Module3'),'Module2,3',
                             ifelse(genes_df$gene_type_module %in% c('Module5','Module6'),'Module5,6',genes_df$gene_type_module))

unique(genes_df$gene_type)
# stat_609 <- get_bootstrap_stat_each_series(genes_df, paste0(save_dir,'pathway_results/'), datatag)

## gprofiler version
gmt_id <- "gp__tSF3_MdRf_Z9o"
stat_609 <- get_gprofiler_pathways(genes_df, save_dir, datatag, 
                      custom_id=gmt_id, pathway_fn=NULL, save_data=F)

save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/tradeSeq_3000/'
genes_df <- data.table::fread(paste0(save_dir, 'SA535_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(genes_df)
datatag <- 'SA535'
# stat_535 <- get_bootstrap_stat_each_series(genes_df, paste0(save_dir,'pathway_results/'), datatag)


## gprofiler version
gmt_id <- "gp__tSF3_MdRf_Z9o"
stat_535 <- get_gprofiler_pathways(genes_df, save_dir, datatag, 
                                   custom_id=gmt_id, pathway_fn=NULL, save_data=F)


## Bootstrap test
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/tradeseq_v2/results/'
# genes_df <- data.table::fread(paste0(save_dir, 'SA1035_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
# dim(genes_df)
# datatag <- 'SA1035'
# stat_1035 <- get_bootstrap_stat_each_series(genes_df, paste0(save_dir,'pathway_results/'), datatag)

## gprofiler statistical test version
gmt_id <- "gp__tSF3_MdRf_Z9o"
stat_1035 <- get_gprofiler_pathways(genes_df, save_dir, datatag, 
                                   custom_id=gmt_id, pathway_fn=NULL, save_data=F)


pw_stats_modules <- dplyr::bind_rows(list(stat_609, stat_535, stat_1035))
save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
data.table::fwrite(pw_stats_modules, paste0(save_fig_dir,'pw_stats_modules_3series.csv'))
data.table::fwrite(pw_stats_modules, paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv'))

pw_stats_modules <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series.csv'))
pw_stats_modules <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv'))
pw_SA609 <- pw_stats_modules %>%
  dplyr::filter(datatag=='SA609')
dim(pw_SA609)


save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv')) %>% as.data.frame()
pathway_stat <- pathway_stat %>%
  # dplyr::filter(datatag=='SA609') %>%
  dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
pathway_stat$gene_type_module <- gsub('Module','M',pathway_stat$gene_type_module)
my_font <- "Helvetica"
pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
pathway_stat1 <- pathway_stat %>%
  dplyr::filter(grepl('SA609',datatag))
# pathway_stat1$datatag <- 'Pt4'
# pathway_stat1$datatag <- paste0(pathway_stat1$datatag, '_', pathway_stat1$gene_type_module)
# pathway_stat1$datatag <- gsub('Pt4_',pathway_stat1$datatag)
pathway_stat1$datatag <- pathway_stat1$gene_type_module
datatag <- 'SA609'
# pw_609 <- viz_pathways(pathway_stat1, datatag, save_fig_dir)
pw_609 <- viz_pathways_dotplot(pathway_stat1, datatag, save_fig_dir)
pw_609

colnames(pathway_stat1)

pathway_stat1 <- pathway_stat1 %>%
  dplyr::select(gene_type_module, pathway)

pathway_stat1$gene_type_module
pathway_stat1 <- pathway_stat1 %>%
  dplyr::left_join(meta_gm_SA609, by=c('gene_type_module')) %>%
  dplyr::rename(gene_type_origin=gene_type_module, gene_type=gm_manuscript_lb)

View(pathway_stat1)

save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv')) %>% as.data.frame()
pathway_stat <- pathway_stat %>%
  # dplyr::filter(datatag=='SA609') %>%
  dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
pathway_stat$gene_type_module <- gsub('Module','M',pathway_stat$gene_type_module)
my_font <- "Helvetica"
pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
pathway_stat1 <- pathway_stat %>%
  dplyr::filter(grepl('SA535',datatag))
# pathway_stat1$datatag <- 'Pt4'
# pathway_stat1$datatag <- paste0(pathway_stat1$datatag, '_', pathway_stat1$gene_type_module)
# pathway_stat1$datatag <- gsub('Pt4_',pathway_stat1$datatag)
pathway_stat1$datatag <- pathway_stat1$gene_type_module
datatag <- 'SA535'
# pw_535 <- viz_pathways(pathway_stat1, datatag, save_fig_dir)
pw_535 <- viz_pathways_dotplot(pathway_stat1, datatag, save_fig_dir)
pw_535
# unique(pathway_stat1$datatag)


save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv')) %>% as.data.frame()
pathway_stat <- pathway_stat %>%
  # dplyr::filter(datatag=='SA609') %>%
  dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
pathway_stat$gene_type_module <- gsub('Module','M',pathway_stat$gene_type_module)
my_font <- "Helvetica"
pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
pathway_stat1 <- pathway_stat %>%
  dplyr::filter(grepl('SA1035',datatag))
# pathway_stat1$datatag <- 'Pt4'
# pathway_stat1$datatag <- paste0(pathway_stat1$datatag, '_', pathway_stat1$gene_type_module)
# pathway_stat1$datatag <- gsub('Pt4_',pathway_stat1$datatag)
datatag <- 'SA1035'
pathway_stat1$datatag <- pathway_stat1$gene_type_module
# pw_1035 <- viz_pathways(pathway_stat1, datatag, save_fig_dir)
pw_1035 <- viz_pathways_dotplot(pathway_stat1, datatag, save_fig_dir)


## New labels 
# save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
# # New labels for gene modules
# meta_gm_SA609 <- data.frame(gene_type_module=c('Module1','Module2','Module3','Module4','Module5','Module6'),
#                             gm_manuscript_lb=c('Module1','Module3','Module5','Module4','Module6','Module2'))
# data.table::fwrite(meta_gm_SA609, paste0(save_fig_dir, 'meta_gene_module_labels_SA609.csv'))
# 
# meta_gm_SA535 <- data.frame(gene_type_module=c('Module1','Module2','Module3','Module4','Module5'),
#                             gm_manuscript_lb=c('Module2','Module1','Module3','Module4','Module5'))
# data.table::fwrite(meta_gm_SA535, paste0(save_fig_dir, 'meta_gene_module_labels_SA535.csv'))



input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
# input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/'
df <- data.table::fread(paste0(input_dir, 'SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz'))
dim(df)
summary(as.factor(df$gene_type_module))


library(gprofiler2)

# Testing pathway analysis
obs_module <- "Module5"
genes_M4 <- df %>%
  dplyr::filter(gene_type_module==obs_module)
dim(genes_M4)
head(genes_M4)
genes_use <- unique(genes_M4$gene_symbol)
length(genes_use)
# genes_use <- unique(genes_M4$ens_gene_id)

cat(genes_use, sep = '\n')

library(gprofiler2)
pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt'
# pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/c2.cp.kegg.v7.1.symbols.gmt'
custom_id <- gprofiler2::upload_GMT_file(pathway_fn)

# pathway_fn: a hallmark genes symbol gmt file
custom_id <- gprofiler2::upload_GMT_file(pathway_fn)
#gp__6TUj_rpgo_WS8 is custom id return from web after uploading
## fdr option
gostres_fdr <- gprofiler2::gost(list(genes_use), organism = custom_id,
                            correction_method='fdr')
## "gSCS" (synonyms: "analytical", "g_SCS") option
gostres_gSCS <- gprofiler2::gost(list(genes_use), organism = custom_id,
                            correction_method='gSCS')
## bonferroni option
gostres_bonferroni <- gprofiler2::gost(list(genes_use), organism = custom_id,
                            correction_method='bonferroni')

gostres <- gprofiler2::gost(genes_use, organism = custom_id,significant = F,
                                correction_method = 'gSCS')


stat <- gostres$result
cols_use <- c('p_value','intersection_size','term_id')
stat <- stat %>%
  dplyr::select(all_of(cols_use)) %>%
  dplyr::filter(p_value<0.5) %>%
  dplyr::rename(reference_set=term_id, nb_signif_genes=intersection_size)
stat
# c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
#   "analytical")

p_pathway <- cowplot::plot_grid(pw_609, pw_535)
p_pathway
save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
df <- data.table::fread(paste0(save_fig_dir, "Gene_Module_Chromatin_Status.csv")) %>% as.data.frame()
dim(df)
head(df)

df$chromatin_status <- gsub(' Promoters','',df$`Promoter Status`)
df$Module <- paste0('M',df$Module)
summary(as.factor(df$chromatin_status))
tag <- 'SA609'
df1 <- df %>%
  dplyr::filter(datatag==tag)
chr_609 <- viz_chromatin_barplot(df1, tag, save_fig_dir)
chr_609$p
chr_609$plg

df$Module <- paste0('Module',df$Module)
meta_gm_SA609$gene_type_module
df <- df %>%
    dplyr::left_join(meta_gm_SA609, by=c('Module'='gene_type_module')) %>%
    dplyr::rename(gene_type_origin=Module, gene_type=gm_manuscript_lb)
data.table::fwrite(df, paste0(save_fig_dir, "Gene_Module_Chromatin_Status_v2.csv"))

View(meta_gm_SA609)
tag <- 'SA535'
df1 <- df %>%
  dplyr::filter(datatag==tag)
chr_535 <- viz_chromatin_barplot(df1, tag, save_fig_dir)
chr_535$p


tag <- 'SA1035'
df1 <- df %>%
  dplyr::filter(datatag==tag)
chr_1035 <- viz_chromatin_barplot(df1, tag, save_fig_dir)
chr_1035$p

df1 <- df %>%
  dplyr::filter(datatag!='SA1035')

chr_plt <- viz_chromatin_barplot(df1, tag, save_fig_dir)
chr_plt$p
p609 <- cowplot::plot_grid(pw_609, chr_609$p, nrow=2, align='v',rel_heights = c(6,1))
p535 <- cowplot::plot_grid(pw_535, chr_535$p, nrow=2, align='v',rel_heights = c(6,1))
ptotal <- cowplot::plot_grid(p609, p535, nrow=1)
ptotal

p_chrom <- ptotal

plegends <-  cowplot::plot_grid(hm_lg, hm_plg_cistrans, pw_total$plg, chr_plt$plg, nrow=4) + theme(plot.background = element_rect(fill = "white", colour = "white"))
plegends
pathway_all <- cowplot::plot_grid(p_chrom, plegends, rel_widths = c(4.2,1), align = 'h')

save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv')) %>% as.data.frame()
pathway_stat <- pathway_stat %>%
  # dplyr::filter(datatag=='SA609') %>%
  dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
pathway_stat$gene_type_module <- gsub('Module','M',pathway_stat$gene_type_module)
my_font <- "Helvetica"
pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
pathway_stat1 <- pathway_stat %>%
  dplyr::filter(grepl('SA535|SA609',datatag))
# View(head(pathway_stat1))
summary(as.factor(pathway_stat1$datatag))
colnames(pathway_stat1)
empty_col <- tibble::tibble(p_value='NULL', nb_signf_genes=0,precision='', recall='',pathway='',
                            gene_type_module=c('M5'),signif_genes='', datatag=c('SA609'))
dim(pathway_stat1)
pathway_stat1 <- rbind(pathway_stat1,empty_col)

# table(pathway_stat1$gene_type_module, pathway_stat1$datatag)
pw_total <- viz_pathways_combined_dotplot(pathway_stat1, 'SA609_SA535', save_fig_dir)
pw_total$p
dev.off()


p_chrom <- cowplot::plot_grid(pw_total$p, chr_plt$p, ncol=1, align = 'v', rel_heights = c(5,1)
                              )#label=c('e','f')
p_chrom
# p_chrom <- ptotal

plegends <-  cowplot::plot_grid(hm_lg, hm_plg_cistrans, pw_total$plg, chr_plt$plg, nrow=4) + theme(plot.background = element_rect(fill = "white", colour = "white"))
plegends
pathway_all <- cowplot::plot_grid(p_chrom, plegends, rel_widths = c(4.2,1), align = 'h')





library(dplyr)
library(ggplot2)
save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv')) %>% as.data.frame()
pathway_stat <- pathway_stat %>%
  # dplyr::filter(datatag=='SA609') %>%
  dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
pathway_stat$gene_type_module <- gsub('Module','M',pathway_stat$gene_type_module)
my_font <- "Helvetica"
pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
pathway_stat1 <- pathway_stat %>%
  dplyr::filter(grepl('SA609',datatag) & p_value <0.05)

head(pathway_stat1)
dim(pathway_stat1)

## TO DO: add new gene module names
p <- viz_pathways_barplot(pathway_stat1)
p


ggsave(paste0(save_fig_dir,"Fig7_trajectory_SA609_pathways.png"),
       plot = p,
       height = 4,
       width = 3.5,
       # useDingbats=F,
       dpi=150)

pathway_stat2 <- pathway_stat %>%
  dplyr::filter(grepl('SA535',datatag) & p_value <0.05)

head(pathway_stat2)
dim(pathway_stat2)
p <- viz_pathways_barplot(pathway_stat2)
ggsave(paste0(save_fig_dir,"Fig7_trajectory_SA535_pathways.png"),
       plot = p,
       height = 4,
       width = 3.5,
       # useDingbats=F,
       dpi=150)

