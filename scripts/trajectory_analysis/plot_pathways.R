save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/pathway_results/'
desc='total_genes_modules'
datatag <- 'SA609'
gsea_609 <- data.table::fread(paste0(save_dir, "signf_pathways_",datatag,"_",desc,".csv")) %>% as.data.frame()
dim(gsea_609)

head(gsea_609)
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
datatag <- 'SA609'
genes_df$gene_type <- ifelse(genes_df$gene_type_module %in% c('Module2','Module3'),'Module2,3',
                             ifelse(genes_df$gene_type_module %in% c('Module5','Module6'),'Module5,6',genes_df$gene_type_module))

unique(genes_df$gene_type)
stat_609 <- get_bootstrap_stat_each_series(genes_df, paste0(save_dir,'pathway_results/'), datatag)

## gprofiler version
gmt_id <- "gp__tSF3_MdRf_Z9o"
stat_609 <- get_gprofiler_pathways(genes_df, save_dir, datatag, 
                      custom_id=gmt_id, pathway_fn=NULL, save_data=F)

save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/tradeSeq_3000/'
genes_df <- data.table::fread(paste0(save_dir, 'SA535_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(genes_df)
datatag <- 'SA535'
stat_535 <- get_bootstrap_stat_each_series(genes_df, paste0(save_dir,'pathway_results/'), datatag)


## gprofiler version
gmt_id <- "gp__tSF3_MdRf_Z9o"
stat_535 <- get_gprofiler_pathways(genes_df, save_dir, datatag, 
                                   custom_id=gmt_id, pathway_fn=NULL, save_data=F)


save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/tradeseq_v2/results/'
genes_df <- data.table::fread(paste0(save_dir, 'SA1035_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(genes_df)
datatag <- 'SA1035'
stat_1035 <- get_bootstrap_stat_each_series(genes_df, paste0(save_dir,'pathway_results/'), datatag)

## gprofiler version
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



library(gprofiler2)

# Testing pathway analysis
obs_module <- "Module4"
genes_M4 <- total_genes %>%
  dplyr::filter(gene_type_module=="Module4")
dim(genes_M4)
head(genes_M4)
genes_use <- unique(genes_M4$gene_symbol)



# pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt'
# custom_id <- gprofiler2::upload_GMT_file(pathway_fn)
# library(gprofiler2)
# # pathway_fn: a hallmark genes symbol gmt file
# custom_id <- gprofiler2::upload_GMT_file(pathway_fn)
# #gp__6TUj_rpgo_WS8 is custom id return from web after uploading
# ## fdr option
# gostres <- gprofiler2::gost(list(genes_use), organism = 'gp__6TUj_rpgo_WS8', 
#                             correction_method='fdr')
# ## "gSCS" (synonyms: "analytical", "g_SCS") option
# gostres <- gprofiler2::gost(list(genes_use), organism = 'gp__6TUj_rpgo_WS8', 
#                             correction_method='gSCS')
# ## bonferroni option
# gostres <- gprofiler2::gost(list(genes_use), organism = 'gp__6TUj_rpgo_WS8', 
#                             correction_method='bonferroni')
