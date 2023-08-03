

base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
source(paste0(base_dir,'scripts/pipeline/utils/plot_utils.R'))
viz_medianCNV <- function(){
  clones <- c('A','G')
  base_name <- 'SA535'
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  save_dir <- paste0(base_dir,'materials/dlp_cnv/')
  save_fig_dir <- paste0(base_dir,'materials/dlp_cnv/trackplots/')
  # dir.create(save_fig_dir)
  df_cnv_fn <- paste0(save_dir, base_name, '_cnv_mat.csv.gz')
  
  # Get genes from filtered genes list, sctransform normalized output
  # meta_genes <- data.frame(ensembl_gene_id=rowData(sce)$ID, gene_symbol=rowData(sce)$Symbol, stringsAsFactors=F)
  # dim(meta_genes)
  # head(meta_genes)
  # data.table::fwrite(meta_genes, paste0(save_dir,base_name, '_meta_filtered_genes.csv'))
  # df_cnv <- data.table::fread(df_cnv_fn, check.names = F, stringsAsFactors = F) %>% as.data.frame()
  # t <- df_cnv %>%
  #   dplyr::filter(ensembl_gene_id=='ENSG00000136997')
  # t
  meta_genes <- data.table::fread(paste0(save_dir,base_name, '_meta_filtered_genes.csv')) %>% as.data.frame()
  dim(meta_genes)
  head(meta_genes)
  pcnv <- plot_CNV(df_cnv_fn, clones, meta_genes)
  saveRDS(pcnv,paste0(save_fig_dir,base_name, '_pcnv.rds'))
  
  pcnv <- readRDS(paste0(save_fig_dir,base_name, '_pcnv.rds'))
  pcnv$cnv_plot
  
  clones <- c('H','E')
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
  base_name <- 'SA1035'
  df_cnv_fn <- paste0(save_dir,base_name, '_cnv_mat.csv.gz')
  meta_genes <- data.table::fread(paste0(save_dir,base_name, '_meta_filtered_genes_10percent.csv')) %>% as.data.frame()
  dim(meta_genes)
  head(meta_genes)
  pcnv_1035 <- plot_CNV(df_cnv_fn, clones, meta_genes)
  pcnv_1035$cnv_plot
  dim(pcnv_1035$df_cnv)
  saveRDS(pcnv_1035, paste0(save_fig_dir,base_name, '_pcnv.rds'))
  pcnv_1035 <- readRDS(paste0(save_fig_dir,base_name, '_pcnv.rds'))
  
}



viz_DEG <- function(){
  ## meta info
  # data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
  # pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
  # pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  # pair_groups <- pair_groups %>%
  #   dplyr::filter(datatag=="SA535" & comp_type=='treated_vs_untreated')
  
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  input_dir <- paste0(base_dir,'materials/cis_trans/signif_genes/')
  save_dir <- paste0(base_dir,'materials/dlp_cnv/')
  save_fig_dir <- paste0(base_dir,'materials/dlp_cnv/trackplots/')
  
  datatag <- 'SA535'
  
  save_dir <- save_fig_dir
  
  # deg_fn <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UUTTTT_T_UUUUUU_J/signif_genes.csv')
  # deg_fn <- paste0(input_dir,datatag,'/SA535_2_SA535_UUTTTT_A_UUUUUU_G/signif_genes_FDR_0.01.csv')
  deg_fn <- paste0(input_dir,'scrande_SA535_3/signif_genes.csv')
  
  # additional_genes_fn <- paste0(base_dir, 'rnaseq_v6/trackplots_v3/SA535_ST_vs_J.csv')
  datatag <- 'Pt5'
  additional_genes_fn <- paste0(save_fig_dir, 'SA535_ST_vs_J_annotated_genes.csv')
  
  subtitle <- ' '
  
  track_incis <- plot_gex_cnv_v2(deg_fn=deg_fn,
                                 df_cnv = pcnv$df_cnv,
                                 clones=NULL,
                                 save_dir=save_fig_dir,
                                 sample_name=NULL,
                                 clone_str=desc,
                                 additional_genes_fn = additional_genes_fn,
                                 n_genes_to_annotate = 10, 
                                 plttitle=subtitle,
                                 gene_type='in-cis')
  # track_incis$track_plot
  subtitle <- 'Differentially expressed genes for Pt5: Rx X10:clone A vs. UnRx X9:clone G'
  track_intrans <- plot_gex_cnv_v2(deg_fn=deg_fn,
                                   df_cnv = pcnv$df_cnv,
                                   clones=NULL,
                                   save_dir=save_fig_dir,
                                   sample_name=NULL,
                                   clone_str=desc,
                                   additional_genes_fn = additional_genes_fn,
                                   n_genes_to_annotate = 10, 
                                   plttitle=subtitle,
                                   gene_type='in-trans')
  
  prop_cistrans_plt <- plot_cis_trans_prop(deg_fn, lg_pos="none")
  
  
  main_plot_535 <- cowplot::plot_grid(
    track_intrans$track_plot,
    track_incis$track_plot,
    prop_cistrans_plt$p ,
    pcnv$cnv_plot + theme(legend.position = 'none'),
    ncol = 1,
    rel_heights = c(1.65,1.6,0.8,1.1),
    align = 'v',
    labels = c('','','','','')
  )
  saveRDS(main_plot_535, paste0(save_fig_dir,datatag, '_trackplt.rds'))
  saveRDS(track_incis, paste0(save_fig_dir,datatag, '_cis_DE_trackplt.rds'))
  saveRDS(track_intrans, paste0(save_fig_dir,datatag, '_trans_DE_trackplt.rds'))
  saveRDS(prop_cistrans_plt, paste0(save_fig_dir,datatag, '_prop_cistrans_plt.rds'))
  datatag <- 'Pt5'
  main_plot_535 <- readRDS(paste0(save_fig_dir, datatag, '_trackplt.rds'))
  
  ### SA1035
  
  ## meta infos: 
  # data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
  # pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
  # pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  # pair_groups <- pair_groups %>%
  #   dplyr::filter(datatag=="SA1035" & comp_type=='treated_vs_untreated')
  # pair_groups
  
  input_dir <- paste0(base_dir,'materials/cis_trans/signif_genes/')
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  save_dir <- paste0(base_dir,'materials/dlp_cnv/')
  save_fig_dir <- paste0(base_dir,'materials/dlp_cnv/trackplots/')
  plottitle <- NULL
  datatag <- 'SA1035'
  save_dir <- save_fig_dir
  # deg_fn <- paste0(input_dir,datatag,'/SA1035_2_SA1035_UTTTT_H_UUUUU_E/signif_genes_FDR_0.01.csv')
  deg_fn <- paste0(input_dir,'scrande_SA1035_2/signif_genes.csv')
  additional_genes_fn <- paste0(save_fig_dir, 'SA1035_H_vs_E_annotated_genes.csv')
  
  
  # subtitle <- paste0(datatag,' cisplatin: in-cis DE genes')
  datatag <- 'Pt6'
  subtitle <- ' '
  track_incis_1035 <- plot_gex_cnv_v2(deg_fn=deg_fn,
                                 df_cnv = pcnv_1035$df_cnv,
                                 clones=NULL,
                                 save_dir=save_fig_dir,
                                 sample_name=NULL,
                                 clone_str=desc,
                                 additional_genes_fn = additional_genes_fn,
                                 n_genes_to_annotate = 10, 
                                 plttitle=subtitle,
                                 gene_type='in-cis')
  
  subtitle <- 'Differentially expressed genes for Pt6: Rx X8:clone H vs. UnRx X8:clone E'
  track_intrans_1035 <- plot_gex_cnv_v2(deg_fn=deg_fn,
                                   df_cnv = pcnv_1035$df_cnv,
                                   clones=NULL,
                                   save_dir=save_fig_dir,
                                   sample_name=NULL,
                                   clone_str=desc,
                                   additional_genes_fn = additional_genes_fn,
                                   n_genes_to_annotate = 10, 
                                   plttitle=subtitle,
                                   gene_type='in-trans')
  
  prop_cistrans_plt_1035 <- plot_cis_trans_prop(deg_fn, lg_pos="none")
  
  
  main_plot_1035 <- cowplot::plot_grid(
    track_intrans_1035$track_plot,
    track_incis_1035$track_plot,
    prop_cistrans_plt_1035$p ,
    pcnv_1035$cnv_plot + theme(legend.position = 'none'),
    ncol = 1,
    rel_heights = c(1.65,1.6,0.8,1.1),
    align = 'v',
    labels = c('','','','','')
  )
  saveRDS(main_plot_1035, paste0(save_fig_dir,datatag, '_trackplt.rds'))
  saveRDS(track_incis_1035, paste0(save_fig_dir,datatag, '_cis_DE_trackplt.rds'))
  saveRDS(track_intrans_1035, paste0(save_fig_dir,datatag, '_trans_DE_trackplt.rds'))
  saveRDS(prop_cistrans_plt_1035, paste0(save_fig_dir,datatag, '_prop_cistrans_plt.rds'))
  
  
  
  trackplts <- cowplot::plot_grid(main_plot_535, main_plot_1035, rel_heights = c(1,1), align = 'v', ncol=1)#labels=c('a','b')
  
  ggsave(paste0(save_dir,"SUPP_Fig5_trackplot_SA535_SA1035.png"),
         plot = trackplts,
         height = 12,
         width = 11,
         # useDingbats=F,
         dpi=150)
  
  ggsave(paste0(save_dir,"SUPP_Fig5_trackplot_SA535_SA1035.svg"),
         plot = trackplts,
         height = 12,
         width = 11,
         # useDingbats=F,
         dpi=250)
  
}





