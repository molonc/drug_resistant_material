# Viz Fig2,3 trackplot and cistrans plot
data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
input_dir <- paste0(data_dir,'signif_genes/')
save_dir <- paste0(data_dir,'signif_genes/')


## Cis trans plot
stat <- data.table::fread(paste0(save_dir,'wholedata_summary.csv')) %>% as.data.frame()
dim(stat)
used_comp <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_treated.csv')) %>% as.data.frame()
# comp_untreated_DH <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_untreated_DH.csv')) %>% as.data.frame()
used_comp1 <- used_comp %>%
  dplyr::filter(file_header!='scrande_SA1035_41') %>% # early time and do not contain cis genes. 
  dplyr::select(file_header, order)

stat <- stat %>%
  dplyr::filter(file_header %in% used_comp1$file_header) %>%
  dplyr::select(-order)%>% 
  left_join(used_comp1, by='file_header')
plot_ls <- plot_intrans_incis_prevalence_v4(stat, save_dir, verbose=F)



stat_cfcr <- data.table::fread(paste0(input_dir,'wholedata_CF_CR_sets.csv'), header=T) %>% as.data.frame()
stat_cfcr <- stat_cfcr %>%
  dplyr::filter(file_header %in% used_comp1$file_header) %>%
  dplyr::select(-order)%>% 
  left_join(used_comp1, by='file_header')
stat$ref_set
prop_plt_ls <- plot_reference_set_proportion_v2(stat_cfcr, save_dir, verbose=F, plttype = 'treated')
prop_plt_ls$prevalence_plt


pw_stat <- data.table::fread(paste0(input_dir, 'wholedata_pathways_CF_CR_sets.csv')) %>% as.data.frame()
dim(pw_stat)

used_comp1 <- used_comp %>%
  dplyr::filter(file_header!='scrande_SA1035_41') %>% # early time and do not contain cis genes. 
  dplyr::select(file_header, order, series)

pw_stat <- pw_stat %>%
  dplyr::filter(file_header %in% used_comp1$file_header) %>%
  # dplyr::select(-order)%>% 
  left_join(used_comp1, by='file_header')

pathway_plt_ls <- viz_reference_pathways(pw_stat, verbose=F, plttype = 'treated')
pathway_plt_ls$pathway_sig_plt
pathway_plt_ls$plg
p_total <- viz_Fig4_v2(plot_ls, prop_plt_ls, pathway_plt_ls, save_dir)


viz_Fig4_v2 <- function(plot_ls, prop_plt_ls, pathway_plt_ls, save_dir){
  rel_ht <- c(5,2,5)
  rx_lg <- cowplot::plot_grid(plot_ls$treated_SA609_c1lb, plot_ls$treated_SA535_c1lb,
                              plot_ls$treated_SA1035_c1lb, rel_heights = rel_ht, ncol=1)
  unrx_lg <- cowplot::plot_grid(plot_ls$treated_SA609_c2lb, plot_ls$treated_SA535_c2lb,
                                plot_ls$treated_SA1035_c2lb, rel_heights = rel_ht, ncol=1)
  
  passage_lg <- cowplot::plot_grid(plot_ls$treated_SA609_passage, plot_ls$treated_SA535_passage,
                                   plot_ls$treated_SA1035_passage, rel_heights = rel_ht, ncol=1)
  
  # prop_plt_ls, pathway_plt_ls
  plt_ct <- cowplot::plot_grid(passage_lg, rx_lg, unrx_lg,
                               plot_ls$treated,plot_ls$treated_proportion,plot_ls$treated_cistrans,
                               prop_plt_ls$prevalence_plt, pathway_plt_ls$pathway_sig_plt,
                               rel_widths = c(0.2,0.2,0.2,
                                              1,1,1,
                                              1,0.3), nrow=1)
  ggsave(paste0(save_dir,"Fig23_cistrans.png"),
         plot = plt_ct,
         height = 4,
         width = 7,
         # useDingbats=F,
         dpi=200)
  return(plt_ct)
  
  
}
