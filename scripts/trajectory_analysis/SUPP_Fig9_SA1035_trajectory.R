
script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/trajectory_analysis/'
source(paste0(script_dir, "slingshot_utils.R"))
source(paste0(script_dir, "tradeseq_utils.R"))

## Heatmap average expression plot
## Extracting heatmap of average genes expression x lineages 
plot_heatmap <- function(){
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  datatag <- 'SA1035'
  genes_df <- data.table::fread(paste0(save_dir, datatag,'_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
  dim(genes_df)
  colnames(genes_df)
  # summary(as.factor(genes_df$gene_type_module))
  genes_df <- genes_df %>%
    dplyr::select(-gene_type) %>%
    dplyr::rename(gene_type=gene_type_module)
  summary(as.factor(genes_df$gene_type))
  
  plttitle <- 'Activated Repressed Transient genes'
  save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/"
  # ts_sce <- readRDS(paste0(save_dir,'tradeSeq_v2/', "fitGAM_out.rds"))
  dim(ts_sce)
  obs_lineages <- c(1,2,3,4) ## lineages list that you want to display
  phm <- viz_heatmap(ts_sce, genes_df, save_dir, datatag, plttitle, obs_lineages)
  output_dir <- save_dir
  exp_mtx <- data.table::fread(paste0(output_dir,"mtx_hm.csv.gz")) %>% as.data.frame()
  obs_genes_df <- data.table::fread(paste0(output_dir,"obs_genes_hm.csv.gz")) %>% as.data.frame()
  obs_cells_df <- data.table::fread(paste0(output_dir,"obs_cells_hm.csv.gz")) %>% as.data.frame()
  rownames(exp_mtx) <- exp_mtx$ens_gene_id
  exp_mtx$ens_gene_id <- NULL
  dim(obs_cells_df)
  head(obs_cells_df)
  unique(obs_cells_df$lineage)
  
  dim(exp_mtx)
  dim(obs_genes_df)
  dim(obs_cells_df)
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/cis_trans_lineages/')
  cistrans_anno <- data.table::fread(paste0(input_dir, datatag,'_genes_cis_trans_lineages.csv')) %>% as.data.frame()
  dim(cistrans_anno)
  # head(cistrans_anno)
  rownames(cistrans_anno) <- NULL
  # cistrans_anno$gene_type <- paste0("in ",cistrans_anno$gene_type)
  summary(as.factor(cistrans_anno$gene_type)) ## To Do: check input data here
  
  
  tag <- 'SA1035'
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  chromatin_df <- data.table::fread(paste0(save_fig_dir, "Gene_Module_Chromatin_Status_v2.csv")) %>% 
    as.data.frame() %>%
    dplyr::filter(datatag==tag)
  
  chromatin_df$Module <- paste0("Module",chromatin_df$Module)
  unique(obs_genes_df$gene_type)
  unique(chromatin_df$Module)
  chr_st <- chromatin_df$chromatin_status
  names(chr_st) <- chromatin_df$Module
  obs_genes_df$chromatin_st <- chr_st[obs_genes_df$gene_type]
  summary(as.factor(obs_genes_df$gene_type))
  
  
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/cis_trans_lineages/')
  meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
  meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]
  meta_clone_lg
  meta_clone_lg$lineage_desc <- gsub(' ','',meta_clone_lg$lineage_desc)
  meta_clone_lg <- meta_clone_lg %>%  # less important, excluded from analysis here
    dplyr::filter(lineage!='Lineage5') %>% 
    dplyr::select(-clone)
  # sum(obs_genes_df$ens_gene_id==rownames(exp_mtx))
  unique(obs_genes_df$gene_type)
  unique(obs_cells_df$lineage)
  head(obs_cells_df)
  # exp_mtx[1:2,1:3]
  obs_cells_df$lineage <- substr(obs_cells_df$lineage,1,2)
  obs_cells_df$lineage <- gsub('L','Lineage',obs_cells_df$lineage)
  gap <- viz_genes_exp_lineages_cistrans_anno_hm(as.matrix(exp_mtx), obs_genes_df, obs_cells_df,
                                                 cistrans_anno, meta_clone_lg,
                                                 output_dir, plttitle)
  plttitle <- gsub(' ','_',plttitle)
  gap <- readRDS(paste0(output_dir,plttitle,".rds"))
  gap
  
}

plot_pathways <- function(){
  ## pathway bootstrap analysis
  # source('/home/htran/Projects/farhia_project/drug_resistant_material/scripts/pipeline/utils/pathway_utils.R')
  script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/'
  source(paste0(script_dir, "pipeline/utils/pathway_utils.R"))
  tag <- 'SA1035'
  # save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  # datatag <- 'SA609'
  # pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series.csv')) %>% as.data.frame()
  # pathway_stat <- pathway_stat %>%
  #   dplyr::filter(datatag=='SA609') %>%
  #   dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
  # pathway_stat$gene_type_module <- gsub('Module','M',pathway_stat$gene_type_module)
  # my_font <- "Helvetica"
  # pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
  # pathway_stat$datatag <- pathway_stat$gene_type_module
  # pw_609 <- viz_pathways(pathway_stat, datatag, save_fig_dir)
  # pw_609
  # save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  # pw_609 <- readRDS(paste0(save_fig_dir, datatag, "_trajectory_pathways.rds"))
  
  # unique(obs_genes_df$gene_type)
  # unique(obs_genes_df$gene_type_origin)
  library(dplyr)
  library(ggplot2)
  
  
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv')) %>% as.data.frame()
  pathway_stat <- pathway_stat %>%
    # dplyr::filter(datatag==tag) %>%
    dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
  unique(pathway_stat$datatag)
  
  
  my_font <- "Helvetica"
  pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
  pathway_stat1 <- pathway_stat %>%
    dplyr::filter(grepl(tag,datatag) & p_value <0.05) %>%
    dplyr::select(-signif_genes)
  View(pathway_stat1)
  
    ## Adding new labels
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  ## For SA1035 series, we use the same gene module names, but to make it consistent with other analysis, we still add this table here
  meta_gm_SA1035 <- tibble::tibble(gene_type_module=unique(pathway_stat1$gene_type_module),
                                   gm_manuscript_lb=unique(pathway_stat1$gene_type_module))
  data.table::fwrite(meta_gm_SA1035, paste0(save_fig_dir,'meta_gene_module_labels_',tag,'.csv.gz'))
  meta_gm_SA1035 <- data.table::fread(paste0(save_fig_dir,'meta_gene_module_labels_',tag,'.csv.gz'))
  # meta_gm_SA609$gm_manuscript_lb
  meta_gm_SA1035
  summary(as.factor(pathway_stat1$gene_type_module))
  pathway_stat1 <- pathway_stat1 %>%
    dplyr::left_join(meta_gm_SA1035, by=c('gene_type_module')) %>%
    dplyr::rename(gene_type_origin=gene_type_module, gene_type_module=gm_manuscript_lb)
  summary(as.factor(pathway_stat1$gene_type_module))
  head(pathway_stat1)
  dim(pathway_stat1)
  
  pathway_stat1$gene_type_module <- gsub('Module','M',pathway_stat1$gene_type_module)
  ## TO DO: add new gene module names
  p_pathway <- viz_pathways_barplot(pathway_stat1)
  p_pathway
  results_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/plots/'
  saveRDS(p_pathway, paste0(results_dir,datatag, "p_pathway.rds"))
  p_pathway_SA1035 <- readRDS(paste0(results_dir,datatag, "p_pathway.rds"))
  dev.off()
} 


plot_lineage_prevalence <- function(){
  
  datatag <- 'SA1035'
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
  output_dir <- save_dir
  nfeatures_use <- 3000
  sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce_v2.rds'))
  dim(sce)
  
  metacell <- colData(sce) %>% as.data.frame()
  data.table::fwrite(metacell, paste0(output_dir, datatag,'_meta_cells.csv.gz'))
  dim(metacell)
  
  
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/cis_trans_lineages/')
  meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
  meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]
  
  meta_clone_lg$lineage_desc <- gsub(' ','',meta_clone_lg$lineage_desc)
  meta_clone_lg <- meta_clone_lg %>%  # less important, excluded from analysis here
    # dplyr::filter(lineage!='Lineage5') %>% 
    dplyr::select(-clone)
  meta_clone_lg
  
  res <- plot_trajectory_clones_prevalence_SA1035(metacell, meta_clone_lg, output_dir, save_dir, datatag)
  res$pts
  res$pc
  tag <- 'treatment_clone'
  p_clone_prevalence_1035 <- readRDS(paste0(paste0(output_dir,'figs_v3/'),datatag, "_",tag,"_trajectory_prevalence.rds"))
  tag <- 'treatment_desc'
  p_treatment_prevalence_1035 <- readRDS(paste0(paste0(output_dir,'figs_v3/'),datatag, "_",tag,"_trajectory_prevalence.rds"))
  
}


plot_smooth_expression <- function(){
  script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/trajectory_analysis/'
  source(paste0(script_dir, "tradeseq_utils.R"))
  datatag <- 'SA1035'
  
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/'
  # input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/'
  output_dir <- paste0(dirname(input_dir),'/figs_v5/')
  output_dir <- paste0(dirname(input_dir),'/figs/')
  res_smooth_exp <- viz_gene_modules_density_plot(xplt='pseudotime',yplt='gene_exp',
                                                  colorplt='lineage',datatag=datatag,input_dir, output_dir)
  
  p <- readRDS('/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/figs_v31/ts_slingshot_out_wholedataset_SA1035.rds')
  
  ggsave(paste0(output_dir,"SUPP_Fig11_trajectory_SA1035_smooth_exp.svg"),
         plot = res_smooth_exp$p_smooth_lineage_exp,
         height = 7,
         width = 2.5,
         # useDingbats=F,
         dpi=150)
  ## Reloading data
  p <- readRDS(paste0(output_dir,"density_plot_",datatag, ".rds"))
  plg <- readRDS(paste0(output_dir,"density_plot_",datatag, "_lg.rds"))
  res_smooth_exp <- list(p_smooth_lineage_exp=p, plg=plg)
  
}

plot_example_genes_SA1035 <- function(){
  library(tidyverse)
  library(ComplexHeatmap)
 
  datatag <- 'SA1035'
  figs_dir <- paste0("/home/htran/storage/datasets/drug_resistance/rna_results/",datatag,"_rna/slingshot_trajectory/figs_v3/")
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  total_genes <- data.table::fread(paste0(save_dir, datatag,'_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
  dim(total_genes)
  obs_genes_ls <- c('S100A6','CST3','NFKBIA','CSTB')
  meta_genes <- total_genes %>%
    dplyr::filter(gene_symbol %in% obs_genes_ls)
  obs_genes_ls <- meta_genes$ens_gene_id
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/'
  ts_sce <- readRDS(paste0(input_dir,'tradeseq_v2/', "fitGAM_out.rds"))
  dim(ts_sce)
  output_fn <- paste0(input_dir,'tradeseq_v2/','sigf_gene_exp.csv.gz')
  avg_gene_exp <- get_average_gene_exp_per_lineage(datatag, ts_sce, obs_genes_ls, output_fn, nEstimatedPoints=100, save_data=T)
  
  ## Loading from file
  avg_gene_exp <- data.table::fread(output_fn) %>% as.data.frame()
  avg_gene_exp$lineage_desc <- NULL
  # save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
  dim(avg_gene_exp)
  unique()
  obs_genes_ls <- c('S100A6','CST3','NFKBIA','CSTB')
  meta_genes <- total_genes %>%
    dplyr::filter(gene_symbol %in% obs_genes_ls)
  
  library(tidyverse)
  meta_genes <- meta_genes %>% 
    remove_rownames %>%
    dplyr::select(ens_gene_id, gene_symbol) %>%
    tibble::column_to_rownames('gene_symbol')
  
  # head(meta_genes)
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/cis_trans_lineages/')
  
  meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
  meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]
  meta_clone_lg
  meta_clone_lg$lineage_desc <- gsub(' ','',meta_clone_lg$lineage_desc)
  meta_clone_lg <- meta_clone_lg %>%  # less important, excluded from analysis here
    dplyr::filter(lineage!='Lineage5') %>% 
    dplyr::select(-clone)
  
  unique(avg_gene_exp$lineage)
  avg_gene_exp <- avg_gene_exp %>%
    mutate(lineage=gsub('Lineage ','L',lineage)) %>%
    filter(lineage!='L5')
    # inner_join(meta_clone_lg, by='lineage') %>%
    # select(-lineage) %>%
    # rename(lineage=lineage_desc)
  dim(avg_gene_exp)
  obs_genes_ls <- c('S100A6','CST3','NFKBIA','CSTB')
  viz_given_gene_exp_lineages(obs_genes_ls, meta_genes, avg_gene_exp, figs_dir, datatag)
  
  
  
  plt_ls <- list()
  for(gsymb in obs_genes_ls){
    pg <- readRDS(paste0(figs_dir,gsymb,'.rds'))  
    plt_ls[[gsymb]] <- pg + theme(legend.position = 'none')
  }  
  pg <- readRDS(paste0(figs_dir,gsymb,'.rds'))  
  plg_lg <- cowplot::get_legend(pg)
  plt_ls[['lg']] <- cowplot::ggdraw(plg_lg)
  
  
  p_obs_genes <- cowplot::plot_grid(plotlist = plt_ls, nrow = 1)
  p_obs_genes
  # gb <- grid.grabExpr(draw(p1, annotation_legend_side = "bottom",
  #                          heatmap_legend_side = "bottom"),
  #                     padding = unit(c(1, 1, 2, 2), "mm"))#,merge_legend = TRUE
  # 
  # ptotal <- cowplot::plot_grid(gb, p_obs_genes, nrow = 2, rel_heights = c(2,1))
  # png(paste0(output_dir,"reverse_33genes_summary_",datatag,".png"), height = 2*650, width=2*400,res = 2*72)
  # print(ptotal)
  # dev.off()
  # 
  # png(paste0(output_dir,"reverse_33genes_summary_",datatag,".png"), height = 2*550, width=2*650,res = 2*72)
  # print(ptotal)
  # dev.off()
  
}

viz_SUPP_Fig9 <- function(){
  ## Heatmap
  datatag <- 'SA1035'
  plttitle <- 'Activated Repressed Transient genes'
  output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/"
  plttitle <- gsub(' ','_',plttitle)
  gap <- readRDS(paste0(output_dir,plttitle,".rds"))
  # gap  
  gap_plt <- grid.grabExpr(ComplexHeatmap::draw(gap, padding = unit(c(10, 1, 5, 1), "mm")))
  
  ## Prevalence
  tag <- 'treatment_clone'
  p_clone_prevalence_1035 <- readRDS(paste0(paste0(output_dir,'figs_v3/'),datatag, "_",tag,"_trajectory_prevalence.rds"))
  tag <- 'treatment_desc'
  p_treatment_prevalence_1035 <- readRDS(paste0(paste0(output_dir,'figs_v3/'),datatag, "_",tag,"_trajectory_prevalence.rds"))
  
  plt_clone_prevalence_1035 <- grid::grid.grabExpr(ComplexHeatmap::draw(p_clone_prevalence_1035, annotation_legend_side = "bottom",
                                                                       heatmap_legend_side = "right"),
                                                  padding = unit(c(10, 1, 1, 20), "mm"))#,merge_legend = TRUE
  p_treatment_prevalence_1035 <- grid::grid.grabExpr(ComplexHeatmap::draw(p_treatment_prevalence_1035, annotation_legend_side = "bottom",
                                                                           heatmap_legend_side = "right"),
                                                      padding = unit(c(10, 1, 1, 20), "mm"))#,merge_legend = TRUE
  
  prevalence_total_1035 <- cowplot::plot_grid(NULL, p_treatment_prevalence_1035, plt_clone_prevalence_1035, 
                                             rel_widths = c(0.05,1,1.1), nrow = 1)
  
  ## To Do: add legend here
  
  ## pathways
  results_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/plots/'
  p_pathway_SA1035 <- readRDS(paste0(results_dir,datatag, "p_pathway.rds"))  
  
  ## Smooth expression
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/'
  # input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/'
  output_dir <- paste0(dirname(input_dir),'/figs/')
  p <- readRDS(paste0(output_dir,"density_plot_",datatag, ".rds"))
  plg <- readRDS(paste0(output_dir,"density_plot_",datatag, "_lg.rds"))
  res_smooth_exp <- list(p_smooth_lineage_exp=p, plg=plg)
  
  
  ## UMAP plots
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/'
  results_dir <- paste0(output_dir,'figs_v4/')
  l1 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage1_',datatag,'.rds'))
  # l1 <- l1 + theme(plot.title = element_text(color="black", size=13, hjust = 0.5, family=my_font))
  l2 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage2_',datatag,'.rds'))
  # l2 <- l2 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))
  l3 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage3_',datatag,'.rds'))
  # l3 <- l3 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))
  l4 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage4_',datatag,'.rds'))
  # l4 <- l4 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))
  
  save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/figs_v4/"
  lg <- readRDS(paste0(save_dir, '_plt_legend.rds'))
  plg_umap <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  
  p_umap <- cowplot::plot_grid(l4, l3, l2, l1, nrow=2)
  p_umap <- cowplot::plot_grid(plg_umap, p_umap, NULL, ncol=1, rel_heights = c(0.15, 1, 0.1)) ## legends
  p1 <- cowplot::plot_grid(NULL, p_umap, plt_treatment_prevalence_1035, NULL, plt_clone_prevalence_1035, NULL,
                           rel_widths = c(0.05, 2.1, 1, 0.08,1.3, 0.2), nrow=1)
  # p1
  
  p2 <- cowplot::plot_grid(res_smooth_exp$p_smooth_lineage_exp, p_pathway_SA1035, nrow=1, rel_widths=c(0.7,1))
  p22 <- cowplot::plot_grid(NULL, p2, res_smooth_exp$plg, ncol=1, rel_heights=c(0.15, 2, 0.15))
  p3 <- cowplot::plot_grid(gap_plt, p22, nrow=1, rel_widths = c(1, 1))
  
  
  p4 <- cowplot::plot_grid(p1, p3, NULL, ncol=1, rel_heights = c(1, 1, 0.1)) + 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  dev.off()
  save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/figs/'
  ggsave(paste0(save_dir,"trajectory_results_SA1035.png"),
         plot = p4,
         height=12,
         width=12,
         # useDingbats=F,
         dpi=250)
  ggsave(paste0(save_dir,"trajectory_results_SA1035_v2.svg"),
         plot = p4,
         height=13,
         width=13,
         # useDingbats=F,
         dpi=200)
  
  ggsave(paste0(save_dir,"SUPP_Fig9_SA1035_trajectory_prevalence_part1.svg"),
         plot = prevalence_total_1035,
         height = 6,
         width = 7,
         # useDingbats=F,
         dpi=250)
  
}




