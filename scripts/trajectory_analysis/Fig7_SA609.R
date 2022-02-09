# SA609
script_dir <- '/home/htran/Projects/farhia_project/rnaseq/trajectory_analysis/'
source(paste0(script_dir, "slingshot_utils.R"))
source(paste0(script_dir, "tradeseq_utils.R"))
library(grid)
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/'
results_dir <- paste0(output_dir,'figs_v3/')
my_font <- "Helvetica"
datatag <- 'SA609'

## Trajectory plot
## See output of res <- plot_all_lingeages(sce, crv_umap_embed, output_dir, datatag)
l1 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage1_SA609.rds'))
l1 <- l1 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))

l2 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage2_SA609.rds'))
l2 <- l2 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))

l3 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage3_SA609.rds'))
l3 <- l3 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))

total <- readRDS(paste0(results_dir,'ts_slingshot_out_wholedataset_SA609.rds'))
total <- total + theme(plot.title = element_text(color="black", size=20, hjust = 0, family=my_font, face="bold"))
ts_color <- readRDS(paste0(results_dir,'treatment_plt.rds'))
# ts_color <- lg2


p1 <- cowplot::plot_grid(l3, l2, l1, ncol=1)
p2 <- cowplot::plot_grid(total, ts_color, ncol=1, rel_heights = c(4,1))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))
ptotal <- cowplot::plot_grid(p2, p1, rel_widths = c(2,1), nrow = 1)
ptotal
dev.off()
png(paste0(results_dir,"slingshot_out_summary_",datatag,"_part11.png"), height = 2*500, width=2*700,res = 2*72)
print(ptotal)
dev.off()


## Prevalence plot
pts <- readRDS(paste0(results_dir,datatag,'_treatment_trajectory_prevalence.rds'))
pc <- readRDS(paste0(results_dir,datatag,'_clone_trajectory_prevalence.rds'))
# gb_heatmap = grid.grabExpr(draw(Heatmap(mat)))
# The value of padding should be a unit vector with length of four. The four values correspond to 
# the space at the bottom, left, top and right sides.
prevalence_ts <- grid.grabExpr(ComplexHeatmap::draw(pts, annotation_legend_side = "bottom",
                                                    heatmap_legend_side = "bottom",
                                                    padding = unit(c(8, 2, 2, 2), "mm")))
prevalence_clone <- grid.grabExpr(ComplexHeatmap::draw(pc, annotation_legend_side = "bottom",
                                                       heatmap_legend_side = "bottom",
                                                       padding = unit(c(8, 2, 2, 2), "mm")))
pvalenc <- cowplot::plot_grid(prevalence_ts, prevalence_clone,rel_widths = c(1,1.1), ncol=2)
ptraj <- cowplot::plot_grid(ptotal, pvalenc,
                            rel_widths = c(3,2.2), ncol=2, labels = c('a','b'))
png(paste0(results_dir,"trajectory_",datatag,"_part1.png"), height = 2*400, width=2*900,res = 2*72)
print(ptraj)
dev.off()


## ## Specific genes in each module plot
## Output of function: viz_given_gene_exp_lineages(obs_genes_ls, meta_genes, ts_sce, figs_dir, datatag)

results_dir2 <- paste0(output_dir,'tradeSeq_3000/SA609_tradeSeq/')
# demo_genes <- c('TOP2A','COX6C','ID1','CCDC26','HIST1H4C')
demo_genes <- c('TOP2A','COX6C','ID1','CCDC26')
# titles <- c('TOP2A: Module 1','COX6C: Module 4','ID1: Module 2','CCDC26: Module 5')
titles <- c('TOP2A: M1','COX6C: M4','ID1: M2','CCDC26: M5')

plt_genes_demo <- list()
for(i in seq_len(length(demo_genes))){
  g <- demo_genes[i]
  exg <- readRDS(paste0(results_dir,g,'.rds'))
  exg <- exg + theme(plot.title = element_text(color="black", size=13, hjust = 0.5, family=my_font))
  exg <- exg + labs(title = titles[i])
  if(i%%2==0){
    exg <- exg + theme(legend.position = 'none')
  }
  plt_genes_demo[[g]] <- exg
}

pgene <- cowplot::plot_grid(plotlist = plt_genes_demo, ncol=2) #, labels = c('a','b','c','d')
pgene

## Heatmap average expression plot
## Extracting heatmap of average genes expression x lineages 
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/'
# genes_df <- data.table::fread(paste0(save_dir, 'SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
# dim(genes_df)
# colnames(genes_df)
# genes_df <- genes_df %>%
#   dplyr::select(-gene_type) %>%
#   dplyr::rename(gene_type=gene_type_module)
# summary(as.factor(genes_df$gene_type))
datatag <- 'SA609'
plttitle <- 'Activated Repressed Transient genes'
# save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
# ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
# save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/"
# phm <- viz_heatmap(ts_sce, genes_df, save_dir, datatag, plttitle)
# output_dir <- save_dir
# exp_mtx <- data.table::fread(paste0(output_dir,"mtx_hm.csv.gz")) %>% as.data.frame()
# obs_genes_df <- data.table::fread(paste0(output_dir,"obs_genes_hm.csv.gz")) %>% as.data.frame()
# obs_cells_df <- data.table::fread(paste0(output_dir,"obs_cells_hm.csv.gz")) %>% as.data.frame()
# rownames(exp_mtx) <- exp_mtx$ens_gene_id
# exp_mtx$ens_gene_id <- NULL
# input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/cis_trans_lineages/')
# cistrans_anno <- data.table::fread(paste0(input_dir, datatag,'_genes_cis_trans_lineages.csv')) %>% as.data.frame()
# dim(cistrans_anno)
# rownames(cistrans_anno) <- NULL
# cistrans_anno$gene_type <- paste0("in ",cistrans_anno$gene_type)
# meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
# meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]
# 
# # sum(obs_genes_df$ens_gene_id==rownames(exp_mtx))
# gap <- viz_genes_exp_lineages_cistrans_anno_hm(as.matrix(exp_mtx), obs_genes_df, obs_cells_df, 
#                                                cistrans_anno,meta_clone_lg,
#                                                output_dir, plttitle)

# p <- viz_genes_exp_lineages_hm(as.matrix(exp_mtx), obs_genes_df, obs_cells_df, output_dir, plttitle)
# p

library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(-4,0,4), c("blue", "white", "red"))
lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "Avg Exp", at = c(-4,-2,0, 2, 4), direction = "horizontal")
hm_plg_exp <- grid.grabExpr(ComplexHeatmap::draw(lgd))
# phm_total <- cowplot::plot_grid(hm_plg, NULL, ncol=2, rel_widths = c(1,0.1),labels = c('c','d'))
# phm_total
output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/"
gap <- readRDS(paste0(output_dir,'Activated_Repressed_Transient_genes.rds'))
gap

dev.off()
gap1 <- grid.grabExpr(ComplexHeatmap::draw(gap, annotation_legend_side = "bottom",
                                           heatmap_legend_side = "bottom",
                                           padding = unit(c(1, 1, 1, 1), "mm")))

hm_plg_cistrans <- readRDS(paste0(output_dir,"gene_type_hm_legend.rds"))

## in case necessary for yaxis label, otherwise, just add a text
# yaxis_label = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_text(c('Pseudotime'), 
#                                                                                 location = 0.2, rot = 0,
#                                                                                 just = "center",which='column'))
# hm_yaxis_title <- grid.grabExpr(ComplexHeatmap::draw(yaxis_label))
# cowplot::plot_grid(hm_yaxis_title, ncol=1)

# gdh1 <- grid.grabExpr(ComplexHeatmap::draw(gdh, annotation_legend_side = "bottom",
#                                            heatmap_legend_side = "bottom",
#                                            padding = unit(c(3, 3, 3, 3), "mm")))
# phm <- cowplot::plot_grid(gap1, gdh1, rel_widths = c(1,1), ncol=2)

# pgene <- cowplot::plot_grid(exg2, exg2, exg2, exg2, exg1, ncol=5, rel_widths = c(1,1,1,1,1.5))
# pgene <- cowplot::plot_grid(plotlist = plt_genes_demo, ncol=2) #, labels = c('a','b','c','d')
# pgene_lg <- cowplot::plot_grid(exg_lg, NULL, ncol=2, rel_widths = c(1.4,0.6)) #,labels = c('f','')




## pathway bootstrap analysis
# save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
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
save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
pw_609 <- readRDS(paste0(save_fig_dir, datatag, "_trajectory_pathways.rds"))

# lg1 <- cowplot::plot_grid(plg, hm_plg, rel_heights = c(3,1), nrow=2)
lg1 <- cowplot::plot_grid(hm_plg_exp, hm_plg_cistrans, rel_heights = c(1,1), nrow=2)+
  theme(plot.background = element_rect(fill = "white", colour = "white"))
pgene_lg <- cowplot::plot_grid(lg1, pw_609, rel_widths = c(0.9,3), nrow=1)
pgenes_ls <- cowplot::plot_grid(pgene, pgene_lg, ncol=1, rel_heights = c(2,1))

hm_total <- cowplot::plot_grid(NULL, NULL, gap1, nrow=3, rel_heights = c(0.2, 0.5, 10),
                               hjust = -0.1,label_fontface="plain",
                               labels=c(' ','Heatmap x: Gene Modules/ y: Pseudotime','')) + 
  theme(plot.background = element_rect(fill = "white", colour = "white"))
phm_total <- cowplot::plot_grid(hm_total, pgenes_ls, ncol=2, rel_widths = c(1,1),labels = c('c','d'))

p609 <- cowplot::plot_grid(ptraj, phm_total, nrow=2, rel_heights = c(1,1.5))
png(paste0(results_dir,"Fig7_trajectory_",datatag,".png"), height = 2*1200, width=2*950,res = 2*72)
print(p609)
dev.off()

ggsave(paste0(results_dir,"Fig7_trajectory_SA609.png"),
       plot = p609,
       height = 13.5,
       width = 12,
       # useDingbats=F,
       dpi=150)


ggsave(paste0(results_dir,"Fig7_trajectory_SA609.pdf"),
       plot = p609,
       height = 13.5,
       width = 12,
       useDingbats=F,
       dpi=150)


