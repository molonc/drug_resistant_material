# df <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/tradeseq/patternRes.csv')
# dim(df)
# colnames(df)
# df <- df %>%
#   dplyr::filter(pvalue<0.05 & gene_symb=='MYC')
# 'MYC' %in% df$gene_symb

source('/home/htran/Projects/farhia_project/rnaseq/pipeline/utils/plot_utils.R')

clones <- c('R','H')
# clones <- c('A','H')
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
base_name <- 'SA609'
df_cnv_fn <- paste0(save_dir,base_name, '_cnv_mat.csv.gz')

# meta_genes <- data.frame(ensembl_gene_id=rownames(norm_sce), stringsAsFactors=F)
# data.table::fwrite(meta_genes, paste0(save_dir,base_name, '_meta_filtered_genes.csv'))
meta_genes <- data.table::fread(paste0(save_dir,base_name, '_meta_filtered_genes.csv')) %>% as.data.frame()
dim(meta_genes)

pcnv_609 <- plot_CNV(df_cnv_fn, clones, meta_genes)
pcnv_609$cnv_plot
pcnv_609$plg
saveRDS(pcnv_609, paste0(save_dir,base_name, '_pcnv.rds'))
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/signif_genes/'
save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/trackplots/'
# dir.create(save_fig_dir)
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# plottitle <- 'SA609 cisplatin: logFC X7-Rx clone A vs X7-UnRx clone H'
plottitle <- NULL
# desc <- 'log2FC X7-Rx clone A vs X7-UnRx clone H'
# desc <- 'SA609: Rx X7:cloneA vs. UnRx X7:cloneH'
desc <- 'Pt4: Rx X7:cloneA vs. UnRx X7:cloneH'
datatag <- 'SA609'
tag <- 'Pt4'

# deg_fn <- paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTTT_R_UUUUU_H/signif_genes.csv')
deg_fn <- paste0(input_dir,'scrande_SA609_4/signif_genes.csv')
# additional_genes_fn <- paste0(base_dir, 'rnaseq_v6/trackplots_v3/SA609_R_vs_H.csv')
additional_genes_fn <- paste0(save_fig_dir, 'SA609_R_vs_H_annotated_genes.csv')
# subtitle <- paste0(datatag,': in cis DE genes, Rx X7:cloneA vs. UnRx X7:cloneH')
# subtitle <- paste0(tag,': Rx X7:cloneA vs. UnRx X7:cloneH, in cis DE genes')
subtitle <- ' '
track_incis <- plot_gex_cnv_v2(deg_fn=deg_fn,
                               df_cnv = pcnv_609$df_cnv,
                               clones=NULL,
                               save_dir=save_fig_dir,
                               sample_name=NULL,
                               clone_str=desc,
                               additional_genes_fn = additional_genes_fn,
                               n_genes_to_annotate = 10, 
                               plttitle=subtitle,
                               gene_type='in-cis')

# subtitle <- paste0(datatag,': in trans DE genes')
# subtitle <- paste0(tag,': Rx X7:cloneA vs. UnRx X7:cloneH, in trans DE genes')
subtitle <- ' '
track_intrans <- plot_gex_cnv_v2(deg_fn=deg_fn,
                                 df_cnv = pcnv_609$df_cnv,
                                 clones=NULL,
                                 save_dir=save_fig_dir,
                                 sample_name=NULL,
                                 clone_str=desc,
                                 additional_genes_fn = additional_genes_fn,
                                 n_genes_to_annotate = 10, 
                                 plttitle=subtitle,
                                 gene_type='in-trans')

prop_cistrans_plt <- plot_cis_trans_prop(deg_fn, lg_pos="none")

# + theme(axis.title.x = element_blank(),
#         strip.text.x = element_blank(),
#         legend.position = 'none')
main_plot <- cowplot::plot_grid(
  track_intrans$track_plot,
  track_incis$track_plot,
  prop_cistrans_plt$p ,
  pcnv_609$cnv_plot,
  ncol = 1,
  rel_heights = c(1.6,1.6,0.8,1.3),
  align = 'v'#,
  # labels = c('a')
)


main_plot_SA609 <- main_plot
dir.create(save_fig_dir)
saveRDS(main_plot_SA609, paste0(save_fig_dir,'main_plot_SA609.rds'))
# track_incis + theme(axis.title.x = element_blank(),
#                     strip.text.x = element_blank()),
save_fn <- paste0(datatag, '_',gsub(' ','_',desc))
save_fn <- gsub(':','_',save_fn)
png(paste0(save_fig_dir,save_fn,"_track_plot.png"), height = 2*490, width=2*1000, res = 2*72)
print(main_plot)
dev.off()
ggsave(paste0(save_dir,"Fig23_part1_logFC_legend.png"),
       plot = track_intrans$plg,
       height = 0.8,
       width = 3,
       # useDingbats=F,
       dpi=200)

ggsave(paste0(save_dir,"Fig23_part1_chr_legend.png"),
       plot = pcnv_609$plg,
       height = 0.8,
       width = 6,
       # useDingbats=F,
       dpi=200)
ggsave(paste0(save_dir,"Fig23_part1_cistrans_chr_legend.png"),
       plot = prop_cistrans_plt$plg,
       height = 0.8,
       width = 6,
       # useDingbats=F,
       dpi=200)

ggsave(paste0(save_dir,"Fig23_part1_track.png"),
       plot = main_plot,
       height = 5,
       width = 7,
       # useDingbats=F,
       dpi=200)



clones <- c('A','G')
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
base_name <- 'SA535'
save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/trackplots/'

df_cnv_fn <- paste0(save_dir,base_name, '_cnv_mat.csv')

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
pcnv_535 <- plot_CNV(df_cnv_fn, clones, meta_genes)
saveRDS(pcnv_535,paste0(save_fig_dir,base_name, '_pcnv.rds'))
pcnv_535 <- readRDS(paste0(save_fig_dir,base_name, '_pcnv.rds'))
# DEG plot
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/trackplots/'
# dir.create(save_fig_dir)
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# desc <- 'log2FC X10-Rx T vs X9-UnRx J'
desc <- 'log2FC X10-Rx A vs. X9-UnRx G'
datatag <- 'SA535'

save_dir <- save_fig_dir

# deg_fn <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UUTTTT_T_UUUUUU_J/signif_genes.csv')
deg_fn <- paste0(input_dir,datatag,'/SA535_2_SA535_UUTTTT_A_UUUUUU_G/signif_genes_FDR_0.01.csv')
# additional_genes_fn <- paste0(base_dir, 'rnaseq_v6/trackplots_v3/SA535_ST_vs_J.csv')
datatag <- 'Pt5'
additional_genes_fn <- paste0(save_fig_dir, 'SA535_ST_vs_J_annotated_genes.csv')
# subtitle <- paste0(datatag,': in cis DE genes, Rx X10:cloneT vs UnRx X9:cloneJ')
subtitle <- paste0(datatag,': in cis DE genes, Rx X10:cloneA vs. UnRx X9:cloneG')
track_incis_535cis <- plot_gex_cnv_v2(deg_fn=deg_fn,
                                      df_cnv=pcnv_535$df_cnv,
                                      clones=NULL,
                                      save_dir=save_fig_dir,
                                      sample_name=paste0(datatag,': in-cis DE genes'),
                                      clone_str=desc,
                                      additional_genes_fn = additional_genes_fn,
                                      n_genes_to_annotate = 10, 
                                      plttitle=subtitle,
                                      gene_type='in-cis')

subtitle <- paste0(datatag,': in trans DE genes, Rx X10:cloneA vs. UnRx X9:cloneG')
track_intrans_535cis <- plot_gex_cnv_v2(deg_fn=deg_fn,
                                        df_cnv=pcnv_535$df_cnv,
                                        clones=NULL,
                                        save_dir=save_fig_dir,
                                        sample_name=subtitle,
                                        clone_str=desc,
                                        additional_genes_fn = additional_genes_fn,
                                        n_genes_to_annotate = 10, 
                                        plttitle=subtitle,
                                        gene_type='in-trans')

main_plot_535cis <- cowplot::plot_grid(
  track_intrans_535cis + theme(axis.title.x = element_blank(),
                               strip.text.x = element_blank()),
  track_incis_535cis + theme(axis.title.x = element_blank(),
                             strip.text.x = element_blank(),
                             legend.position = 'none'),
  pcnv_535$cnv_plot,
  ncol = 1,
  rel_heights = c(1.8,1.6,1.1),
  align = 'v'
)

# track_incis + theme(axis.title.x = element_blank(),
#                     strip.text.x = element_blank()),
save_fn <- paste0(datatag, '_',gsub(' ','_',desc))
png(paste0(save_dir,save_fn,"_track_plot.png"), height = 2*490, width=2*1000, res = 2*72)
print(main_plot_535cis)
dev.off()
saveRDS(main_plot_535cis, paste0(save_fig_dir,'main_plot_',base_name,'.rds'))




# length(intersect(df_track$ensembl_gene_id, df_cnv$ensembl_gene_id))
clones <- c('H','E')
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
base_name <- 'SA1035'
df_cnv_fn <- paste0(save_dir,base_name, '_cnv_mat.csv')
# sce <- readRDS('/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA1035-v6/total_sce_v3.rds')
# dim(sce)
# meta_genes <- data.frame(ensembl_gene_id=rowData(sce)$ID, gene_symbol=rowData(sce)$Symbol, stringsAsFactors=F)
# data.table::fwrite(meta_genes, paste0(save_dir,base_name, '_meta_filtered_genes_10percent.csv'))
meta_genes <- data.table::fread(paste0(save_dir,base_name, '_meta_filtered_genes_10percent.csv')) %>% as.data.frame()
dim(meta_genes)
head(meta_genes)
pcnv_1035 <- plot_CNV(df_cnv_fn, clones, meta_genes)
pcnv_1035$cnv_plot
dim(pcnv_1035$df_cnv)
saveRDS(pcnv_1035, paste0(save_fig_dir,base_name, '_pcnv.rds'))
pcnv_1035 <- readRDS(paste0(save_fig_dir,base_name, '_pcnv.rds'))

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/trackplots/'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
plottitle <- NULL
desc <- 'log2FC X8-Rx clone H vs. X8-UnRx clone E'
datatag <- 'SA1035'
save_dir <- save_fig_dir

# deg_fn <- paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/SA1035_UTTTT_H_UUUUU_E/signif_genes.csv')
deg_fn <- paste0(input_dir,datatag,'/SA1035_2_SA1035_UTTTT_H_UUUUU_E/signif_genes_FDR_0.01.csv')
additional_genes_fn <- paste0(save_fig_dir, 'SA1035_H_vs_E_annotated_genes.csv')


# subtitle <- paste0(datatag,' cisplatin: in-cis DE genes')
datatag <- 'Pt6'
subtitle <- paste0(datatag,': in cis DE genes, Rx X8:cloneH vs. UnRx X8:cloneE')
track_incis_1035 <- plot_gex_cnv_v2( deg_fn=deg_fn,
                                     df_cnv = pcnv_1035$df_cnv,
                                     clones=NULL,
                                     save_dir=save_fig_dir,
                                     sample_name=subtitle,
                                     clone_str=desc,
                                     additional_genes_fn = additional_genes_fn,
                                     n_genes_to_annotate = 10, 
                                     plttitle=subtitle,
                                     gene_type='in-cis')

# subtitle <- paste0(datatag,' cisplatin: in-trans DE genes')
subtitle <- paste0(datatag,': in trans DE genes, Rx X8:cloneH vs. UnRx X8:cloneE')
track_intrans_1035 <- plot_gex_cnv_v2(deg_fn=deg_fn,
                                      df_cnv = pcnv_1035$df_cnv,
                                      clones=NULL,
                                      save_dir=save_fig_dir,
                                      sample_name=subtitle,
                                      clone_str=desc,
                                      additional_genes_fn = additional_genes_fn,
                                      n_genes_to_annotate = 10, 
                                      plttitle=subtitle,
                                      gene_type='in-trans')

main_plot_1035 <- cowplot::plot_grid(
  track_intrans_1035 + theme(axis.title.x = element_blank(),
                             strip.text.x = element_blank()),
  track_incis_1035 + theme(axis.title.x = element_blank(),
                           strip.text.x = element_blank(),
                           legend.position = 'none'),
  pcnv_1035$cnv_plot,
  ncol = 1,
  rel_heights = c(1.8,1.6,1.1),
  align = 'v'
)
save_fn <- paste0(datatag, '_',gsub(' ','_',desc))
png(paste0(save_dir,save_fn,"_track_plot.png"), height = 2*490, width=2*1000, res = 2*72)
print(main_plot_1035)
dev.off()

saveRDS(main_plot_1035, paste0(save_fig_dir,'main_plot_',base_name,'.rds'))


# track_plot <- cowplot::plot_grid(
#   main_plot_SA609,
#   main_plot_535cis,
#   # main_plot_1035,
#   ncol = 1,
#   align = 'v'
# )
# 
# ggsave(paste0(save_fig_dir,"SUPP_fig21_track_plots.pdf"),
#        plot = track_plot,
#        height = 10.5,
#        width = 9.5,
#        useDingbats=F,
#        dpi = 150
#        )
# 
# ggsave(paste0(save_fig_dir,"SUPP_fig21_track_plots_2plots.png"),
#        plot = track_plot,
#        height = 9.5,
#        width = 11,
#        # useDingbats=F,
#        type = "cairo-png",
#        dpi = 200
# )











clones <- c('T','J')
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/clonealign/whole_data/combined_clones/'
base_name <- 'SA535_total'
df_cnv_fn <- paste0(save_dir,base_name, '_cnv_mat_pivot.csv')
df_cnv <- read.csv(df_cnv_fn, check.names = F, stringsAsFactors = F)
pcnv_535_cis <- plot_CNV(df_cnv_fn, clones)

# pcnv_535_cis$cnv_plot

clones <- c('U','J')
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/clonealign/whole_data/combined_clones/'
base_name <- 'SA535_total'
df_cnv_fn <- paste0(save_dir,base_name, '_cnv_mat_pivot.csv')
# df_cnv_fn <- paste0(save_dir,base_name, '_cnv_mat.csv')
# df_cnv_1 <- read.csv(df_cnv_fn, check.names = F, stringsAsFactors = F)
# df_cnv_1 <- df_cnv_1[,clones]
# View(df_cnv_1)
pcnv_535_cx <- plot_CNV(df_cnv_fn, clones)
# pcnv_535_cx$cnv_plot
# 
# 

# cnv_mat <- read.csv(paste0(save_dir,base_name,'_cnv_mat.csv'), check.names = F, stringsAsFactors = F)
# var_genes <- apply(cnv_mat[,colnames(cnv_mat)!='ensembl_gene_id'], 1, var)
# cnv_mat <- cnv_mat[var_genes!=0,]
# cnv_mat <- cnv_mat %>%
#   pivot_longer(!ensembl_gene_id, names_to = "clone", values_to = "cnv")
# 
# print(dim(cnv_mat))
# write.csv(cnv_mat, paste0(save_dir,base_name,'_cnv_mat_pivot_filtered.csv'), row.names = F, quote=F)


clones <- c('R','H')
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/added_segments/clonealign/whole_data/'
base_name <- 'SA609'
df_cnv_fn <- paste0(save_dir, base_name, '_cnv_mat_pivot.csv')
pcnv_609_cis <- plot_CNV(df_cnv_fn, clones)
pcnv_609_cis$cnv_plot   

# pcnv_609_cis
# 
clones <- c('H','E')
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/clonealign/whole_data/'
base_name <- 'SA1035'
df_cnv_fn <- paste0(save_dir,base_name, '_cnv_mat_pivot.csv')
pcnv_1035_cis <- plot_CNV(df_cnv_fn, clones)
# pcnv_1035_cis
# pcnv_1035_cis$cnv_plot
# save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/trackplots/'
# dir.create(save_fig_dir, showWarnings = F)
# saveRDS(list(pcnv_609_cis, pcnv_1035_cis, pcnv_535_cis, pcnv_535_cx), file=paste0(save_fig_dir, 'cnv_ls.rds'))




# input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/rnaseq/differential_expression/results/'
# # logFC_ls <- c('SA609-v6/comps/pathway9_SA609_UTTTT_R_UUUUU_H','SA1035-v6/comps/ressens8_SA1035_UTTTT_H_UUUUU_E',
# #               'SA535-v6/comps/pathway4Tvs4U_SA535_UUTTTT_S_T_UUUUUU_J_Q','SA535-v6/comps/ressens18_SA535_UXXXX_U_UUUUU_J')
# logFC_ls <- c('SA609-v6/comps/pathway9_SA609_UTTTT_R_UUUUU_H','SA1035-v6/comps/ressens8_SA1035_UTTTT_H_UUUUU_E',
#               'SA535-v6/comps/ressens11_SA535_UUTTTT_T_UUUUUU_J','SA535-v6/comps/ressens18_SA535_UXXXX_U_UUUUU_J')
# 
# 
# tag_ls <- c('UTTTT-A vs UUUUU-H', 'UTTTT-H vs UUUUU-E', 'UUTTTT-T vs UUUUUU-J', 'UXXXX-U vs UUUUU-J')
# datatag_ls <- c('SA609 Cisplatin','SA1035 Cisplatin','SA535 Cisplatin','SA535 CX5461')
# additional_genes_ls <- c('SA609_R_vs_H.csv','SA1035_H_vs_E.csv','SA535_ST_vs_J.csv','SA535_U_vs_J.csv')
# df_cnv_ls <- list(as.data.frame(pcnv_609_cis$df_cnv), as.data.frame(pcnv_1035_cis$df_cnv),
#                   as.data.frame(pcnv_535_cis$df_cnv), as.data.frame(pcnv_535_cx$df_cnv))
# 
# df_cnv_plot <- list(pcnv_609_cis$cnv_plot, pcnv_1035_cis$cnv_plot,
#                     pcnv_535_cis$cnv_plot, pcnv_535_cx$cnv_plot)
# 
# 
# n_genes_to_annotate = 10

# for(i in seq(length(logFC_ls))){
#   # deg_fn <- paste0(input_dir,logFC_ls[i],'_logfc_results.csv')
#   # clone_str <- tag_ls[i]
#   # sample_name <- datatag_ls[i]
#   # additional_genes_fn <- paste0(save_dir, 'trackplot/',additional_genes_ls[i])
#   plot_gex_cnv(paste0(input_dir,logFC_ls[i],'_logfc_results.csv'),
#                df_cnv_ls[[i]],
#                df_cnv_plot[[i]],
#                NULL,
#                save_fig_dir,
#                datatag_ls[i],
#                tag_ls[i],
#                paste0(save_fig_dir, additional_genes_ls[i]),
#                n_genes_to_annotate)
# 
# }







# p_vol <- cowplot::plot_grid(vol_609, vol_1035, ncol=2,  align = 'h', rel_widths = c(1,1.25))
# 
# lbs <- c('a','b','c')

# p_total <- cowplot::plot_grid(p_vol, main_plot, main_plot_1035, ncol=1,  rel_heights = c(1.4,2,2), labels = lbs) # align = 'v',

p_total <- cowplot::plot_grid(main_plot, main_plot_1035, ncol=1,  rel_heights = c(2,2), labels = c('a','b'), align = 'v')
png(paste0(save_fig_dir, "SA609_SA1035_track.png"), height =2*1000, width= 2*1000,res = 2*72)
print(p_total)
dev.off()

pdf(paste0(save_fig_dir, "SA609_SA1035_track.pdf"), height =12, width= 10)
print(p_total)
dev.off()

ggsave(paste0(save_fig_dir,"SA609_SA1035_track.pdf"),
       plot = p_total,
       height = 7,
       width = 6.5,
       useDingbats=F,
       dpi = 50)

saveRDS(p_total, paste0(save_fig_dir,"SA609_SA1035_track.rds"))





# desc <- 'log2FC X8-Rx U vs X9-UnRx J'
# datatag <- 'SA535:CX5461'
# save_dir <- save_fig_dir
# 
# deg_fn <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UXXXX_U_UUUUU_J/signif_genes.csv')
# additional_genes_fn <- paste0(base_dir, 'rnaseq_v6/trackplots_v3/SA535_U_vs_J.csv')
# 
# subtitle <- paste0(datatag,': in cis DE genes, Rx X8:cloneU vs UnRx X9:cloneJ')
# track_incis_535cx <- plot_gex_cnv_v2(deg_fn=deg_fn,
#                                      df_cnv=pcnv_535_cx$df_cnv,
#                                       clones=NULL,
#                                       save_dir=save_fig_dir,
#                                       sample_name=paste0(datatag,': in-cis DE genes'),
#                                       clone_str=desc,
#                                       additional_genes_fn = additional_genes_fn,
#                                       n_genes_to_annotate = 10, 
#                                       plttitle=subtitle,
#                                       gene_type='in-cis')
# subtitle <- paste0(datatag,': in trans DE genes, Rx X8:cloneU vs UnRx X9:cloneJ')
# track_intrans_535cx <- plot_gex_cnv_v2(deg_fn=deg_fn,
#                                        df_cnv=pcnv_535_cx$df_cnv,
#                                         clones=NULL,
#                                         save_dir=save_fig_dir,
#                                         sample_name=paste0(datatag,': in-trans DE genes'),
#                                         clone_str=desc,
#                                         additional_genes_fn = additional_genes_fn,
#                                         n_genes_to_annotate = 10, 
#                                         plttitle=plottitle,
#                                         gene_type='in-trans')
# 
# main_plot_535cx <- cowplot::plot_grid(
#   track_intrans_535cx + theme(axis.title.x = element_blank(),
#                                strip.text.x = element_blank()),
#   track_incis_535cx + theme(axis.title.x = element_blank(),
#                              strip.text.x = element_blank()),
#   pcnv_535_cx$cnv_plot,
#   ncol = 1,
#   rel_heights = c(1.8,1.55,0.9),
#   align = 'v'
# )
# 
# # track_incis + theme(axis.title.x = element_blank(),
# #                     strip.text.x = element_blank()),
# save_fn <- paste0(datatag, '_',gsub(' ','_',desc))
# png(paste0(save_dir,save_fn,"_track_plot.png"), height = 2*550, width=2*1100, res = 2*72)
# print(main_plot_535cx)
# dev.off()
# # write.csv(df_track_obs, file=paste0(save_dir,save_fn,"_chr_13_22.csv"), quote=F, row.names = F)
# saveRDS(main_plot_535cx, paste0(save_fig_dir,save_fn,"_track_plot_total.rds"))
# 
# 
# p_vol535 <- cowplot::plot_grid(vol_535_cis, vol_535_cx, ncol=2,  align = 'h', rel_widths = c(1,1.25))
# 
# lbs <- c('a','b')
# 
# p_total_535 <- cowplot::plot_grid(p_vol535, main_plot_535cis, main_plot_535cx, ncol=1,  rel_heights = c(1.5,2,2), labels = lbs) # align = 'v',
# png(paste0(save_fig_dir, "SA535_cis_cx_track.png"), height =2*1500, width= 2*1200,res = 2*72)
# print(p_total_535)
# dev.off()
# 
# 
# p_total_535 <- cowplot::plot_grid(main_plot_535cis, main_plot_535cx, ncol=1,  rel_heights = c(2,2), labels = c('a','b'), align = 'v')
# png(paste0(save_fig_dir, "SA535_cis_cx_track.png"), height =2*1000, width= 2*1000,res = 2*72)
# print(p_total_535)
# dev.off()
# 
# pdf(paste0(save_fig_dir, "SA535_cis_cx_track.pdf"), height =12, width= 10)
# print(p_total_535)
# dev.off()
# 
# ggsave(paste0(save_fig_dir,"SA609_SA1035_track.pdf"),
#        plot = p_total,
#        height = 7,
#        width = 6.5,
#        useDingbats=F,
#        dpi = 50)
# 
# saveRDS(p_total_535, paste0(save_fig_dir,"SA535_cis_cx_track.rds"))
# 
