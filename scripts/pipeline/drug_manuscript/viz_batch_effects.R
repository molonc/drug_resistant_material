# res_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/'
# datatag <- 'SA1035'
# p1035_hk <- readRDS(paste0(res_dir,datatag,'/',datatag,'_HK_HK_eval_v2.rds'))
# p1035_raw <- readRDS(paste0(res_dir,datatag,'/',datatag,'_raw_data_features_v2.rds'))
# 
# datatag <- 'SA609'
# p609_hk <- readRDS(paste0(res_dir,datatag,'/',datatag,'_HK_HK_eval_v2.rds'))
# p609_raw <- readRDS(paste0(res_dir,datatag,'/',datatag,'_raw_data_features_v2.rds'))
# 
# datatag <- 'SA535_cisplatin'
# p535_cis_hk <- readRDS(paste0(res_dir,datatag,'/','SA535_CIS_HK_HK_eval_v2.rds'))
# p535_cis_raw <- readRDS(paste0(res_dir,datatag,'/',datatag,'_raw_data_features_v2.rds'))
# 
# datatag <- 'SA535_CX5461'
# p535_cx_hk <- readRDS(paste0(res_dir,datatag,'/','SA535_CX5461_HK_HK_eval_v2.rds'))
# p535_cx_raw <- readRDS(paste0(res_dir,datatag,'/',datatag,'_raw_data_features_v2.rds'))
# 
# p_total <- cowplot::plot_grid(p609_hk, p1035_hk, p535_cis_hk, p535_cx_hk, nrow = 4, align='v',
#                               labels = c('a','b','c','d'))
# 
# png(paste0(res_dir,"diagnostic_plot_3series_batch_effect_eval.png"), height = 2*1200, width=2*1350,res = 2*72)
# print(p_total)
# dev.off()
# 
# # options("scipen"=-100, "digits"=4)
# p_raw <- cowplot::plot_grid(p609_raw, p1035_raw, p535_cis_raw, p535_cx_raw, nrow = 4, 
#                             align='v', labels = c('a','b','c','d'))
# 
# png(paste0(res_dir,"diagnostic_plot_3series_raw_data_features_08_Feb.png"), height = 2*1150, width=2*1100,res = 2*72)
# print(p_raw)
# dev.off()



source('/home/htran/Projects/farhia_project/rnaseq/pipeline/utils/normalize_utils.R')
# meta_ls <- p_raw$meta_data + p_scran$meta_data + p_seurat$meta_data + p_sctransform$meta_data
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA1035'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
meta_SA1035 <- readRDS(paste0(save_dir,datatag,'_HK_eval_metadata_v2.rds'))
dim(meta_SA1035)
p1035_hk <- readRDS(paste0(save_dir,"batch_effect_HK_plots.rds"))
# p1035_hk <- plot_batch_effects(meta_SA1035, xstring="series", ystring="mean_exp", 
#                               plottype="batch_label", plottitle=paste0(datatag),
#                               xlabel=NULL, ylabel=paste0('456 HK genes avg exp'),
#                               lg_pos="bottom", save_dir)

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA609'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
meta_SA609 <- readRDS(paste0(save_dir,datatag,'_HK_eval_metadata_v2.rds'))
dim(meta_SA609)
View(head(meta_SA609))
meta_SA609 <- get_metadata(meta_SA609)
meta_SA609$batch
p609_hk <- readRDS(paste0(save_dir,"batch_effect_HK_plots.rds"))
p609_hk <- plot_batch_effects(meta_SA609, xstring="series", ystring="mean_exp", 
                              plottype="batch_label", plottitle=paste0(datatag),
                              xlabel=NULL, ylabel=paste0('512 HK genes avg exp'),
                              lg_pos="bottom", save_dir)


input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/others_dataset/SA535_cisplatin/'
datatag <- 'SA535'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/SA535_cisplatin/')
meta_SA535_cis <- readRDS(paste0(save_dir,datatag,':Cisplatin_HK_eval_metadata_v2.rds'))
p535_cis_hk <- readRDS(paste0(save_dir,"batch_effect_HK_plots.rds"))
p535_cis_hk <- plot_batch_effects(meta_SA535_cis, xstring="series", ystring="mean_exp", 
                                  plottype="batch_label", plottitle=paste0(datatag,":Cisplatin"),
                                  xlabel=NULL, ylabel=paste0('362 HK genes avg exp'),
                                  lg_pos="bottom", save_dir)



input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/others_dataset/SA535_cisplatin/'
datatag <- 'SA535'
tag <- 'CX5461'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'_',tag,'/')
meta_SA535_cx <- readRDS(paste0(save_dir,datatag,':CX5461_HK_eval_metadata_v2.rds'))
dim(meta_SA535_cx)
# p535_cx_hk <- readRDS(paste0(save_dir,"batch_effect_HK_plots.rds"))
meta_SA535_cx$batch <- ifelse(meta_SA535_cx$batch=='CHIP001','CHIP0047',
                              ifelse(meta_SA535_cx$batch=='CHIP003' | meta_SA535_cx$batch=='CHIP004','CHIP0192',meta_SA535_cx$batch))

p535_cx_hk <- plot_batch_effects(meta_SA535_cx, xstring="series", ystring="mean_exp", 
                                 plottype="batch_label", plottitle=paste0(datatag,":CX5461"),
                                 xlabel=NULL, ylabel=paste0('359 HK genes avg exp'),
                                 lg_pos="bottom", save_dir)




phk <- cowplot::plot_grid(p609_hk, p1035_hk, p535_cis_hk, p535_cx_hk, ncol=4, 
                          labels = c('a','b','c','d'),
                          align = 'h', rel_widths = c(1,0.8,1,1))

save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/'
ggsave(paste0(save_dir,"batch_effect_HK_diagnostic_plots_v3.pdf"),
       plot = phk,
       height = 8.5,
       width = 8,
       units = 'in',
       useDingbats=F
)#dpi = 150


ggsave(filename = paste0(save_dir,"batch_effect_HK_plots_v3.png"),
       plot = phk,
       height = 8,
       width = 8.5,
       # useDingbats=F,
       type = "cairo-png"
)

ggsave(filename = paste0(save_dir,"batch_effect_HK_plots_SA535_cis.png"),
       plot = p535_cis_hk,
       height = 8,
       width = 4,
       # useDingbats=F,
       type = "cairo-png"
)

ggsave(filename = paste0(save_dir,"batch_effect_HK_plots_SA535_cx.png"),
       plot = p535_cx_hk,
       height = 8,
       width = 4,
       # useDingbats=F,
       type = "cairo-png"
)

ggsave(filename = paste0(save_dir,"batch_effect_HK_plots_SA1035.png"),
       plot = p1035_hk,
       height = 8,
       width = 4,
       # useDingbats=F,
       type = "cairo-png"
)

ggsave(filename = paste0(save_dir,"batch_effect_HK_plots_SA609.png"),
       plot = p609_hk,
       height = 8,
       width = 4,
       # useDingbats=F,
       type = "cairo-png"
)
# ggsave(
#   filename = paste0(save_dir, 'LOH_profile_cairo.png'),
#   plot = cnv_plot,
#   height = 3,
#   width = w,
#   # useDingbats=F
#   type = "cairo-png")

# saveRDS(p535_cis_raw, paste0(output_dir,datatag,"_raw_data_features_v2.rds"))
# png(paste0(output_dir,datatag,"_raw_data_features.png"), height = 2*420, width=2*860,res = 2*72)
# print(p535_cis_raw)
# dev.off()



# Get statistic output
get_hk_stat <- function(df, datatag, nbhk){
  df <- get_metadata(df)
  df <- df %>%
    dplyr::group_by(series, tag) %>%
    dplyr::summarise(avg_exp=round(mean(mean_exp),2),
                     med_exp=round(median(mean_exp),2),
                     iqr_exp=round(IQR(mean_exp),2))%>%
    ungroup()
  
  df2 <- df %>%
    dplyr::group_by(tag) %>%
    dplyr::summarise(min_med_exp=min(med_exp),
                     max_med_exp=max(med_exp),
                     med_avg_exp=round(median(avg_exp),2),
                     min_iqr_exp=round(min(iqr_exp),2),
                     max_iqr_exp=round(max(iqr_exp),2))%>%
    ungroup() %>%
    as.data.frame()
  desc <- paste0("Series ",datatag, ': ',nbhk, ' HK genes average expression over all samples: ')
  rownames(df2) <- df2$tag
  for(tag in df2$tag){
    desc <- paste0(desc,tag,': median range: [',
                   df2[tag,'min_med_exp'],'-',df2[tag,'max_med_exp'],']  ',
                   ' interquartile: [',df2[tag,'min_iqr_exp'],'-',df2[tag,'max_iqr_exp'],']    ')
  }
  
  return(list(datatag=datatag,
              desc=desc, 
              df_series=df,
              df_summary=df2))
}  


hstat_609 <- get_hk_stat(meta_SA609, 'SA609', 512)
hstat_609$desc

hstat_1035 <- get_hk_stat(meta_SA1035, 'SA1035', 456)
hstat_1035$desc

hstat_SA535_cis <- get_hk_stat(meta_SA535_cis, 'SA535:Cisplatin', 362)
hstat_SA535_cis$desc

hstat_SA535_cx <- get_hk_stat(meta_SA535_cx, 'SA535:CX5461', 359)
hstat_SA535_cx$desc


# p_hk <- plot_batch_estimation_v2(sce_scMerge, sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'HK', T)

# ref <- res$filtered_gene_attr_scSEG
# ref <- ref %>%
#   dplyr::filter(abs(log_var)<0.2)
# stable_genes <- ref$gene_id
# length(stable_genes)
# p_ref_scSEG <- plot_batch_estimation_v2(sce_scMerge, sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'ref scSEG', T)
# 
# ref <- res$filtered_our_gene_attr_scSEG
# ref <- ref %>%
#   dplyr::filter(abs(log_var)<0.2)
# stable_genes <- ref$gene_id
# length(stable_genes)
# 
# p_custom_scSEG <- plot_batch_estimation_v2(sce_scMerge, sce_raw, sce_scran, sce_seurat, sce_transform, stable_genes, save_dir, datatag, 'custom scSEG', T)
# 
# p_total <- cowplot::plot_grid(p_hk, p_ref_scSEG, p_custom_scSEG, ncol = 3, align='hv', rel_widths = c(1.45,1.2,1.2))
# png(paste0(save_dir,datatag,"_total_eval.png"), height = 2*970, width=2*1200,res = 2*72)
# print(p_total)
# dev.off()


