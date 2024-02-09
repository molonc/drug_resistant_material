t <- umap_df %>%
  dplyr::filter(sample_id=='SA1035X5XB03015')
dim(t)
ass_10 <- data.table::fread(paste0(output_dir, 'clonealign_result_fit_SA1035X5XB03015_min_counts_per_cell_10.csv'))
ass_50 <- data.table::fread(paste0(output_dir, 'clonealign_result_fit_SA1035X5XB03015_min_counts_per_cell_50.csv'))
dim(ass_10)
dim(ass_50)
colnames(ass_10)
sum(ass_10$cell_id %in% t$cell_id)
t10 <- t %>%
  dplyr::select(-clone) %>%
  dplyr::left_join(ass_10, by='cell_id') %>%
  dplyr::mutate(
    clone=case_when(
      !is.na(clone) ~ clone,
      TRUE ~ 'None'
    )
  )
t50 <- t %>%
  dplyr::select(-clone) %>%
  dplyr::left_join(ass_50, by='cell_id') %>%
  dplyr::mutate(
    clone=case_when(
      !is.na(clone) ~ clone,
      TRUE ~ 'None'
    )
  )
base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
source(paste0(base_dir, 'scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))

input_dir <- paste0(base_dir, 'materials/umap_figs/')
output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
# output_dir <- paste0(base_dir, 'materials/umap_figs/testing/')
# dir.create(output_dir)
datatag <- 'SA1035'
basename <- datatag
umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap_filtered_outliers.csv.gz')) %>% as.data.frame()
dim(umap_df)
unique(umap_df$clone)
summary(as.factor(umap_df$treatmentSt))
sids <- sapply(strsplit(umap_df$cell_id,'_'), function(x){
  return(x[1])
})  
umap_df$sample_id <- as.character(sids)

umap_df <- umap_df %>%
  dplyr::mutate(clone=get_unique_clone_id(clone))
color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign

summary(as.factor(t10$clone))
summary(as.factor(t50$clone))
dim(t10)
obs_treatment <- 'Rx'
obs_passage <- 'X5'
plottitle= paste0(obs_treatment,': ',obs_passage)
plottitle <- ''
res10 <- viz_umap_obs_clones_testing(umap_df, t10, cols_use, datatag, 
                                     output_dir, obs_treatment, obs_passage, 
                                     plottitle, T)

res50 <- viz_umap_obs_clones_testing(umap_df, t50, cols_use, datatag, 
                                     output_dir, obs_treatment, obs_passage, 
                                     plottitle, T)



p_total <- cowplot::plot_grid(res10$p, res50$p, labels = c('Thres 10 - 176Un - 8None', 
                                                           'Thres 50 - 116Un - 125None'))
png(paste0(output_dir,"SA1035X5XB03015_clonealign_thresholds.png"), height = 2*450, width=2*850,res = 2*72)
print(p_total)
dev.off()
