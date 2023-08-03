
script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/'
source(paste0(script_dir,'drug_manuscript/cnv_viz_utils.R'))

clones <- c('R','H')
# clones <- c('A','H')
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
base_name <- 'SA609'
df_cnv_fn <- paste0(save_dir,base_name, '_cnv_mat.csv')


input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609/'
# input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/dlp_cnv/Fig6_Mirela_cnv/'
df_cnv_fn <- paste0(input_dir, 'total_merged_filtered_states.csv.gz') #



base_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/cell_clones/'
save_dir <- paste0(input_dir, 'fig_medianCNV/')


datatag <- 'SA501'
copynumber_fn <- paste0(base_dir,'SA501_v2/total_merged_filtered_states.csv')
cellclone_fn <- paste0(input_dir, datatag, '_cell_clones.csv.gz')
library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv.gz') 
df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir,
                                 cellclone_fn, library_grouping_fn) 
dim(df_cnv)
res_501 <- plot_CNV_profile(df_cnv, clones= levels(df_cnv$clone),plttitle='Pt1')
res_501$cnv_plot
res_501$plg
dev.off()




datatag <- 'SA530'
copynumber_fn <- paste0(base_dir,'SA530_v2/total_merged_filtered_states.csv')
cellclone_fn <- paste0(input_dir, datatag, '_cell_clones.csv')
df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir,
                                 cellclone_fn, library_grouping_fn=NULL) 
dim(df_cnv)
res_530 <- plot_CNV_profile(df_cnv, clones= levels(df_cnv$clone),plttitle='Pt2')
res_530$cnv_plot
dev.off()


datatag <- 'SA604'
copynumber_fn <- paste0(base_dir,'SA604_v2/total_merged_filtered_states.csv')
cellclone_fn <- paste0(input_dir, datatag, '_cell_clones.csv')
df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir,
                                 cellclone_fn, library_grouping_fn=NULL) 
dim(df_cnv)
res_604 <- plot_CNV_profile(df_cnv, clones= levels(df_cnv$clone),plttitle='Pt3')
res_604$cnv_plot
res_604$plg
dev.off()


datatag <- 'SA609'
copynumber_fn <- paste0(base_dir,datatag,'/total_merged_filtered_states_original.csv')
cellclone_fn <- paste0(input_dir, datatag, '_cell_clones.csv')
df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir,
                                 cellclone_fn, library_grouping_fn=NULL) 
dim(df_cnv)
res_609 <- plot_CNV_profile(df_cnv, clones= levels(df_cnv$clone),plttitle='Pt4')
res_609$cnv_plot
res_609$plg


datatag <- 'SA535'
copynumber_fn <- paste0(base_dir,'SA535/SA535_fitness/total_merged_filtered_states.csv')
cellclone_fn <- paste0(input_dir, datatag, '_cell_clones.csv')
df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir,
                                 cellclone_fn, library_grouping_fn=NULL) 
dim(df_cnv)
colnames(df_cnv)
unique(df_cnv$clone)
res_535 <- plot_CNV_profile(df_cnv, clones= levels(df_cnv$clone),plttitle='Pt5')
res_535$cnv_plot
dev.off()

datatag <- 'SA1035'
copynumber_fn <- paste0(base_dir,'SA1035_new_encode/SA1035_Tyler_v2/total_merged_filtered_states.csv')
cellclone_fn <- paste0(input_dir, datatag, '_cell_clones.csv')
df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir,
                                 cellclone_fn, library_grouping_fn=NULL) 
dim(df_cnv)
res_1035 <- plot_CNV_profile(df_cnv, clones= levels(df_cnv$clone),plttitle='Pt6')
res_1035$cnv_plot
dev.off()

# rel_heights = c(4.5,4,
#                 4.5,10.5,
#                 12.5,11.5)
plt_total <- cowplot::plot_grid(res_501$cnv_plot, res_609$cnv_plot, 
                                res_530$cnv_plot, res_535$cnv_plot, 
                                res_604$cnv_plot, res_1035$cnv_plot,
                                # align = 'v',
                                ncol = 2
                                ) + #labels = c('SA501','SA530','SA604','SA609','SA535','SA1035')
             theme(plot.background = element_rect(fill = "white", colour = "white"))

plt_total2 <- cowplot::plot_grid(plt_total, res_501$plg, ncol = 1, rel_heights=c(15,1)) + theme(plot.background = element_rect(fill = "white", colour = "white"))

plt_total
ggsave(paste0(save_dir,"Fig2_DLPtree_clonealign_correlation.pdf"),
       plot = plt_total2,
       height = 12,
       width = 7.5,
       useDingbats=F,
       dpi = 150
)

ggsave(paste0(save_dir,"Fig_SUPP2_segmentCNV.png"),
       plot = plt_total2,
       height = 8,
       width = 9,
       type = "cairo-png",
       dpi=300
)





save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
results_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/cell_clones/'
datatag <- 'SA501'

save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/dlp_cnv/medianCNV_MainFig2/'
if(!dir.exists(save_dir)){
        dir.create(save_dir)
}
copynumber_fn <- paste0(results_dir, datatag, '_total_merged_filtered_states.csv.gz')
cellclone_fn <- paste0(results_dir, datatag, '_cell_clones.csv.gz')
res_501 <- get_median_genotype_v2(datatag, save_dir, copynumber_fn, cellclone_fn,
                                  calcul_distance=T, distance_type='Manhattan')
# res_501 <- readRDS(paste0(save_dir, datatag, '_results_CNA_distance.rds'))

# saveRDS(res_501$hm, paste0(save_dir, datatag, '_results_CNA_distance_plt.rds'))
min(res_501$out_mtx$CNA_Distance)
max(res_501$out_mtx$CNA_Distance)
dim(res_501$dist_to_median)

results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA604_v2/'
datatag <- 'SA604'
res_SA604 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                    calcul_distance=T, distance_type='Manhattan')
saveRDS(res_SA604, paste0(save_dir, datatag, '_results_CNA_distance.rds'))
res_SA604 <- readRDS(paste0(save_dir, datatag, '_results_CNA_distance.rds'))
# min(res_SA604$out_mtx$CNA_Distance)
# max(res_SA604$out_mtx$CNA_Distance)

results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA530_v2/'
datatag <- 'SA530'
res_SA530 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                    calcul_distance=T, distance_type='Manhattan')
saveRDS(res_SA530, paste0(save_dir, datatag, '_results_CNA_distance.rds'))
res_SA530 <- readRDS(paste0(save_dir, datatag, '_results_CNA_distance_plt.rds'))


results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/'
datatag <- 'SA1035'
res_SA1035 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                     calcul_distance=T, distance_type='Manhattan')
saveRDS(res_SA1035, paste0(save_dir, datatag, '_results_CNA_distance_plt.rds'))
res_SA1035 <- readRDS(paste0(save_dir, datatag, '_results_CNA_distance_plt.rds'))
# min(res_SA1035$out_mtx$CNA_Distance)
# max(res_SA1035$out_mtx$CNA_Distance)


# chrs <- data.table::fread('/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/CN_profile/median_cnv.csv')
# length(unique(chrs$chr_desc))
results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609/'
datatag <- 'SA609'
base_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones/'
df <- data.table::fread(paste0(base_dir,'SA609_cisplatin_segments_total.csv')) #filtered_reads_total.csv.gz
dim(df)
head(df)
unique(df$clone_id)
df$chr_desc <- paste0(df$chr,'_',df$start,'_',df$end)
df$median_cn <- df$state
df$clone_id <- gsub('clone_','',df$clone_id)
stat_cnv <- df
res_SA609 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                    calcul_distance=T, distance_type='Manhattan')
# res_SA609 <- res
saveRDS(res_SA609$hm, paste0(save_dir, datatag, '_results_CNA_distance_plt.rds'))
res_SA609 <- readRDS(paste0(save_dir, datatag, '_results_CNA_distance_plt.rds'))
min(res_SA609$out_mtx$CNA_Distance)
max(res_SA609$out_mtx$CNA_Distance)


results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_fitness/'
datatag <- 'SA535'
res_SA535 <- get_median_genotype_v2(datatag, results_dir, copynumber_fn=NULL, cellclone_fn=NULL,
                                    calcul_distance=T, distance_type='Manhattan')
saveRDS(res_SA535, paste0(save_dir, datatag, '_results_CNA_distance.rds'))
res_SA535 <- readRDS(paste0(save_dir, datatag, '_results_CNA_distance.rds'))



res_501$out_mtx$datatag <- 'Pt1'
untreated_st <- '-Rx'
treated_st <- '+Rx +Rx'
res_501$out_mtx$treatment_status <- untreated_st
res_SA530$out_mtx$datatag <- 'Pt2'
res_SA530$out_mtx$treatment_status <- untreated_st
res_SA604$out_mtx$datatag <- 'Pt3'
res_SA604$out_mtx$treatment_status <- untreated_st

res_SA609$out_mtx$datatag <- 'Pt4'
res_SA609$out_mtx$treatment_status <- treated_st
res_SA535$out_mtx$datatag <- 'Pt5'
res_SA535$out_mtx$treatment_status <- treated_st
res_SA1035$out_mtx$datatag <- 'Pt6'
res_SA1035$out_mtx$treatment_status <- treated_st

stat <- tibble::tibble()
stat <- dplyr::bind_rows(list(res_501$out_mtx, res_SA530$out_mtx, res_SA604$out_mtx,
                              res_SA609$out_mtx, res_SA535$out_mtx, res_SA1035$out_mtx))
dim(stat)
head(stat)
unique(stat$treatment_status)
data.table::fwrite(stat, paste0(save_dir,'total_stat_median_CN_distance_6series.csv'))
save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609/CN_profile/'
stat <- data.table::fread(paste0(save_dir,'total_stat_median_CN_distance_6series.csv')) %>% as.data.frame()


stat$CNA_Distance <- as.numeric(stat$CNA_Distance)
# series <- c('SA501','SA530','SA604','SA609','SA535','SA1035')
# stat$datatag <- factor(stat$datatag, levels = series)
dim(stat)
untreated_st <- '-Rx'
treated_st <- '+Rx +Rx'
unique(stat$treatment_status)
stat$treatment_status <- ifelse(stat$treatment_status=="UnRx",untreated_st,treated_st)
p_manhattan <- ggplot(stat, aes(x = datatag, y=CNA_Distance, color = datatag)) +
        geom_boxplot() + #alpha = 0.6, size = 0.3
        facet_wrap(~factor(treatment_status, levels = c(untreated_st,treated_st)), scales = "free_x",
                   strip.position = "top",
                   # switch = "x",
                   nrow = 2, drop = T) + 
        theme_bw()+
        theme(legend.position = 'none',
              axis.title.x = element_blank(),
              strip.text.x = element_blank(),
              axis.text.x = element_text(size=9, family='Helvetica', angle = 30)) + 
        labs(x=NULL, y='CNA distance between clones')

# p_manhattan
saveRDS(p_manhattan, paste0(save_dir,"median_CN_distance_6series.rds"))
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
p_manhattan <- readRDS(paste0(save_dir,"median_CN_distance_6series.rds"))

png(paste0(save_dir,"median_CN_distance_6series.png"), height = 2*300, width=2*500, res = 2*72)
print(p)
dev.off()

# hm_ls <- res_501$hm + res_SA530$hm + res_SA604$hm + res_SA609$hm + res_SA535$hm + res_SA1035$hm
stat1 <- stat %>%
        # select(-SourceClone,-TargetClone)%>%
        group_by(datatag)%>%
        summarise(mean_distance=round(mean(CNA_Distance),2),
                  sd=round(sd(CNA_Distance),2))
stat1
series <- c('SA501','SA530','SA604','SA609','SA535','SA1035')
stat1 <- as.data.frame(stat1)
stat1$datatag <- factor(stat1$datatag, levels = series)
rownames(stat1) <- stat1$datatag
stat1 <- stat1[series,]
for(d in stat1$datatag){
        print(paste0('Series: ',d,'; mean CNA distance: ',stat1[d,'mean_distance'], '; sd CNA distance: ',stat1[d,'sd']))
}

stat2 <- stat1 %>%
        filter(datatag %in% c('SA609','SA535','SA1035'))
View(stat2)
stat2$idx <- paste0(stat2$datatag, stat2$SourceClone, stat2$TargetClone)
stat3 <- stat2 %>%
        filter(idx %in% c('SA609AH','SA609HA','SA535AG','SA535GA','SA1035HE','SA1035EH'))
View(stat3)

