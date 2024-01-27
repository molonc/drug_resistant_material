library(dplyr)

script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/'
source(paste0(script_dir,'drug_manuscript/cnv_viz_utils.R'))


results_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/cell_clones/'

save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/dlp_cnv/medianCNV_MainFig2/'
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}


res_total <- list()

series <- c('SA501','SA530','SA604','SA535','SA1035')

for(datatag in series){
  copynumber_fn <- paste0(results_dir, datatag, '_total_merged_filtered_states.csv.gz')
  copynumber_fn2 <- paste0(results_dir, datatag, '_total_merged_filtered_states_longformat.csv.gz')
  
  cellclone_fn <- paste0(results_dir, datatag, '_cell_clones.csv.gz')
  tmp <- get_median_genotype_v2(datatag, save_dir, copynumber_fn, cellclone_fn,
                                    calcul_distance=T, distance_type='Manhattan')
  # res_501 <- readRDS(paste0(save_dir, datatag, '_results_CNA_distance.rds'))
  
  # saveRDS(res_501$hm, paste0(save_dir, datatag, '_results_CNA_distance_plt.rds'))
  # print(min(res$out_mtx$CNA_Distance))
  # print(max(res$out_mtx$CNA_Distance))
  # print(dim(res$dist_to_median))
  res_total[[datatag]] <- tmp
}
rm(tmp)

datatag <- 'SA609'
copynumber_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/cnv/sa609_total_merged_filtered_states.csv'
# copynumber_fn <- paste0(results_dir, datatag, '_total_merged_filtered_states.csv.gz')
cellclone_fn <- paste0(results_dir, datatag, '_cell_clones.csv.gz')
res <- get_median_genotype_v2(datatag, save_dir, copynumber_fn, cellclone_fn,
                              calcul_distance=T, distance_type='Manhattan')
# res_501 <- readRDS(paste0(save_dir, datatag, '_results_CNA_distance.rds'))

# saveRDS(res_501$hm, paste0(save_dir, datatag, '_results_CNA_distance_plt.rds'))
print(min(res$out_mtx$CNA_Distance))
print(max(res$out_mtx$CNA_Distance))
print(dim(res$dist_to_median))
pt4_df <- res$out_mtx
dim(pt4_df)
head(pt4_df)
pt4_df$CNA_Distance <- as.numeric(pt4_df$CNA_Distance)
median_distance=round(median(pt4_df$CNA_Distance),2)
mean_distance=round(mean(pt4_df$CNA_Distance),2)
sd_distance=round(sd(pt4_df$CNA_Distance),2)
stat1 <- pt4_df %>%
  # select(-SourceClone,-TargetClone)%>%
  # group_by(datatag)%>%
  summarise(median_distance=round(median(CNA_Distance),2),
            sd=round(sd(CNA_Distance),2))
stat1


res_total[[datatag]] <- res

# res$hm
series <- c('SA501','SA530','SA604','SA609','SA535','SA1035')
pts_lb <- c('Pt1','Pt2','Pt3','Pt4','Pt5','Pt6') 
names(pts_lb) <- series

ts <- c('UnRx','UnRx','UnRx','Rx','Rx','Rx') 
names(ts) <- series
stat <- tibble::tibble()
for(datatag in series){
  res_total[[datatag]]$out_mtx$datatag <- pts_lb[datatag]
  res_total[[datatag]]$out_mtx$treatment_status <- ts[datatag]
  stat <- dplyr::bind_rows(stat, res_total[[datatag]]$out_mtx)
}

dim(stat)
head(stat)
unique(stat$treatment_status)
unique(stat$datatag)
data.table::fwrite(stat, paste0(save_dir,'total_stat_median_CN_distance_6series.csv.gz'))

# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# p <- readRDS(paste0(save_dir,'/manuscript/clonealign_plot/median_CN_distance_6series.rds'))
# stat <- p$data %>% as.data.frame()
# # stat <- data.table::fread(paste0(save_dir,'total_stat_median_CN_distance_6series.csv')) %>% as.data.frame()
# dim(stat)
# head(stat)
# t <- stat %>%
#   dplyr::filter(datatag=='Pt6') %>%
#   # dplyr::group_by(datatag) %>%
#   dplyr::summarise(median_dist=round(median(CNA_Distance),2),
#                    sd_dist=round(sd(CNA_Distance),2))
# dim(t)
# t
unique(t$datatag)
stat$CNA_Distance <- as.numeric(stat$CNA_Distance)
# series <- c('SA501','SA530','SA604','SA609','SA535','SA1035')
# stat$datatag <- factor(stat$datatag, levels = series)
dim(stat)
unique(stat$treatment_status)
# stat$treatment_status <- ifelse(stat$treatment_status=="UnRx",untreated_st,treated_st)
p_manhattan <- ggplot(stat, aes(x = datatag, y=CNA_Distance, color = datatag)) +
  geom_boxplot() + #alpha = 0.6, size = 0.3
  facet_wrap(~factor(treatment_status, levels = c("UnRx","Rx")), scales = "free_x",
             strip.position = "top",
             # switch = "x",
             nrow = 2, drop = T) + 
  theme_bw()+
  theme(legend.position = 'none',
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=11, family='Helvetica', color='black', angle = 20),
        axis.title.x = element_text(size=11, family='Helvetica', color='black')) + #
  labs(x=NULL, y='CNA distances between clones')

p_manhattan
output_dir <- "/home/htran/Projects/farhia_project/drug_resistant_material/materials/umap_figs/main_fig2/"
saveRDS(p_manhattan, paste0(output_dir,"median_CN_distance_6series.rds"))
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
# p_manhattan <- readRDS(paste0(save_dir,"median_CN_distance_6series.rds"))

# png(paste0(save_dir,"median_CN_distance_6series.png"), height = 2*300, width=2*500, res = 2*72)
# print(p)
# dev.off()

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


