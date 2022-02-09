
# remove.packages("schnapps")
# library(devtools)
# devtools::install_github("shahcompbio/schnapps")
suppressPackageStartupMessages({
  library(schnapps)
  library(dplyr)
  library(cowplot)
  library(RColorBrewer)
})

source_dir <- '/home/htran/Projects/hakwoo_project/rscript/chnapp_haplo/'
source(paste0(source_dir,'schnapp_utils.R'))
save_dir <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/schnapps_out/"
prefix <- 'total'
prefix_filtered <- 'total_filtered'
# results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'

tyler_results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler_left_padding/'
# cell_clones <- paste0(results_dir,'tree_cut_out/cell_clones_if_0.02_af_0.75_p0.75_e0.04.csv')

cell_clones <- paste0(tyler_results_dir, 'cell_clones.csv')
clones_labels_df <- read.csv(cell_clones, check.names = F, stringsAsFactors = FALSE)

# Read clones
# clones_labels_df <- read.csv(paste0(results_dir,'cell_clones.csv'),header=T,check.names=F)
# clones_labels_df <- read.csv(paste0(results_dir,'tree_cut_out/cell_clones_if_0.02_af_0.75_p0.75_e0.04.csv'),header=T,check.names=F)
# unique(clones_labels_df$clone_id)
# ex_clones <- c('I','J')
# clones_labels_df <- clones_labels_df[!clones_labels_df$clone_id %in% ex_clones, ,drop = TRUE]
print(dim(clones_labels_df))
print(unique(clones_labels_df$clone_id))
rownames(clones_labels_df) <- clones_labels_df$cell_id


# CNBAF <- combine_BAF(save_dir, prefix, prefix_filtered)


# Run allele specific copy number
# run_alleleSC(CNBAF=NULL, save_dir, prefix)
# fn <- paste0(save_dir,prefix,'_CNBAF.rds')
# CNBAF <- readRDS(file = fn)
# p1 <- plotCNBAF(CNBAF)
# png(paste0(save_dir,paste0(prefix,"_BAF_CNstate.png")), height = 2*460, width=2*900,res = 2*72)
# print(p1)
# dev.off()


# Visualization
# output of CNbins contain width?
fn <- paste0(save_dir,prefix,'_alleleCNbins_out.rds')
CNbins <- readRDS(fn)
main_source <- '/home/htran/storage/install_software/schnapps/R/'
source(paste0(main_source,'heatmap_plot.R'))
# visualize_hm_clustering(CNbins,save_dir,prefix)
# plots <- visualize_hm_each_clone(CNbins,clones_labels_df,save_dir,prefix)
save_dir <- paste0(tyler_results_dir,'schnapps_out/')
if (!file.exists(save_dir)){
  dir.create(save_dir)
}
# visualize_hm_each_clone_plot(plots,save_dir,prefix)
visualize_hm_tree(CNbins, clones_labels_df, tyler_results_dir, save_dir, prefix='total')

sum(is.na(CNbins$copy))

dim(CNbins)







# Code example from Marc
# CNclones <- ascn_SA609 %>%
#   left_join(., select(cells_SA609, clone_id, cell_id)) %>%
#   group_by(chr, start, end, clone_id) %>%
#   summarise(state = schnapps::Mode(state),
#             state_phase = schnapps::Mode(state_phase),
#             state_min = schnapps::Mode(state_min),
#             BAF = median(BAF),
#             copy = median(copy)) %>%
#   ungroup() %>%
#   mutate(cell_id = paste0("Clone ", clone_id))
# 
# for (cl in unique(CNclones$cell_id)){
#   cat('###', cl,' \n')
#   print(plotCNprofileBAF(CNclones, cellid = cl))
#   cat(' \n \n')
# }


