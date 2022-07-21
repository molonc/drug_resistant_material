source('/home/htran/Projects/farhia_project/rnaseq/drug_manuscript/viz_umap_figs/viz_umaps.R')
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/'
output_dir <- paste0(base_dir,'umap_figs/figs/')

datatag <- 'SA501'
library_grouping_fn <- paste0(base_dir,'cell_clones/',datatag, '_DLP_library_groupings.csv')
cell_clones_fn <- paste0(base_dir,'cell_clones/',datatag, '_cell_clones.csv')
cell_clones <- data.table::fread(cell_clones_fn) %>% as.data.frame()
dim(cell_clones)
cols_use <- make_clone_palette(unique(cell_clones$clone_id))

obs_treatment <- 'UnRx'
obs_passage <- 'X2'
cell_clones$timepoint <- 'X2'

## version 1
# plottitle <- obs_treatment
# res_501 <- plot_fill_barplot(cell_clones, cols_use, 
#                              plottitle, output_dir, datatag, 
#                              plotlegend=F)



## version 2
plottitle <- 'Pt1-DLP+'
cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
res_barDLP_501 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                            datatag, plottitle=NULL, plotlegend=F)

res_barDLP_501$p


## Building color maps here
# colorcode <- tibble::tibble()
# df <- tibble(datatag=datatag, clone_id=names(cols_use), colour=cols_use)
# colorcode <- dplyr::bind_rows(colorcode, df)
# data.table::fwrite(colorcode, paste0(output_dir,'colorcode_total.csv'))



datatag <- 'SA530'
cell_clones_fn <- paste0(base_dir,'cell_clones/',datatag, '_cell_clones.csv')
library_grouping_fn <- paste0(base_dir,'cell_clones/',datatag, '_DLP_library_groupings.csv')
cell_clones <- data.table::fread(paste0(base_dir,'cell_clones/',datatag, '_cell_clones.csv')) %>% as.data.frame()
dim(cell_clones)
head(cell_clones)
cols_use <- make_clone_palette(unique(cell_clones$clone_id))
obs_treatment <- 'UnRx'
obs_passage <- 'X3'
cell_clones$timepoint <- 'X3'

# plottitle <- obs_treatment
# cell_clones_extracted <- cell_clones %>% 
#   dplyr::filter(timepoint %in% obs_passage)
# res_530 <- plot_fill_barplot(cell_clones, cols_use, 
#                              plottitle, output_dir, datatag, 
#                              plotlegend=F)

# df <- tibble(datatag=datatag, clone_id=names(cols_use), colour=cols_use)
# colorcode <- dplyr::bind_rows(colorcode, df)
# data.table::fwrite(colorcode, paste0(output_dir,'colorcode_total.csv'))
# unique(colorcode$datatag)

## version 2
plottitle <- 'Pt2-DLP+'
cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
res_barDLP_530 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                          datatag, plottitle=NULL, plotlegend=F)

res_barDLP_530$p
res_barDLP_530$plg


datatag <- 'SA604'
cell_clones_fn <- paste0(base_dir,'cell_clones/',datatag, '_cell_clones.csv')
library_grouping_fn <- paste0(base_dir,'cell_clones/',datatag, '_DLP_library_groupings.csv')
cell_clones <- data.table::fread(paste0(base_dir,'cell_clones/',datatag, '_cell_clones.csv')) %>% as.data.frame()
dim(cell_clones)
head(cell_clones)
unique(cell_clones$timepoint)
# cols_use <- make_clone_palette(unique(cell_clones$clone_id))
cols_use <- get_color_clones(datatag)
# cell_clones$timepoint <- stringr::str_sub(cell_clones$cell_id,6,7)
# unique(cell_clones$timepoint)
# cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
# unique(cell_clones$sample_id)
# unique(cell_clones$clone_id)
# 
# data.table::fwrite(cell_clones, cell_clones_fn)
# obs_treatment <- 'UnRx'
# plottitle <- obs_treatment
# cell_clones_extracted <- cell_clones %>% 
#   dplyr::filter(timepoint %in% obs_passage)
# res_604 <- plot_fill_barplot(cell_clones, cols_use, 
#                              plottitle, output_dir, datatag, 
#                              plotlegend=F)

## version 2
plottitle <- 'Pt3-DLP+'
cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
res_barDLP_604 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                          datatag, plottitle=NULL, plotlegend=F)
res_barDLP_604$p
res_barDLP_604$plg
## Building color code
# df <- tibble(datatag=datatag, clone_id=names(cols_use), colour=cols_use)
# colorcode <- dplyr::bind_rows(colorcode, df)
# data.table::fwrite(colorcode, paste0(output_dir,'colorcode_total.csv'))
# unique(colorcode$datatag)


datatag <- 'SA609'
cell_clones <- data.table::fread(paste0(base_dir,'cell_clones/',datatag, '_cell_clones_total.csv')) %>% as.data.frame()
dim(cell_clones)
cols_use <- make_clone_palette(unique(cell_clones$clone_id))
# cols_use <- get_color_clones(datatag, obs_clones=NULL)
# cell_clones$timepoint <- stringr::str_sub(cell_clones$cell_id,6,8)
# unique(cell_clones$timepoint)
# cell_clones$timepoint <- gsub('X$','',cell_clones$timepoint)
# cell_clones$mouse_id <- get_sample_id(cell_clones$cell_id)

# if(datatag=='SA609'){
obs_samples <- c('SA609X3XB01584','SA609X4XB03080','SA609X5XB03223',
                 'SA609X6XB03447','SA609X7XB03554','SA609X4XB003083',
                 'SA609X5XB03230','SA609X6XB03404','SA609X7XB03505',
                 'SA609X5XB03231','SA609X6XB03401','SA609X7XB03510')
  # sids <- unique(cell_clones$sample)
  # sum(sids %in% obs_samples)
#   # t <- sids[!sids %in% obs_samples]
#   # x <- grep('SA609X3*',t,value=T)
#   
#   cell_clones <- cell_clones %>%
#     dplyr::filter(mouse_id %in% obs_samples)
#   print('Replicates excluded!!!')
#   print(dim(cell_clones))
# }

# unique(cell_clones$clone_id)

# meta_cells <- data.table::fread(paste0(dirname(base_dir),'/SA609_rna/snakemake_10x/SA609_10x.csv')) %>% as.data.frame()
# dim(meta_cells)
# meta_cells$mouse_id
# head(meta_cells)
# cell_clones <- cell_clones %>% left_join(meta_cells, by=c('mouse_id'))
unique(cell_clones$treatmentstr) #"R" "H" "U" "M"
obs_treatment <- 'UnRx'
plottitle <- obs_treatment
cell_clones_extracted <- cell_clones %>%
  dplyr::filter(treatmentstr=="U" & timepoint %in% c('X4','X5','X6','X7'))
res_SA609_UnRx <- plot_fill_barplot(cell_clones_extracted, cols_use, 
                                    plottitle, output_dir, datatag, 
                                    plotlegend=F)


obs_treatment <- 'Rx'
plottitle <- obs_treatment
cell_clones_extracted <- cell_clones %>%
  dplyr::filter(treatmentstr=="R" & timepoint %in% c('X4','X5','X6','X7'))
res_SA609_Rx <- plot_fill_barplot(cell_clones_extracted, cols_use, 
                                  plottitle, output_dir, datatag, 
                                  plotlegend=F)

obs_treatment <- 'RxH'
plottitle <- obs_treatment
cell_clones_extracted <- cell_clones %>%
  dplyr::filter(treatmentstr=="H" & timepoint %in% c('X4','X5','X6','X7'))
res_SA609_RxH <- plot_fill_barplot(cell_clones_extracted, cols_use, 
                                   plottitle, output_dir, datatag, 
                                   plotlegend=F)

df <- tibble(datatag=datatag, clone_id=names(cols_use), colour=cols_use)
colorcode <- dplyr::bind_rows(colorcode, df)
data.table::fwrite(colorcode, paste0(output_dir,'colorcode_total.csv'))
unique(colorcode$datatag)


datatag <- 'SA535'
cell_clones <- data.table::fread(paste0(base_dir,'cell_clones/',datatag, '_cell_clones.csv')) %>% as.data.frame()
dim(cell_clones)
head(cell_clones)
grouping_df <- data.table::fread(paste0(base_dir,'cell_clones/',datatag, '_DLP_library_groupings.csv')) %>% as.data.frame()
dim(grouping_df)
colnames(grouping_df)
head(grouping_df)
pdx <- unique(grouping_df$PDX)
pdx <- pdx[!grepl('CX5461', pdx)]
grouping_df <- grouping_df %>%
  dplyr::filter(PDX %in% pdx)%>%
  dplyr::rename(library_id=grouping, treatmentSt=treatment_st)
grouping_df <- grouping_df %>%
  dplyr::select(library_id, treatmentSt, passage, sample)%>%
  dplyr::rename(timepoint=passage)

cell_clones$library_id <- get_library_id(cell_clones$cell_id)

cell_clones <- cell_clones %>% left_join(grouping_df, by=c('library_id','sample'))


cols_use <- make_clone_palette(unique(cell_clones$clone_id))


obs_treatment <- 'UnRx'
plottitle <- obs_treatment
conds <- c('U',grep('*UU$',unique(cell_clones$treatmentSt), value=T))

cell_clones_extracted <- cell_clones %>%
  dplyr::filter(treatmentSt %in% conds)
dim(cell_clones_extracted)
res_SA535_UnRx <- plot_fill_barplot(cell_clones_extracted, cols_use, 
                                    plottitle, output_dir, datatag, 
                                    plotlegend=F)
res_SA535_UnRx$p

obs_treatment <- 'Rx'
plottitle <- obs_treatment
conds <- grep('*T$',unique(cell_clones$treatmentSt), value=T)
cell_clones_extracted <- cell_clones %>%
  dplyr::filter(treatmentSt %in% conds)
dim(cell_clones_extracted)
res_SA535_Rx <- plot_fill_barplot(cell_clones_extracted, cols_use, 
                                  plottitle, output_dir, datatag, 
                                  plotlegend=F)

obs_treatment <- 'RxH'
plottitle <- obs_treatment
conds <- grep('*TU$',unique(cell_clones$treatmentSt), value=T)
cell_clones_extracted <- cell_clones %>%
  dplyr::filter(treatmentSt %in% conds)
dim(cell_clones_extracted)
res_SA535_RxH <- plot_fill_barplot(cell_clones_extracted, cols_use, 
                                   plottitle, output_dir, datatag, 
                                   plotlegend=F)





base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/'
output_dir <- paste0(base_dir,'umap_figs/figs/')
datatag <- 'SA1035'
cell_clones <- data.table::fread(paste0(base_dir,'cell_clones/',datatag, '_cell_clones.csv')) %>% as.data.frame()
dim(cell_clones)
head(cell_clones)
grouping_df <- data.table::fread(paste0(base_dir,'cell_clones/',datatag, '_DLP_library_groupings.csv')) %>% as.data.frame()
dim(grouping_df)
colnames(grouping_df)
# head(grouping_df)
pdx <- unique(grouping_df$PDX)
grouping_df$passage[1]
grouping_df <- grouping_df %>%
  dplyr::rename(library_id=grouping, treatmentSt=treatment_st, timepoint=passage)%>%
  dplyr::select(library_id, treatmentSt, sample, timepoint)

cell_clones$library_id <- get_library_id(cell_clones$cell_id)

cell_clones <- cell_clones %>% left_join(grouping_df, by=c('library_id','sample'))


# cols_use <- make_clone_palette(unique(cell_clones$clone_id))
cols_use <- get_color_clones(datatag)

obs_treatment <- 'UnRx'
plottitle <- obs_treatment
conds <- c('U',grep('*UU$',unique(cell_clones$treatmentSt), value=T))
conds
cell_clones_extracted <- cell_clones %>%
  dplyr::filter(treatmentSt %in% conds)
dim(cell_clones_extracted)
unique(cell_clones_extracted$treatmentSt)
res_SA1035_UnRx <- plot_fill_barplot(cell_clones_extracted, cols_use, 
                                     plottitle, output_dir, datatag, 
                                     plotlegend=F)

res_SA1035_UnRx$p
obs_treatment <- 'Rx'
plottitle <- obs_treatment
conds <- grep('*T$',unique(cell_clones$treatmentSt), value=T)
cell_clones_extracted <- cell_clones %>%
  dplyr::filter(treatmentSt %in% conds)
dim(cell_clones_extracted)
conds
unique(cell_clones_extracted$treatmentSt)
res_SA1035_Rx <- plot_fill_barplot(cell_clones_extracted, cols_use, 
                                   plottitle, output_dir, datatag, 
                                   plotlegend=F)
res_SA1035_Rx$p
obs_treatment <- 'RxH'
plottitle <- obs_treatment
conds <- grep('*TU$',unique(cell_clones$treatmentSt), value=T)
cell_clones_extracted <- cell_clones %>%
  dplyr::filter(treatmentSt %in% conds)
dim(cell_clones_extracted)
conds
unique(cell_clones_extracted$treatmentSt)
res_SA1035_RxH <- plot_fill_barplot(cell_clones_extracted, cols_use, 
                                    plottitle, output_dir, datatag, 
                                    plotlegend=F)

res_SA1035_RxH$p




##########################
### 10x 



base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/'
output_dir <- paste0(base_dir,'umap_figs/hvg_umaps/')
datatag <- 'SA609'
plt_ls <- list()
for(nb_hvg in c(5000,10000)){
  umap_df <- data.table::fread(paste0(output_dir, datatag, "_norm_umap_",nb_hvg,".csv")) %>% as.data.frame()
  base_name <- paste0(datatag,': ',nb_hvg,' input genes, 50 input PCs, UMAP plot')
  p <- viz_rd(umap_df, base_name, xstring='UMAP_1', ystring='UMAP_2', plottype='treatmentSt')
  plt_ls[[gsub(' ','',base_name)]] <- p
  
  pca_df <- data.table::fread(paste0(output_dir, datatag, "_norm_pca_",nb_hvg,".csv")) %>% as.data.frame()
  base_name <- paste0(datatag,': ',nb_hvg,' input genes, 50 input PCs, PCA plot')
  p <- viz_rd(pca_df, base_name, xstring='PC_1', ystring='PC_2', plottype='treatmentSt')
  plt_ls[[gsub(' ','',base_name)]] <- p
}
nb_hvg <- 10000
umap_df <- data.table::fread(paste0(output_dir, datatag, "_norm_umap_",nb_hvg,".csv")) %>% as.data.frame()
base_name <- paste0(datatag,': ',nb_hvg,' input genes, 30 input PCs, UMAP plot')
p <- viz_rd(umap_df, base_name, xstring='UMAP_1', ystring='UMAP_2', plottype='treatmentSt')
plt_ls[[gsub(' ','',base_name)]] <- p
dim(umap_df)
pca_df <- data.table::fread(paste0(output_dir, datatag, "_norm_pca_",nb_hvg,".csv")) %>% as.data.frame()
base_name <- paste0(datatag,': ',nb_hvg,' input genes, 30 input PCs, PCA plot')
p <- viz_rd(pca_df, base_name, xstring='PC_1', ystring='PC_2', plottype='treatmentSt')
plt_ls[[gsub(' ','',base_name)]] <- p


prd <- cowplot::plot_grid(plt_ls$`SA609:5000inputgenes,50inputPCs,UMAPplot`, 
                          plt_ls$`SA609:10000inputgenes,50inputPCs,UMAPplot`, 
                          plt_ls$`SA609:10000inputgenes,30inputPCs,UMAPplot`,
                          plt_ls$`SA609:5000inputgenes,50inputPCs,PCAplot`, 
                          plt_ls$`SA609:10000inputgenes,50inputPCs,PCAplot`, 
                          plt_ls$`SA609:10000inputgenes,30inputPCs,PCAplot`, 
                          align = 'hv', ncol=3)
png(paste0(output_dir,"UMAPs_PCAs_SA609.png"), height = 2*700, width=2*1200,res = 2*72)
print(prd)
dev.off()

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/'
output_dir <- paste0(base_dir,'umap_figs/hvg_umaps/')
datatag <- 'SA535'
plt_ls <- list()
for(nb_hvg in c(5000,10000)){
  umap_df <- data.table::fread(paste0(output_dir, datatag, "_norm_umap_",nb_hvg,".csv")) %>% as.data.frame()
  base_name <- paste0(datatag,': ',nb_hvg,' input genes, 30 input PCs, UMAP plot')
  p <- viz_rd(umap_df, base_name, xstring='UMAP_1', ystring='UMAP_2', plottype='treatmentSt')
  plt_ls[[gsub(' ','',base_name)]] <- p
  
  pca_df <- data.table::fread(paste0(output_dir, datatag, "_norm_pca_",nb_hvg,".csv")) %>% as.data.frame()
  base_name <- paste0(datatag,': ',nb_hvg,' input genes, 30 input PCs, PCA plot')
  p <- viz_rd(pca_df, base_name, xstring='PC_1', ystring='PC_2', plottype='treatmentSt')
  plt_ls[[gsub(' ','',base_name)]] <- p
}

prd <- cowplot::plot_grid(plt_ls$`SA535:5000inputgenes,30inputPCs,UMAPplot`, 
                          plt_ls$`SA535:10000inputgenes,30inputPCs,UMAPplot`,
                          plt_ls$`SA535:5000inputgenes,30inputPCs,PCAplot`, 
                          plt_ls$`SA535:10000inputgenes,30inputPCs,PCAplot`, 
                          align = 'hv', ncol=2)
png(paste0(output_dir,"UMAPs_PCAs_",datatag,".png"), height = 2*700, width=2*850,res = 2*72)
print(prd)
dev.off()

# Just for debug
viz_rd <- function(df, base_name, xstring='UMAP_1', ystring='UMAP_2', plottype='treatmentSt'){
  # df_text <- df %>% 
  #   group_by(treatmentSt) %>% 
  #   summarise(avg_x = mean(PC_1),avg_y = mean(PC_2)) %>%
  #   ungroup()
  if(!plottype %in% colnames(df) & 'treat' %in% colnames(df)){
    df[,plottype] <- df[,'treat']
  }
  df_text <- df %>% 
    group_by(!!sym(plottype)) %>% 
    summarise(avg_x = mean(!!sym(xstring)),avg_y = mean(!!sym(ystring))) %>%
    ungroup()
  
  p <- ggplot2::ggplot(df, aes_string(x=xstring, y=ystring, color=plottype)) + 
    geom_point(size=0.5) + 
    annotate("text", 
             x = df_text$avg_x, 
             y = df_text$avg_y,
             label = df_text$treatmentSt) +
    theme_bw() + 
    labs(title = base_name)
  # p
  return(p)
}  

