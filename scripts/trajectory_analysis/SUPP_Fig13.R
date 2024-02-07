datatag <- 'SA609'
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
save_figs_dir <- paste0(save_dir,'figs_revision/')
legend_label='clone'
# color_lg_609 <- readRDS(paste0(save_figs_dir, gsub(' ','_',legend_label),'_plt.rds'))
pseudo_clone_609 <- readRDS(paste0(save_figs_dir,"ts_slingshot_out_wholedataset_clones_",datatag,".rds"))

datatag <- 'SA535'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
save_figs_dir <- paste0(save_dir,'figs_revision/')
legend_label='clone'
# color_lg_535 <- readRDS(paste0(save_figs_dir, gsub(' ','_',legend_label),'_plt.rds'))
pseudo_clone_535 <- readRDS(paste0(save_figs_dir,"ts_slingshot_out_wholedataset_clones_",datatag,".rds"))


datatag <- 'SA1035'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
save_figs_dir <- paste0(save_dir,'figs_revision/')
legend_label='clone'
# color_lg_1035 <- readRDS(paste0(save_figs_dir, gsub(' ','_',legend_label),'_plt.rds'))
pseudo_clone_1035 <- readRDS(paste0(save_figs_dir,"ts_slingshot_out_wholedataset_clones_",datatag,".rds"))

p_609 <- cowplot::plot_grid(pseudo_clone_609, NULL, ncol=2, rel_widths = c(1,0.3), hjust = 0)
p_1035 <- cowplot::plot_grid(pseudo_clone_1035, NULL, ncol=2, rel_widths = c(1,0.2), hjust = 0)
p_total <- cowplot::plot_grid(p_609, pseudo_clone_535, p_1035, 
                              ncol=1, rel_heights = c(1.1,0.9, 0.95))

# p_total <- cowplot::plot_grid(pseudo_clone_609, pseudo_clone_535, 
#                          color_lg_609, color_lg_535,
#                          label_size=12, hjust = 0,
#                          labels = c('a  Pt4 lineages - inferred clones',' b  Pt5 lineages - inferred clones'),
#                          ncol=2, rel_heights = c(1,0.6),
#                          rel_widths = c(1,1.2))
p_total
library(ggplot2)
ggsave(paste0(save_figs_dir,'revision_question38_Pt4_Pt5_Pt6_pseudo_clones.svg'),
       plot = p_total,
       height = 7.5,
       width = 5.5,
       # useDingbats=F,
       dpi=250)


# library(dplyr)
# base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/clonealign_plot/clonealign/SA1035-v6/'
# aligned <- data.table::fread(paste0(base_dir,'clonealign_plot/clonealign/SA1035-v6/SA1035X4XB02879.csv.gz'))
# dim(aligned)
# summary(as.factor(aligned$clone))
# 
# input_dir <- paste0(base_dir, 'materials/umap_figs/')
# datatag <- 'SA1035'
# umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap_filtered_outliers.csv.gz')) %>% as.data.frame()
# dim(umap_df)
# umap_df <- umap_df %>%
#    dplyr::filter(timepoint=='X4')
# 
# fns <- list.files(base_dir)
# total_cells <- 0
# for(f in fns){
#         tmp <- data.table::fread(paste0(base_dir, f))
#         total_cells <- total_cells + dim(tmp)[1]
# }
# total_cells



plot_colors_full <- function(cols, save_dir, datatag, legend_label='', nrow_plt=5){
        # clusters_use <- gtools::mixedsort(names(cols))
        # clusters_use <- c('UU','UUU','UUUU','UUUUU','UUT','UUTT','UUTTT','UUTTTTT','UUTU','UUTTU','UUTTTTU')
        # cols <- cols[clusters_use]
        my_font <- "Helvetica"
        color_df <- data.frame(color_code=cols, cluster=names(cols), vals=c(rep(1:length(cols),1)))
        # data.table::fwrite(color_df, paste0(save_dir, legend_label,'_metacolors.csv'))
        # color_df$cluster <- factor(color_df$cluster, levels = lvs)
        color_df$cluster <- factor(color_df$cluster, levels = color_df$cluster)
        p <- ggplot(color_df, aes(x=cluster, y=vals, fill=cluster)) +
                # geom_line(aes_string(color=plottype)) +  #,color=colorcode
                # geom_bar(size=2) + 
                geom_bar(stat="identity", width = 0.2)+
                # geom_point(size=7.5) +
                scale_fill_manual(values = cols) + labs(colour='') + 
                theme(legend.text=element_text(size=11, hjust = 0, family=my_font),
                      legend.title=element_text(size=11, hjust = 0.5, family=my_font),
                      panel.spacing = unit(c(0, 0, 0, 0), "null"),
                      legend.key=element_blank())#panel.margin = unit(c(0, 0, 0, 0), "null")
        p <- p + guides(fill = guide_legend(title = legend_label, 
                                             # nrow=nrow_plt, 
                                             ncol=1))
                                             # override.aes = list(size=1))) 
        
        lg <- cowplot::get_legend(p)
        pc <- cowplot::ggdraw() + cowplot::draw_plot(lg) #+ theme(plot.background = element_rect(fill = "white", colour = "white"))
        saveRDS(pc, paste0(save_dir, datatag,'_plt_legend.rds'))
       
        return(pc)
}
get_color_clones <- function(tag, color_fn){
        color_df <- data.table::fread(color_fn)
        # color_df <- data.table::fread(paste0(output_dir,'colorcode_total.csv'))
        color_df <- color_df %>%
                dplyr::filter(datatag==tag)
        
        if(dim(color_df)[1] > 0){
                cols_use <- color_df$colour
                names(cols_use) <- color_df$clone_id
                # cols_use['None'] <- '#D3D3D3'
                # if('None' %in% names(cols_use)){
                #     
                # }
                
        }else{
                # cols_use <- make_clone_palette(obs_clones)
                stop('Error, check color mapping')
        }
        return(cols_use)
}
save_dir <- save_figs_dir
base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
datatag <- 'SA609'
cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
lg_609 <- plot_colors_full(cols_use, save_dir, datatag, legend_label='Clone')

datatag <- 'SA535'
cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
lg_535 <- plot_colors_full(cols_use, save_dir, datatag, legend_label='Clone')

datatag <- 'SA1035'
cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
lg_1035 <- plot_colors_full(cols_use, save_dir, datatag, legend_label='Clone')

p_lg <- cowplot::plot_grid(lg_609, lg_535, lg_1035, ncol=3)
p_lg
ggsave(paste0(save_dir,'revision_question38_clone_lgs.svg'),
       plot = p_lg,
       height = 3,
       width = 2,
       # useDingbats=F,
       dpi=100)
