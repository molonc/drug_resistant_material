## Plot raw data features 
ylabel <- '# genes non-zero counts'
suffix1 <- ": raw data \n Genes with non-zero counts"
suffix2 <- ": raw data \n Total counts"
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA1035'
output_dir <- paste0(input_dir,'rnaseq_v6/',datatag, '-v6/')
p1035_raw <- readRDS(paste0(output_dir,datatag,"_raw_data_features_v2.rds"))


meta_cells <- read.csv(paste0(output_dir,datatag,"_raw_data_features_v2.csv"), check.names=F, stringsAsFactors=F)
meta_cells$sample[1]
# Process treatment status first
meta_cells$tm_st <- ifelse(grepl('TU$',meta_cells$treatmentSt),'RxH',
                           ifelse(grepl('T$',meta_cells$treatmentSt),'Rx','UnRx'))
unique(meta_cells$tm_st)
meta_cells$series <- paste0(meta_cells$timepoint,' ',meta_cells$tm_st)
meta_cells <- meta_cells %>%
  dplyr::group_by(series) %>%
  dplyr::mutate(nb_cells = n())
meta_cells$sample_feature <- paste0(meta_cells$timepoint,' ',meta_cells$tm_st, ': ',meta_cells$nb_cells,' cells')

# Process chip infos
meta_cells$series <- factor(meta_cells$series, levels = gtools::mixedsort(unique(meta_cells$series)))
meta_cells$sample_feature <- factor(meta_cells$sample_feature, levels = gtools::mixedsort(unique(meta_cells$sample_feature)))
meta_cells_1035 <- meta_cells
p <- plot_variation_function_v3(meta_data=meta_cells, xstring="series", ystring="total_features_by_counts",
                                plottype="sample_feature", plottitle=paste0(datatag, suffix1),
                                xlabel=NULL, ylabel=NULL, lg_pos="right")


p_c <- plot_variation_function_v3(meta_data=meta_cells, xstring="series", ystring="total_counts",
                                  plottype="sample_feature", plottitle=paste0(datatag, suffix2),
                                  xlabel=NULL, ylabel=NULL,lg_pos='none')

p1035_raw <- cowplot::plot_grid(p, p_c, ncol=2, align = 'h', rel_widths = c(4.2,3))
p1035_raw
datatag <- gsub(' ','_',datatag)
saveRDS(p1035_raw, paste0(output_dir,datatag,"_raw_data_features_v2.rds"))
png(paste0(output_dir,datatag,"_raw_data_features.png"), height = 2*420, width=2*860,res = 2*72)
print(p1035_raw)
dev.off()


input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA609'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
output_dir <- save_dir

meta_cells <- read.csv(paste0(save_dir,datatag,"_raw_data_features_v2.csv"), check.names=F, stringsAsFactors=F)
meta_cells$sample[1]
dim(meta_cells)
# Process treatment status first
meta_cells$tm_st <- ifelse(grepl('TU$',meta_cells$treatmentSt),'RxH',
                           ifelse(grepl('T$',meta_cells$treatmentSt),'Rx','UnRx'))

meta_cells$series <- paste0(meta_cells$timepoint,' ',meta_cells$tm_st)
meta_cells <- meta_cells %>%
  dplyr::group_by(series) %>%
  dplyr::mutate(nb_cells = n())
meta_cells$sample_feature <- paste0(meta_cells$timepoint,' ',meta_cells$tm_st, ': ',meta_cells$nb_cells,' cells')

# Process chip infos
meta_cells$series <- factor(meta_cells$series, levels = gtools::mixedsort(unique(meta_cells$series)))
meta_cells$sample_feature <- factor(meta_cells$sample_feature, levels = gtools::mixedsort(unique(meta_cells$sample_feature)))
meta_cells_609 <- meta_cells
p <- plot_variation_function_v3(meta_data=meta_cells, xstring="series", ystring="total_features_by_counts",
                                plottype="sample_feature", plottitle=paste0(datatag, suffix1),
                                xlabel=NULL, ylabel=NULL,lg_pos="right")


p_c <- plot_variation_function_v3(meta_data=meta_cells, xstring="series", ystring="total_counts",
                                  plottype="sample_feature", plottitle=paste0(datatag, suffix2),
                                  xlabel=NULL, ylabel=NULL,lg_pos='none')

p609_raw <- cowplot::plot_grid(p, p_c, ncol=2, align = 'h', rel_widths = c(4.2,3))

datatag <- gsub(' ','_',datatag)
# p609_raw <- readRDS(paste0(output_dir,datatag,"_raw_data_features_v2.rds"))
saveRDS(p609_raw, paste0(output_dir,datatag,"_raw_data_features_v2.rds"))
png(paste0(output_dir,datatag,"_raw_data_features.png"), height = 2*420, width=2*860,res = 2*72)
print(p609_raw)
dev.off()

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA535_cisplatin'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
output_dir <- save_dir
rm(meta_cells)
meta_cells <- read.csv(paste0(save_dir,datatag,"_raw_data_features_v2.csv"), check.names=F, stringsAsFactors=F)
meta_cells$sample[1]
dim(meta_cells)
# Process treatment status first
meta_cells$tm_st <- ifelse(grepl('TU$',meta_cells$treatmentSt),'RxH',
                           ifelse(grepl('T$',meta_cells$treatmentSt),'Rx','UnRx'))
meta_cells$series <- paste0(meta_cells$timepoint,' ',meta_cells$tm_st)

meta_cells <- meta_cells %>%
  dplyr::group_by(series) %>%
  dplyr::mutate(nb_cells = n())
meta_cells$sample_feature <- paste0(meta_cells$timepoint,' ',meta_cells$tm_st, ': ',meta_cells$nb_cells,' cells')

meta_cells$series <- factor(meta_cells$series, levels = gtools::mixedsort(unique(meta_cells$series)))
meta_cells$sample_feature <- factor(meta_cells$sample_feature, levels = gtools::mixedsort(unique(meta_cells$sample_feature)))
meta_cells_535_cis <- meta_cells
p <- plot_variation_function_v3(meta_data=meta_cells, xstring="series", ystring="total_features_by_counts",
                                plottype="sample_feature", plottitle=paste0('SA535:Cisplatin', suffix1),
                                xlabel=NULL, ylabel=NULL,lg_pos="right")


p_c <- plot_variation_function_v3(meta_data=meta_cells, xstring="series", ystring="total_counts",
                                  plottype="sample_feature", plottitle=paste0('SA535:Cisplatin', suffix2),
                                  xlabel=NULL, ylabel=NULL,lg_pos='none')

p535_cis_raw <- cowplot::plot_grid(p, p_c, ncol=2, align = 'h', rel_widths = c(4.2,3))
datatag <- gsub(' ','_',datatag)
saveRDS(p535_cis_raw, paste0(output_dir,datatag,"_raw_data_features_v2.rds"))
p535_cis_raw <- readRDS(paste0(output_dir,datatag,"_raw_data_features_v2.rds"))
png(paste0(output_dir,datatag,"_raw_data_features.png"), height = 2*420, width=2*860,res = 2*72)
print(p535_cis_raw)
dev.off()

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA535_CX5461'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
output_dir <- save_dir
rm(meta_cells)
meta_cells <- read.csv(paste0(save_dir,datatag,"_raw_data_features_v2.csv"), check.names=F, stringsAsFactors=F)
meta_cells$sample[1]
dim(meta_cells)
# Process treatment status first
meta_cells$tm_st <- ifelse(grepl('XU$',meta_cells$treatmentSt),'RxH',
                           ifelse(grepl('X$',meta_cells$treatmentSt),'Rx','UnRx'))
meta_cells$tm_st[1]


# t <- unique(meta_cells$treatmentSt)
# t[grepl('TU$',t)]
# t[grepl('T$',t)]
# t[!grepl('TU$',t) & !grepl('T$',t)]

meta_cells$series <- paste0(meta_cells$timepoint,' ',meta_cells$tm_st)
meta_cells <- meta_cells %>%
  dplyr::group_by(series) %>%
  dplyr::mutate(nb_cells = n())
meta_cells$sample_feature <- paste0(meta_cells$timepoint,' ',meta_cells$tm_st, ': ',meta_cells$nb_cells,' cells')
meta_cells$series <- factor(meta_cells$series, levels = gtools::mixedsort(unique(meta_cells$series)))
meta_cells$sample_feature <- factor(meta_cells$sample_feature, levels = gtools::mixedsort(unique(meta_cells$sample_feature)))
meta_cells_535_cx <- meta_cells
p <- plot_variation_function_v3(meta_data=meta_cells, xstring="series", ystring="total_features_by_counts",
                                plottype="sample_feature", plottitle=paste0('SA535:CX5461', suffix1),
                                xlabel=NULL, ylabel=NULL,lg_pos="right")


p_c <- plot_variation_function_v3(meta_data=meta_cells, xstring="series", ystring="total_counts",
                                  plottype="sample_feature", plottitle=paste0('SA535:CX5461', suffix2),
                                  xlabel=NULL, ylabel=NULL,lg_pos='none')
# p_c
p535_cx_raw <- cowplot::plot_grid(p, p_c, ncol=2, align = 'h', rel_widths = c(4.2,3))
datatag <- gsub(' ','_',datatag)
saveRDS(p535_cx_raw, paste0(output_dir,datatag,"_raw_data_features_v2.rds"))
png(paste0(output_dir,datatag,"_raw_data_features.png"), height = 2*420, width=2*860,res = 2*72)
print(p535_cx_raw)
dev.off()

res_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/'
p_raw <- cowplot::plot_grid(p609_raw, p1035_raw, p535_cis_raw, p535_cx_raw, nrow = 4, 
                            align='v', labels = c('a','b','c','d'))

png(paste0(res_dir,"diagnostic_plot_3series_raw_data_features.png"), height = 2*1150, width=2*1100,res = 2*72)
print(p_raw)
dev.off()

ggsave(paste0(res_dir,"diagnostic_plot_3series_raw_data_features_23_Feb.png"),
       plot = p_raw,
       height = 8.5,
       width = 6.5,
       # useDingbats=F,
       type = "cairo-png"
)
ggsave(paste0(res_dir,"diagnostic_plot_3series_raw_data_features_23_Feb.pdf"),
       plot = p_raw,
       height = 8.5,
       width = 6.5,
       useDingbats=F
)#dpi = 150


meta_cells_609
stat_609 <- get_rawdata_stat(meta_cells_609, "SA609")

meta_cells_1035
stat_1035 <- get_rawdata_stat(meta_cells_1035, "SA1035")

meta_cells_535_cis
stat_535_cis <- get_rawdata_stat(meta_cells_535_cis, "SA535:Cisplatin")

meta_cells_535_cx
stat_535_cx <- get_rawdata_stat(meta_cells_535_cx, "SA535:CX5461")

stat_609$desc
stat_1035$desc
stat_535_cis$desc
stat_535_cx$desc

get_rawdata_stat <- function(df, datatag){
  df <- df %>%
    dplyr::group_by(series) %>%
    dplyr::summarise(avg_tc=round(mean(total_counts),2),
                     avg_fc=round(mean(total_features_by_counts),2),
                     med_tc=round(median(total_counts),2),
                     med_fc=round(median(total_features_by_counts),2),
                     iqr_tc=IQR(total_counts),
                     iqr_fc=IQR(total_features_by_counts))
  
  df2 <- df %>%
    dplyr::summarise(min_med_tc=min(avg_tc),
                     max_med_tc=max(avg_tc),
                     med_med_tc=round(median(avg_tc),1),
                     min_med_fc=min(med_fc),
                     max_med_fc=max(avg_tc),
                     med_med_fc=round(median(avg_fc),1),
                     min_iqr_tc=round(min(iqr_tc),1),
                     max_iqr_tc=round(max(iqr_tc),1),
                     min_iqr_fc=round(min(med_fc),1),
                     max_iqr_fc=round(max(med_fc),1))
  desc <- paste0("Series ",datatag, ': total counts - median over all samples: ',df2$med_med_tc,
                 ' interquartile range: [',df2$min_iqr_tc,'-',df2$max_iqr_tc,']  ',
                 ' number of genes with non-zero counts - median over all samples: ',df2$med_med_fc,
                 ' interquartile range: [',df2$min_iqr_fc,'-',df2$max_iqr_fc,']  ')
  return(list(datatag=datatag,
              desc=desc, 
              df_series=df,
              df_summary=df2))
}  


