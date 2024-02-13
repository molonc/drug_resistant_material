source("/home/htran/Projects/farhia_project/drug_resistant_material/scripts/cis_trans/in_cis_trans_utils.R")
source("/home/htran/Projects/farhia_project/drug_resistant_material/scripts/cis_trans/viz_utils.R")

data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
datatag <- 'SA604'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
res <- load_data(pair_groups_fn, datatag, cnv_fn)
dim(res$pair_groups)
res <- get_gene_type_edgeR_v2(res$pair_groups, res$cnv_mat, 
                              datatag, input_dir, save_dir,
                              FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)


data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
pair_groups_fn <- paste0(input_dir,'comparisons_across_treatments_revision35.csv')
save_dir <- paste0(data_dir,'signif_genes_revision/')
datatag <- 'SA609'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
res <- load_data(pair_groups_fn, datatag, cnv_fn)
dim(res$pair_groups)
# head(res$cnv_mat)
res <- get_gene_type_edgeR_v2(res$pair_groups, res$cnv_mat, 
                              datatag, input_dir, save_dir,
                              FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)


data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
pair_groups_fn <- paste0(input_dir,'comparisons_across_treatments_revision35.csv')
datatag <- 'SA535'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
save_dir <- paste0(data_dir,'signif_genes_revision/')
res <- load_data(pair_groups_fn, datatag, cnv_fn)
dim(res$pair_groups)
res <- get_gene_type_edgeR_v2(res$pair_groups, res$cnv_mat, 
                              datatag, input_dir, save_dir,
                              FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)


data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
pair_groups_fn <- paste0(input_dir,'comparisons_across_treatments_revision35.csv')
datatag <- 'SA1035'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
save_dir <- paste0(data_dir,'signif_genes_revision/')

res <- load_data(pair_groups_fn, datatag, cnv_fn)
# View(res$pair_groups)
pair_groups <- res$pair_groups
dim(pair_groups)
# pair_groups <- pair_groups %>%
#   dplyr::filter(file_header=='scrande_SA1035_52')
# View(pair_groups)
summary(as.factor(pair_groups$comp_type))
pair_groups <- pair_groups %>% 
  dplyr::filter(clone1!='' & clone2 !='')
res <- get_gene_type_edgeR_v2(pair_groups, res$cnv_mat, 
                              datatag, input_dir, save_dir,
                              FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)


## Summary cis/trans genes proportion in whole dataset
library(dplyr)
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
pair_groups_fn <- paste0(input_dir,'comparisons_across_treatments_revision35.csv')
save_dir <- paste0(data_dir,'signif_genes_revision/')
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
# pair_groups <- pair_groups[pair_groups$order!=0,] 
dim(pair_groups)
# unique(pair_groups$datatag)
# View(pair_groups)
stat <- get_gene_type_stat_v3(pair_groups, save_dir, save_dir)
dim(stat)

stat <- data.table::fread(paste0(save_dir,'wholedata_summary.csv')) %>% as.data.frame()
dim(stat)

head(stat)
t <- stat %>%
  dplyr::filter(file_header=="scrande_SA609_202")  %>%
  dplyr::select(gene_type, pct_genes)
                  
sum(t$pct_genes)
## Summary results
# 3.5 part 1: Overall, what are the percentages of copy number in-cis transcripts 
# in treated groups comparing to un-treated groups?
unique(stat$comp_type)
unique(stat$datatag)
stat1 <- stat %>%
   dplyr::filter(gene_type!="unmapped" & 
                comp_type %in% c("UnRx_vs_UnRx","Rx_expand_vs_UnRx","RxH_vs_UnRx")) %>%
   dplyr::mutate(
     group=case_when(
       comp_type=="UnRx_vs_UnRx" ~ "untreated groups",
       TRUE ~ "treated groups"),
     gt=case_when(
       grepl("In_cis",gene_type) ~ 'cis',
       TRUE ~ "trans"
     ))  
unique(stat$comp_type)
stat1 <- stat %>%
  dplyr::filter(gene_type!="unmapped") %>% #  & comp_type !="Rx_expand_vs_Rx_shrink"
  dplyr::mutate(
    gt=case_when(
      grepl("In_cis",gene_type) ~ 'cis gene',
      TRUE ~ "trans gene"
    ))  
dim(stat1)

stat11 <- stat1 %>%
  dplyr::group_by(gt, labels_detailed, comp_type, series) %>%
  dplyr::summarise(total_pct_genes=sum(pct_genes))
dim(stat11)
stat11$comp_type <- gsub('_',' ',stat11$comp_type)
library(ggplot2)
p <- ggplot(stat11, aes(x=comp_type, y=total_pct_genes, shape=comp_type)) + 
  geom_jitter(position=position_jitter(0.2), cex=2.5) + 
  facet_wrap(gt~series,nrow = 2,scales = "free_y") + 
  labs(x=NULL,y='% gene type',title='Gene dosage effects in treated, drug holiday, untreated groups')
p <- p + stat_summary(fun.y=median, geom="point", shape=18,
                        size=2, color="red")
my_font <- "Helvetica"
thesis_theme <- ggplot2::theme(
  text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
  axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),  
  axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
  axis.text.x = element_text(color="black",size=10, hjust = 0.5, family=my_font, angle = 90),
  # axis.text.x = element_blank(),
  axis.text.y = element_text(color="black",size=7, hjust = 0.5, family=my_font),
  plot.title = element_text(color="black",size=10, face="bold", hjust=0, family=my_font),
  legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font), 
  legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
  strip.text.x = element_text(color="black",size=10, family=my_font),
  strip.text.y = element_text(color="black",size=10, family=my_font),
  legend.spacing.x = unit(0.1, 'mm'),
  legend.spacing.y = unit(0.1, 'mm'),
  legend.key.height=unit(1,"line"),
  legend.position = 'none',
  # panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                  colour = "lightgrey"),
  panel.background = element_rect(fill = "white", #,
                                  # size = 0.5, linetype = "solid"
                                  colour = "lightgrey")
)
p <- p + thesis_theme
p

png(paste0(save_dir,'revision_question3.5.png'), height = 2*550, width=2*680,res = 2*72)
print(p)
dev.off()

ggsave(paste0(save_dir,'revision_question3.5.svg'),
       plot = p,
       height = 5.5,
       width = 7,
       # useDingbats=F,
       dpi=250)


png(paste0(save_dir,'revision_question3.5_2.png'), height = 2*600, width=2*700,res = 2*72)
print(p)
dev.off()

ggsave(paste0(save_dir,'revision_question3.5_2.svg'),
       plot = p,
       height = 6,
       width = 7,
       # useDingbats=F,
       dpi=250)

dim(stat11)  

stat12 <- stat11 %>%
  dplyr::filter(comp_type !="Rx_expand_vs_Rx_shrink") %>%
  dplyr::group_by(gt, comp_type, series) %>%
  dplyr::summarise(avg_pct_genes=round(mean(total_pct_genes),2),
                   sd_pct_genes=round(sd(total_pct_genes),2))
stat12 <- as.data.frame(stat12)
rownames(stat12) <- paste0(stat12$series, '_',stat12$comp_type,stat12$gt)
for(s in unique(stat12$series)){
  tmp <- stat12 %>%
    dplyr::filter(series==s)
  for(c in rownames(tmp)){
    print(paste0(tmp[c,'series'],'; ',tmp[c,'comp_type'],' ',
                 tmp[c,'gt'],': ',tmp[c,'avg_pct_genes'],'+-',tmp[c,'sd_pct_genes']))
  }
}
dim(stat12)
# stat13 <- stat12 %>%
#   tidyr::pivot_wider(names_from = 'gt', values_from = 'avg_pct_genes')
# View(stat13)


data.table::fwrite(stat12, paste0(input_dir,'revision_35_summary.csv'))



stat12 <- stat11 %>%
  dplyr::group_by(datatag, gt, group) %>%
  dplyr::summarise(total_pct_genes=total_pct_genes)
stat12

data.table::fwrite(stat1, paste0(input_dir,'revision_35_part1_summary.csv'))

## Summary results
# 3.5 part 2:  Would dosage effect stronger or weaker 
# in expanded clones compared to other shrinking clones?
unique(stat$comp_type)
stat2 <- stat %>%
  dplyr::filter(gene_type!="unmapped" & 
                  comp_type %in% c("Rx_expand_vs_UnRx","Rx_shrink_vs_UnRx")) %>%
  dplyr::mutate(
    gt=case_when(
      grepl("In_cis",gene_type) ~ 'cis',
      TRUE ~ "trans"
    ))  
dim(stat2) 
unique(stat2$labels_detailed)
stat21 <- stat2 %>%
  dplyr::group_by(datatag, gt, comp_type) %>%
  dplyr::summarise(total_pct_genes=mean(pct_genes))
dim(stat21)

stat21
data.table::fwrite(stat2, paste0(input_dir,'revision_35_part2_summary.csv'))


stat23 <- stat2 %>%
  dplyr::group_by(datatag, gt, labels_detailed) %>%
  dplyr::summarise(total_pct_genes=mean(pct_genes))
dim(stat23)
View(stat23)


# checking
df <- stat2 %>%
  dplyr::filter(datatag=="SA535" & gt=="cis") %>%
  dplyr::group_by(labels_detailed, comp_type, passage) %>%
  dplyr::summarise(total_pct_genes=sum(pct_genes))
dim(df)
View(df)
df
