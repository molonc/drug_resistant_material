

# source("/home/htran/Projects/farhia_project/rnaseq/cis_trans/in_cis_trans_utils.R")
# source("/home/htran/Projects/farhia_project/rnaseq/cis_trans/viz_utils.R")

source("/home/htran/Projects/farhia_project/drug_resistant_material/scripts/cis_trans/in_cis_trans_utils.R")
source("/home/htran/Projects/farhia_project/drug_resistant_material/scripts/cis_trans/viz_utils.R")

      

data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
save_dir <- paste0(data_dir,'signif_genes/')
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
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
datatag <- 'SA535'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
save_dir <- paste0(data_dir,'signif_genes/')
res <- load_data(pair_groups_fn, datatag, cnv_fn)
dim(res$pair_groups)
res <- get_gene_type_edgeR_v2(res$pair_groups, res$cnv_mat, 
                       datatag, input_dir, save_dir,
                       FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)


data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
datatag <- 'SA1035'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
save_dir <- paste0(data_dir,'signif_genes/')

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


data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
save_dir <- paste0(data_dir,'signif_genes/')
datatag <- 'SA501'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
res <- load_data(pair_groups_fn, datatag, cnv_fn)
# View(res$pair_groups)
pair_groups <- res$pair_groups
dim(pair_groups)
colnames(pair_groups)
res <- get_gene_type_edgeR_v2(pair_groups, res$cnv_mat, 
                       datatag, input_dir, save_dir,
                       FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)

datatag <- 'SA530'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
res <- load_data(pair_groups_fn, datatag, cnv_fn)
dim(res$pair_groups)
res <- get_gene_type_edgeR_v2(res$pair_groups, res$cnv_mat, 
                       datatag, input_dir, save_dir,
                       FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)


datatag <- 'SA604'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
res <- load_data(pair_groups_fn, datatag, cnv_fn)
dim(res$pair_groups)
res <- get_gene_type_edgeR_v2(res$pair_groups, res$cnv_mat, 
                       datatag, input_dir, save_dir,
                       FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)


## Summary cis/trans genes proportion in whole dataset
data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
input_dir <- paste0(data_dir,'signif_genes/')
save_dir <- paste0(data_dir,'signif_genes/')
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
# pair_groups <- pair_groups[pair_groups$order!=0,] 
dim(pair_groups)
# unique(pair_groups$datatag)
# View(pair_groups)
stat <- get_gene_type_stat_v3(pair_groups, input_dir, save_dir)
data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
save_dir <- paste0(data_dir,'signif_genes/')
stat <- data.table::fread(paste0(save_dir,'wholedata_summary.csv')) %>% as.data.frame()
dim(stat)
# head(stat)
# View(pair_groups)


## Early time where we can not get good mapping between DLP+ clones and RNAseq
# excluded_early_times <- c('scrande_SA535_4','scrande_SA535_5','scrande_SA535_52')
# stat <- stat %>%
#   dplyr::filter(!file_header %in% excluded_early_times)
# dim(stat)

# unique(stat$comp_type)
# unique(stat$datatag)

## Plotting cistrans proportion for treated comparisons
stat1 <- stat
# stat <- stat1
comp_treated <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_treated.csv')) %>% as.data.frame()
# comp_untreated_DH <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_untreated_DH.csv')) %>% as.data.frame()
comp_treated <- comp_treated %>%
  dplyr::filter(file_header!='scrande_SA1035_41') %>% # early time and do not contain cis genes. 
  dplyr::select(file_header, order)

stat <- stat %>%
  dplyr::filter(file_header %in% comp_treated$file_header) %>%
  dplyr::select(-order)%>% 
  left_join(comp_treated, by='file_header')
stat1 <- stat
# plot_ls <- plot_intrans_incis_prevalence_v3(stat, verbose=F)
plot_ls <- plot_intrans_incis_prevalence_v4(stat, save_dir, verbose=F)

plot_ls$treated


comp <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_untreated_DH.csv')) %>% as.data.frame()
dim(comp)
summary(as.factor(comp$comp_type))
table(comp$comp_type, comp$datatag)
comp_RxH <- comp %>%
  dplyr::filter(comp_type=='drugHoliday_vs_untreated') %>% 
  dplyr::select(file_header, order)
dim(comp_RxH)
# View(comp_RxH)
stat <- data.table::fread(paste0(save_dir,'wholedata_summary.csv')) %>% as.data.frame()
dim(stat)
stat <- stat %>%
  dplyr::filter(file_header %in% comp_RxH$file_header) %>%
  dplyr::select(-order)%>% 
  left_join(comp_RxH, by='file_header')
dim(stat)
plot_ls <- plot_intrans_incis_prevalence_v4(stat, save_dir, verbose=F)

plot_ls$treated
plot_ls$treated_proportion
plot_ls$treated_cistrans
plot_ls$treated_SA609_c1lb
plot_ls$treated_SA609_passage


## Plotting untreated stat values
stat <- data.table::fread(paste0(save_dir,'wholedata_summary.csv')) %>% as.data.frame()
dim(stat)
comp_UnRx <- comp %>%
  dplyr::filter(comp_type=='untreated') 
comp_UnRx$passage <- c('X7','X6','X3','X6','X9','X8')
comp_UnRx <- comp_UnRx %>%
  dplyr::select(file_header, order, passage)
dim(comp_UnRx)
# View(comp_UnRx)
stat <- stat %>%
  dplyr::filter(file_header %in% comp_UnRx$file_header) %>%
  dplyr::select(-order, -passage)%>% 
  left_join(comp_UnRx, by='file_header')
dim(stat)
summary(as.factor(stat$passage))

plot_untreated_ls <- plot_intrans_incis_prevalence_untreated_v4(stat, save_dir, verbose=F)
plot_untreated_ls$UnRx
plot_untreated_ls$UnRx_proportion
plot_untreated_ls$UnRx_cistrans
plot_untreated_ls$UnRx_SA501_c2lb
plot_untreated_ls$UnRx_SA501_passage


# plot_ls$untreated_proportion
# plot_ls$SA609
# plot_ls$SA609_proportion
# plot_ls$SA535
# plot_ls$SA535_proportion
# plot_ls$SA1035
# plot_ls$SA1035_proportion
# dev.off()
# plot_ls$untreated_cistrans
# plot_ls$SA1035_cistrans
# plot_ls$untreated_c1lb
# plot_ls$untreated_c2lb
# plot_ls$untreated_passage
# plot_ls$SA609_c1lb
# plot_ls$SA609_c2lb
# plot_ls$SA609_passage
# plot_ls$SA609

## Get gene type that belong to reference gene sets
data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
input_dir <- paste0(data_dir,'signif_genes/')
save_dir <- input_dir
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
# pair_groups <- pair_groups[pair_groups$order!=0,] 
dim(pair_groups)
unique(pair_groups$comp_type)
pair_groups <- pair_groups %>%
  dplyr::filter(comp_type !='untreated')
res <- get_proportion_cistrans_reference_set(pair_groups, input_dir, save_dir)

# dim(res$stat_cf)
# dim(res$stat_cr)
# View(res$stat_cf)
## Pathway analysis
data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
input_dir <- paste0(data_dir,'signif_genes/')
save_dir <- input_dir
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
# pair_groups <- pair_groups[pair_groups$order!=0,] 
pair_groups <- pair_groups %>%
  dplyr::filter(comp_type !='untreated')

pathway_stat <- get_pathways_reference_set(pair_groups, input_dir, save_dir)
dim(pathway_stat)
sum(pathway_stat$padj<0.05)
pathway_stat1 <- pathway_stat %>%
  dplyr::filter(padj<0.05)

## Plotting here 
data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
input_dir <- paste0(data_dir,'signif_genes/')
save_dir <- input_dir
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
pair_groups <- pair_groups[pair_groups$order!=0,] 
pair_groups <- pair_groups %>%
  dplyr::filter(comp_type !='untreated')
stat_cf <- data.table::fread(paste0(input_dir, 'wholedata_summary_CF.csv')) %>% as.data.frame()
stat_cr <- data.table::fread(paste0(input_dir, 'wholedata_summary_CR.csv')) %>% as.data.frame()
stat <- dplyr::bind_rows(stat_cf, stat_cr)
data.table::fwrite(stat, paste0(input_dir, 'wholedata_CF_CR_sets.csv'))
rm(stat_cf)
rm(stat_cr)

stat_cfcr <- data.table::fread(paste0(input_dir,'wholedata_CF_CR_sets.csv'), header=T) %>% as.data.frame()
dim(stat_cfcr)
# unique(stat$ref_set)
# excluded_early_times <- c('scrande_SA535_4','scrande_SA535_5','scrande_SA535_52')
# stat <- stat %>%
#   dplyr::filter(!file_header %in% excluded_early_times)
# dim(stat)
# unique(stat$ref_set)
## TO DO: summary this plot
# prop_plt_ls <- plot_reference_set_proportion_total(stat, save_dir, verbose=F)
unique(stat$ref_set)
comp_treated <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_treated.csv')) %>% as.data.frame()

comp_treated <- comp_treated %>%
  dplyr::filter(file_header!='scrande_SA1035_41') %>% # early time and do not contain cis genes. 
  dplyr::select(file_header, order)

stat_cfcr <- stat_cfcr %>%
  dplyr::filter(file_header %in% comp_treated$file_header) %>%
  dplyr::select(-order)%>% 
  left_join(comp_treated, by='file_header')
stat$ref_set
prop_plt_ls <- plot_reference_set_proportion_v2(stat_cfcr, save_dir, verbose=verbose, plttype = 'treated')
prop_plt_ls$prevalence_plt
# prop_plt_ls$SA609
# prop_plt_ls$SA535
# prop_plt_ls$SA1035
# prop_plt_ls$plg
# dev.off()


comp <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_untreated_DH.csv')) %>% as.data.frame()
dim(comp)
summary(as.factor(comp$comp_type))
table(comp$comp_type, comp$datatag)
comp_RxH <- comp %>%
  dplyr::filter(comp_type=='drugHoliday_vs_untreated') %>% 
  dplyr::select(file_header, order) #, labels_detailed
dim(comp_RxH)
head(comp_RxH)
stat_cfcr <- stat_cfcr %>%
  dplyr::filter(file_header %in% comp_RxH$file_header) %>%
  dplyr::select(-order)%>% 
  left_join(comp_RxH, by='file_header')
dim(stat_cfcr)
colnames(stat_cfcr)
prop_plt_ls <- plot_reference_set_proportion_v2(stat_cfcr, save_dir, verbose=F, plttype = 'RxH')
prop_plt_ls$prevalence_plt


data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
input_dir <- paste0(data_dir,'signif_genes/')
save_dir <- input_dir
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
# pair_groups <- pair_groups[pair_groups$order!=0,] 

pw_stat <- data.table::fread(paste0(input_dir, 'wholedata_pathways_CF_CR_sets.csv')) %>% as.data.frame()
dim(pw_stat)

comp_treated <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_treated.csv')) %>% as.data.frame()
comp_treated <- comp_treated %>%
  dplyr::filter(file_header!='scrande_SA1035_41') %>% # early time and do not contain cis genes. 
  dplyr::select(file_header, order, series)

pw_stat <- pw_stat %>%
  dplyr::filter(file_header %in% comp_treated$file_header) %>%
  dplyr::filter(padj<0.05)%>%
  # dplyr::select(-order)%>% 
  left_join(comp_treated, by='file_header')
dim(pw_stat)
# excluded_early_times <- c('scrande_SA535_4','scrande_SA535_5','scrande_SA535_52')
# pw_stat <- pw_stat %>%
#   dplyr::filter(!file_header %in% excluded_early_times)
View(pw_stat)

unique(pw_stat$series)
# pathway_plt_ls <- plot_pathway_custom_sets_total(pw_stat, pair_groups, save_dir, verbose=F)
pathway_plt_ls <- viz_reference_pathways(pw_stat, verbose=F, plttype = 'treated')
pathway_plt_ls$pathway_sig_plt
pathway_plt_ls$plg


## For drug holiday
pw_stat <- data.table::fread(paste0(input_dir, 'wholedata_pathways_CF_CR_sets.csv')) %>% as.data.frame()
dim(pw_stat)

comp_RxH <- comp %>%
  dplyr::filter(comp_type=='drugHoliday_vs_untreated') %>% 
  dplyr::select(file_header, order, series) #, labels_detailed
dim(comp_RxH)

pw_stat <- pw_stat %>%
  dplyr::filter(file_header %in% comp_RxH$file_header) %>%
  # dplyr::select(-order)%>% 
  dplyr::filter(padj<0.05)%>%
  left_join(comp_RxH, by='file_header')
# excluded_early_times <- c('scrande_SA535_4','scrande_SA535_5','scrande_SA535_52')
# pw_stat <- pw_stat %>%
#   dplyr::filter(!file_header %in% excluded_early_times)
# dim(pw_stat)
dim(pw_stat)
pw_stat$pathway
unique(pw_stat$series)
# View(pw_stat)

# pathway_plt_ls <- plot_pathway_custom_sets_total(pw_stat, pair_groups, save_dir, verbose=F)
pathway_plt_ls <- viz_reference_pathways(pw_stat, verbose=F, plttype = 'RxH')
pathway_plt_ls$pathway_sig_plt
pathway_plt_ls$plg


p_total <- viz_Fig4_v2(plot_ls, prop_plt_ls, pathway_plt_ls, save_dir)
save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/cis_trans/SUPPFig6_UnRx_RxH_cistrans/'
p_total <- viz_SUPP_Fig6_DrugHolidayRxH(plot_ls, plot_untreated_ls, prop_plt_ls, NULL, save_dir)
# p_total <- viz_Fig4(plot_ls, prop_plt_ls, pathway_plt_ls, save_dir)
# p_total
dev.off()
# View(pair_groups)
# ## Compare old and new output
# fn1 <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/SA609/SA609_3_SA609_UTTT_R_UUUU_H/signif_genes_FDR_0.01.csv'
# s1 <- data.table::fread(fn1) %>% as.data.frame()
# dim(s1)
# 
# fn2 <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/signif_genes/scrande_SA609_3/signif_genes.csv'
# s2 <- data.table::fread(fn2) %>% as.data.frame() %>%
#   dplyr::filter(abs(logFC)>0.5 & FDR<0.01 & PValue<0.05)
# dim(s2)
# 
# df <- table(s1$Gene_Type, s2$Gene_Type)%>% as.data.frame()
# df <- df[df$Freq>5,]
# View(df)
# ([0-9]+|X|Y)_[0-9]+_[0-9]+.   ([a-z])+_([AZ][0-9])+_[0-9]+



obs_comps <- c('rstmm1_SA609_UT_R_UU_H','rstmm2_SA609_UTT_R_UUU_H','rstmm3_SA609_UTTT_R_UUUU_H','rstmm4_SA609_UTTTT_R_UUUUU_H','fitness31_SA609_UUUU_H_UUUU_C')
for(comp in obs_comps){
  de_genes <- data.table::fread(paste0(save_dir, comp,'/signif_genes_FDR0.01.csv')) %>% as.data.frame()
  # dim(de_genes)
  # summary(as.factor(de_genes$Gene_Type))
  de_genes <- de_genes %>%
    dplyr::filter(abs(logFC)>=0.5)
  cis_pos <- de_genes %>%
    dplyr::filter(Gene_Type %in% c('In_cis_Decrease_DownRegulated','In_cis_Increase_UpRegulated'))%>%
    dplyr::pull(ensembl_gene_id)
  cis_neg <- de_genes %>%
    dplyr::filter(Gene_Type %in% c('In_cis_Decrease_UpRegulated','In_cis_Increase_DownRegulated'))%>%
    dplyr::pull(ensembl_gene_id)
  nb_cis <- length(cis_pos) + length(cis_neg)
  print(paste0(comp,' :   pos(DD,IU): ',round(length(cis_pos)*100/nb_cis,1),'  neg(DU, ID): ',round(length(cis_neg)*100/nb_cis,1)))  
  
}


## Summary cis/trans genes proportion in whole dataset, put values into manuscript
data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
input_dir <- paste0(data_dir,'signif_genes/')
save_dir <- paste0(data_dir,'signif_genes/')
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
pair_groups <- pair_groups[pair_groups$order!=0,] 
# dim(pair_groups)
# unique(pair_groups$datatag)
# View(pair_groups)

excluded_early_times <- c('scrande_SA535_4','scrande_SA535_5','scrande_SA535_52')
pair_groups <- pair_groups %>%
  dplyr::filter(!file_header %in% excluded_early_times)

stat <- get_gene_type_stat_manuscript(pair_groups, input_dir, save_dir)
# stat <- data.table::fread(paste0(save_dir,'wholedata_summary_manuscript.csv.gz')) %>% as.data.frame()
dim(stat)
stat <- stat %>%
  dplyr::filter(gene_type!='unmapped')

head(stat)
unique(stat$gene_type)

stat_treated <- stat %>%
  dplyr::filter(datatag %in% c('SA609','SA535','SA1035') &
                  comp_type=='treated_vs_untreated' &
                  gene_type %in% c('cis_neg'))%>% #,'cis_pos'
  dplyr::group_by(datatag)%>%
  dplyr::summarise(avgCis=round(mean(pct_genes),2),
                   minCis=round(min(pct_genes),2),
                   maxCis=round(max(pct_genes),2),
                   sdLogFC=round(sd(pct_genes),2))
stat_treated <- stat %>%
  dplyr::filter(datatag %in% c('SA609','SA535','SA1035') &
                  comp_type=='treated_vs_untreated' &
                  gene_type %in% c('cis'))%>% #,'cis_pos'
  dplyr::group_by(datatag)%>%
  dplyr::summarise(avgCis=round(mean(pct_genes),2),
                   minCis=round(min(pct_genes),2),
                   maxCis=round(max(pct_genes),2),
                   sdLogFC=round(sd(pct_genes),2))
stat$file_header
s1 <- stat %>%
  dplyr::filter(datatag %in% c('SA609','SA535','SA1035') &
                comp_type=='treated_vs_untreated' &
                gene_type %in% c('cis_neg','cis_pos'))%>% #,'cis_pos'
  dplyr::group_by(datatag, file_header)%>%
  dplyr::summarise(nb_totalcis=sum(nb_genes))

s1 <- stat %>%
  dplyr::filter(datatag %in% c('SA609','SA535','SA1035') &
                  comp_type=='treated_vs_untreated' &
                  gene_type %in% c('cis_neg'))%>% #,'cis_pos'
  dplyr::group_by(datatag)%>%
  dplyr::summarise(var_pct=sd(pct_genes))

s1

s1 <- stat %>%
  dplyr::filter(datatag %in% c('SA609','SA535','SA1035') & #
                  comp_type=='treated_vs_untreated' &
                  gene_type %in% c('In_cis_Increase_UpRegulated','In_trans_UpRegulated'))%>% #,'cis_pos'
  dplyr::group_by(datatag, gene_type)%>%
  dplyr::summarise(var_pct=sd(pct_genes))

s2 <- s1 %>%
  dplyr::group_by(gene_type)%>%
  dplyr::summarise(var_var_pct=sd(var_pct))
s2

s1 <- stat %>%
  dplyr::filter(datatag %in% c('SA609','SA535','SA1035') & #
                  comp_type=='treated_vs_untreated' &
                  gene_type %in% c('In_cis_Decrease_DownRegulated','In_trans_DownRegulated'))%>% #,'cis_pos'
  dplyr::group_by(datatag, gene_type)%>%
  dplyr::summarise(var_pct=sd(pct_genes))

s2 <- s1 %>%
  dplyr::group_by(gene_type)%>%
  dplyr::summarise(var_var_pct=sd(var_pct))
s2

dim(s1)                   
s1                   
s1 <- stat %>%
  dplyr::filter(datatag %in% c('SA609') & #,'SA535','SA1035'
                  comp_type=='treated_vs_untreated' &
                  gene_type %in% c('trans'))%>% #,'cis_pos'
  dplyr::group_by(datatag)%>%
  dplyr::summarise(var_pct=sd(pct_genes),
                   rg=max(pct_genes)-min(pct_genes))
dim(s1)


s2 <- stat %>%
  dplyr::filter(datatag %in% c('SA609','SA535','SA1035') &
                  comp_type=='treated_vs_untreated' &
                  gene_type %in% c('cis_pos'))%>% #,'cis_neg',
  dplyr::left_join(s1, by=c('file_header','datatag'))%>%
  dplyr::group_by(datatag, file_header)%>%
  dplyr::summarise(pct_genes=round(nb_genes/nb_totalcis*100,2))%>%
  # dplyr::group_by(datatag)%>%
  dplyr::ungroup()%>%
  dplyr::summarise(avg_pct_genes=round(mean(pct_genes),2),
                   sd_pct_genes=round(sd(pct_genes),2))
s2 <- stat %>%
  dplyr::filter(datatag %in% c('SA609','SA535','SA1035') &
                  comp_type=='treated_vs_untreated' &
                  gene_type %in% c('cis_neg'))%>% #,'cis_neg',
  dplyr::left_join(s1, by=c('file_header','datatag'))%>%
  dplyr::group_by(datatag, file_header)%>%
  dplyr::summarise(pct_genes=round(nb_genes/nb_totalcis*100,2))%>%
  # dplyr::group_by(datatag)%>%
  dplyr::ungroup()%>%
  dplyr::summarise(avg_pct_genes=round(mean(pct_genes),2),
                   sd_pct_genes=round(sd(pct_genes),2))

dim(s2)
s2 <- stat %>%
  dplyr::filter(datatag %in% c('SA535') &
                  comp_type=='treated_vs_untreated' &
                  gene_type %in% c('cis_pos'))%>% #,'cis_neg',
  dplyr::select(datatag,gene_type, order, nb_genes, pct_genes)%>%
  dplyr::arrange(order) 
  
s2
View(s2)
dim(stat)
colnames(s2)
s2$nb_totalcis
stat_treated <- stat %>%
  dplyr::filter(datatag %in% c('SA609','SA535','SA1035'))%>% # & gene_type=='cis_pos'
  dplyr::filter(gene_type %in% c('cis_neg','cis_pos'))%>% #,'cis_pos'
  tidyr::pivot_wider(names_from = 'gene_type', values_from = 'nb_genes')
  
dim(stat_treated)
View(head(stat_treated))






data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/SA535_v2/'
pair_groups_fn <- paste0(data_dir,'tested_SA535_early/tested_SA535_early.csv')
datatag <- 'SA535'
cnv_fn <- paste0(dlp_dir,datatag,'_cnv_mat.csv.gz')
save_dir <- paste0(data_dir,'tested_SA535_early/')
# dir.create(save_dir)


res <- load_data(pair_groups_fn, datatag, cnv_fn)
# View(res$pair_groups)
df <- get_gene_type_edgeR_v2(res$pair_groups, res$cnv_mat, 
                       datatag, input_dir, save_dir,
                       FDR_cutoff=0.05, minLogFC=0.25, pValueThrs=0.05)




data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
input_dir <- paste0(data_dir,'signif_genes/')
# View(pair_groups)
pair_groups <- pair_groups %>%
  dplyr::filter(comp_type=="treated_vs_untreated" & order!=0)
save_dir <- paste0(input_dir,'ks_cistrans_test/')
res <- get_stat_logFC(pair_groups, input_dir, save_dir)
dim(res$combined_DE)
# length(unique(res$combined_DE$file_header))
# tag <- 'SA535'
# out_SA535 <- out %>%
#   dplyr::filter(datatag==tag) %>%
#   dplyr::arrange(order) %>%
#   dplyr::select(datatag, labels_detailed, ks_stat, ks_signf)

stat_df <- res$pair_groups
de_df <- res$combined_DE


stat_df <- data.table::fread(paste0(save_dir, 'pair_groups_ks_test.csv.gz'))  
de_df <- data.table::fread(paste0(save_dir, 'combined_DE_logFC.csv.gz'))  
datatag <- 'all'
p1 <- get_boxplot_ksStat(stat_df, de_df, save_dir, datatag='allSeries', plottype='boxplot')
p11 <- get_boxplot_ksStat(stat_df, de_df, save_dir, datatag='allSeries', plottype='violinplot')
saveRDS(p1, paste0(save_dir,datatag,'_ks_test_boxplot.rds'))
saveRDS(p11, paste0(save_dir,datatag,'_ks_test_violinplot.rds'))


datatag='all'
pdis_logFC <- readRDS(paste0(save_dir,datatag,'_ks_test_boxplot.rds'))
