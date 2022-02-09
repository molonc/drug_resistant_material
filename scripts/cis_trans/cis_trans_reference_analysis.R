library(dplyr)
library(ggplot2)

# TO DO 
## Top 30 up, and 10 down genes check first
## Chr distribution
## Pathway analysis: with current genes first
## Then go back to cis trans computation with all cis genes, check threshold used
## Hallmark, Kegg pathways
## Provide list of up, down-regulated genes at different passages
##
##
##
##

# ref_genes_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/reference_CNV_paper/'
# au_ref_genes <- data.table::fread(paste0(ref_genes_dir,'top30_AU_genes.csv'), header=T) %>% as.data.frame()
# length(au_ref_genes$Oncogenes)
# dd_ref_genes <- data.table::fread(paste0(ref_genes_dir,'top_10_DD_genes_.csv'), header=T) %>% as.data.frame()
# length(dd_ref_genes$Tumor_suppressor_genes)


input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
pair_df <- data.table::fread(pair_groups_fn) %>% as.data.frame()
rownames(pair_df) <- pair_df$desc

save_dir <- input_dir
# View(head(pair_df))
pair_df <- pair_df %>% 
  dplyr::filter(comp_type=='treated_vs_untreated')#%>% 
  # dplyr::select(-result_fn)
dim(pair_df)
rownames(pair_df) <- pair_df$desc
pair_df$datatag
xresults_dir <-  "/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/pathways_IU_DD/"

pathways_stat <- tibble::tibble()
obs_pw <- 'hallmark'
for(dtag in unique(pair_df$datatag)){
  print("------------------------")
  print(dtag)
  
  comp_ls <- pair_df %>% 
    dplyr::filter(datatag==dtag) 
  dim(comp_ls)
  if(dtag=='SA535'){
    comp_ls <- comp_ls %>% 
      dplyr::filter(plt_desc!="Res(X9:A) vs. Sen(X9:G)") 
  }
  for(de in comp_ls$desc){
    de_df <- data.table::fread(paste0(results_dir,de,'/signif_genes_FDR_0.05_minLogFC_0.05.csv')) %>% as.data.frame()
    de_df <- de_df %>%
      inner_join(meta_genes, by='ensembl_gene_id') %>%
      dplyr::filter(chr %in% chrs)  %>%
      dplyr::select(-gene_symbol) %>%
      dplyr::rename(gene_symbol=ref_gene_symbol)
    
    obs_au <- de_df %>%
      dplyr::filter(Gene_Type=='In_cis_Increase_UpRegulated')
      
    obs_dd <- de_df %>%
      dplyr::filter(Gene_Type=='In_cis_Decrease_DownRegulated')
    print("____IU___")
    print(dim(obs_au))
    gsea_pw_up <- get_custom_pathway_results(obs_au,                      # named vector of statistical significance
                                             desc=de, base_name = paste0(pair_df[de,'datatag'],'_IU'),
                                             pathway_name=obs_pw, #'cosmic' or 'cisplatin_resistance', or 'metastasis'
                                             save_dir = paste0(results_dir,de,'/'))
    gsea_pw_up$padj
    if(!is.null(gsea_pw_up)){
      gsea_pw_up <- gsea_pw_up %>%
        dplyr::filter(padj<0.3)
      if(!is.null(gsea_pw_up)){
        pathways_stat <- dplyr::bind_rows(pathways_stat, gsea_pw_up)
      }  
      print(gsea_pw_up$pathway)
    }
    
    # gsea_pw_up$signf_genes
    print("____DD___")
    print(dim(obs_dd))
    
    gsea_pw_down <- get_custom_pathway_results(obs_dd,                      # named vector of statistical significance
                                          desc=de, base_name = paste0(pair_df[de,'datatag'],'_DD'),
                                          pathway_name=obs_pw, #'cosmic' or 'cisplatin_resistance', or 'metastasis'
                                          save_dir = paste0(results_dir, de,'/'))
    
    if(!is.null(gsea_pw_down)){
      gsea_pw_down <- gsea_pw_down %>%
        dplyr::filter(padj<0.3)
      if(!is.null(gsea_pw_down)){
        pathways_stat <- dplyr::bind_rows(pathways_stat, gsea_pw_down)
      }  
      print(gsea_pw_down$pathway)
    }
  }
}  

View(pathways_stat)
data.table::fwrite(pathways_stat, paste0(results_dir,'pathways_summary.csv'))

get_common_genes <- function(){
  chrs <- c(as.character(1:22), "X")
  meta_genes <- annotables::grch38 %>%
    dplyr::select(ensembl_gene_id=ensgene, ref_gene_symbol=symbol, chr) 
  meta_genes <- meta_genes[!duplicated(meta_genes$ensembl_gene_id), ]
  # dim(meta_genes)
  results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
  
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
  pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
  pair_df <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  rownames(pair_df) <- pair_df$desc
  
  save_dir <- input_dir
  View(pair_df)
  pair_df <- pair_df %>% 
    dplyr::filter(desc!='SA535_2_SA535_UUTTTT_A_UUUUUU_G')
  pair_df <- pair_df %>% 
    dplyr::filter(comp_type=='treated_vs_untreated')%>% 
    dplyr::select(-result_fn)
  dim(pair_df)
  # pair_df_backup <- pair_df
  plt_desc_df <- data.table::fread(paste0(input_dir,'plt_desc_v2.csv')) %>% as.data.frame()
  dim(plt_desc_df)
  
  results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
  
  csv_dir <- paste0(results_dir, 'chr_distribution/')
  dim(pair_df)
  rownames(pair_df) <- pair_df$desc
  pair_df <- pair_df %>%
    dplyr::inner_join(plt_desc_df, by='desc') 
  unique(pair_df$comp_type)
  
  pair_Rx_UnRx <- pair_df %>% 
    dplyr::filter(comp_type=='treated_vs_untreated')%>% 
    dplyr::select(-result_fn)
  
  pair_UnRx <- pair_df %>% 
    dplyr::filter(comp_type=="untreated")%>% 
    dplyr::select(-result_fn)
  dim(pair_UnRx)
  dim(pair_Rx_UnRx)
  de_genes <- tibble::tibble()
  for(de in pair_df$desc){
    de_df <- data.table::fread(paste0(results_dir,pair_df[de,'datatag'],'/',de,'/signif_genes_FDR_0.01.csv')) %>% as.data.frame()
    # dim(de_df)
    de_df <- de_df %>%
      dplyr::filter(grepl('In_cis',Gene_Type)) %>%
      # inner_join(meta_genes, by='ensembl_gene_id') %>%
      # dplyr::filter(chr %in% chrs)  %>%
      # dplyr::group_by(Gene_Type)%>%
      dplyr::summarise(nb_genes=round(n()*100/dim(de_df)[1],1))%>%
      dplyr::mutate(datatag=pair_df[de,'datatag'], desc=de)
    # de_df <- de_df %>%
    #   dplyr::filter(grepl('In_trans',Gene_Type)) %>% 
    #   # inner_join(meta_genes, by='ensembl_gene_id') %>%
    #   # dplyr::filter(chr %in% chrs)  %>%
    #   # dplyr::group_by(Gene_Type)%>% 
    #   dplyr::summarise(nb_genes=round(n()*100/dim(de_df)[1],1))%>% 
    #   dplyr::mutate(datatag=pair_df[de,'datatag'], desc=de)
    
    de_genes <- dplyr::bind_rows(de_genes, de_df)
  }
  unique(de_df$Gene_Type)
  dim(de_genes)
  de_genes_UnRx <- de_genes %>%
    dplyr::inner_join(pair_UnRx, by='desc')
  de_genes_Rx_UnRx <- de_genes %>%
    dplyr::inner_join(pair_Rx_UnRx, by='desc')
  get_stat(de_genes_UnRx$nb_genes, 'in cis genes', 'Untreated series comparison')
  get_stat(de_genes_Rx_UnRx$nb_genes, 'in cis genes', 'Treated vs. untreated comparisons')  
    
  de_genes <- tibble::tibble()
  for(de in pair_df$desc){
    de_df <- data.table::fread(paste0(results_dir,pair_df[de,'datatag'],'/',de,'/signif_genes_FDR_0.01.csv')) %>% as.data.frame()
    
    de_cis <- de_df %>%
      dplyr::filter(grepl('In_cis',Gene_Type))
    same_direction <- c('In_cis_Increase_UpRegulated','In_cis_Decrease_DownRegulated')
    opposite_direction <- c('In_cis_Increase_DownRegulated','In_cis_Decrease_UpRegulated')
    # de_df <- de_df %>%
    #   dplyr::filter(Gene_Type %in% same_direction) %>%
    #   # inner_join(meta_genes, by='ensembl_gene_id') %>%
    #   # dplyr::filter(chr %in% chrs) %>%
    #   dplyr::summarise(pct_genes=round(n()*100/dim(de_cis)[1],1)) %>%
    #   # dplyr::select(ref_gene_symbol, logFC, Gene_Type) %>%
    #   dplyr::mutate(desc=de, datatag=pair_df[de,'datatag']) #logFC=round(logFC, 1),
    
    de_df <- de_df %>%
      dplyr::filter(Gene_Type %in% opposite_direction) %>%
      # inner_join(meta_genes, by='ensembl_gene_id') %>%
      # dplyr::filter(chr %in% chrs) %>%
      dplyr::summarise(pct_genes=round(n()*100/dim(de_cis)[1],1)) %>%
      # dplyr::select(ref_gene_symbol, logFC, Gene_Type) %>%
      dplyr::mutate(desc=de, datatag=pair_df[de,'datatag']) #logFC=round(logFC, 1),
    
    de_genes <- dplyr::bind_rows(de_genes, de_df)
  }
  dim(de_genes)
  get_stat(de_genes$pct_genes, 'in cis same direction, positive linear influence genes', 'In all comparisons')
  get_stat(de_genes$pct_genes, 'in cis opposite direction, negative linear influence genes', 'In all comparisons')
  de_genes_UnRx <- de_genes %>%
    dplyr::inner_join(pair_UnRx, by='desc')
  de_genes_Rx_UnRx <- de_genes %>%
    dplyr::inner_join(pair_Rx_UnRx, by='desc')
  get_stat(de_genes_UnRx$pct_genes, 'in cis same direction, positive linear influence genes', 'In untreated comparisons')
  get_stat(de_genes_Rx_UnRx$pct_genes, 'in cis same direction, positive linear influence genes', 'In treated vs. untreated comparisons')
  
  get_stat(de_genes_UnRx$pct_genes, 'in cis opposite direction, negative linear influence genes', 'In untreated comparisons')
  get_stat(de_genes_Rx_UnRx$pct_genes, 'in cis opposite direction, negative linear influence genes', 'In treated vs. untreated comparisons')
  
  
  for(dtag in unique(de_genes1$datatag)){
    print("------------------------")
    print(dtag)
    tmp <- de_genes1 %>% 
      dplyr::filter(datatag==dtag) 
    if(dtag=='SA535'){
      tmp <- tmp %>% 
        dplyr::filter(desc!="Res(X9:A) vs. Sen(X9:G)") 
    }
    signf_genes <- tmp %>% 
      dplyr::group_by(ref_gene_symbol)%>% 
      dplyr::summarise(nb_occurence=n())%>%
      dplyr::filter(nb_occurence>1)%>%
      dplyr::pull(ref_gene_symbol)  
    
    obs_gene_type <- 'In_cis_Increase_UpRegulated'
    print(obs_gene_type)
    iu_genes <- tmp %>% 
      dplyr::mutate(desc=paste0(desc, ': logFC')) %>% 
      dplyr::filter(ref_gene_symbol %in% signf_genes & Gene_Type==obs_gene_type) %>% 
      dplyr::select(-datatag, -Gene_Type)%>% 
      tidyr::pivot_wider(names_from='desc', values_from='logFC', values_fill = NA)%>%
      tibble::column_to_rownames('ref_gene_symbol')
    iu_genes <- iu_genes[order(rowSums(abs(iu_genes)), decreasing=T),]
    iu_genes$gene_symbol <- rownames(iu_genes)
    data.table::fwrite(iu_genes, paste0(csv_dir, dtag,'_',obs_gene_type,'.csv'))
    print(dim(iu_genes))  
    
    obs_gene_type <- 'In_cis_Decrease_DownRegulated'
    print(obs_gene_type)
    dd_genes <- tmp %>% 
      dplyr::mutate(desc=paste0(desc, ': logFC')) %>% 
      dplyr::filter(ref_gene_symbol %in% signf_genes & Gene_Type==obs_gene_type) %>% 
      dplyr::select(-datatag, -Gene_Type)%>% 
      tidyr::pivot_wider(names_from='desc', values_from='logFC', values_fill = NA) %>%
      tibble::column_to_rownames('ref_gene_symbol')
  
    dd_genes <- dd_genes[order(rowSums(abs(dd_genes)), decreasing=T),]
    dd_genes$gene_symbol <- rownames(dd_genes)
    data.table::fwrite(dd_genes, paste0(csv_dir, dtag,'_',obs_gene_type,'.csv'))
    print(dim(dd_genes))
    print("------------------------")
  }   
  
  
  obs_gene_types <- c('In_cis_Increase_UpRegulated', 'In_cis_Decrease_DownRegulated')
  dtag <- 'SA609'
  sa609_iu <- data.table::fread(paste0(csv_dir, dtag,'_',obs_gene_types[1],'.csv')) %>% as.data.frame()
  sa609_iu$Gene_Type_Pt3 <- obs_gene_types[1]
  genes_SA609 <- sa609_iu$gene_symbol
  sa609_dd <- data.table::fread(paste0(csv_dir, dtag,'_',obs_gene_types[2],'.csv')) %>% as.data.frame()
  sa609_dd$Gene_Type_Pt3 <- obs_gene_types[2]
  genes_SA609 <- unique(c(genes_SA609, sa609_dd$gene_symbol))
  length(genes_SA609)
  dtag <- 'SA535'
  sa535_iu <- data.table::fread(paste0(csv_dir, dtag,'_',obs_gene_types[1],'.csv')) %>% as.data.frame()
  sa535_iu$Gene_Type_Pt4 <- obs_gene_types[1]
  genes_SA535 <- sa535_iu$gene_symbol
  sa535_dd <- data.table::fread(paste0(csv_dir, dtag,'_',obs_gene_types[2],'.csv')) %>% as.data.frame()
  sa535_dd$Gene_Type_Pt4 <- obs_gene_types[2]
  genes_SA535 <- unique(c(genes_SA535, sa535_dd$gene_symbol))
  length(genes_SA535)
  common_genes <- intersect(genes_SA609, genes_SA535)    
  length(common_genes)
  sa609_iu <- sa609_iu %>%
    dplyr::filter(gene_symbol %in% common_genes)
  sa609_dd <- sa609_dd %>%
    dplyr::filter(gene_symbol %in% common_genes)
  sa535_iu <- sa535_iu %>%
    dplyr::filter(gene_symbol %in% common_genes)
  sa535_dd <- sa535_dd %>%
    dplyr::filter(gene_symbol %in% common_genes)
  
  common_genes_df1 <- sa609_iu %>% inner_join(sa535_dd, by='gene_symbol')
  common_genes_df2 <- sa609_dd %>% inner_join(sa535_dd, by='gene_symbol')
  View(common_genes_df1)
  View(common_genes_df2)
  data.table::fwrite(common_genes_df1,paste0(csv_dir, 'common_genes_Pt3_Pt4_part1.csv'))
  data.table::fwrite(common_genes_df2,paste0(csv_dir, 'common_genes_Pt3_Pt4_part2.csv'))
}

get_stat <- function(vals, tag, datatag){
  print(paste0(datatag, ': ',tag,':  mean:',round(mean(vals),1),' sd:',round(sd(vals),1),' IQR:',round(IQR(vals),1)))
}
get_chrs_distribution <- function(){
  # library(dplyr)
  # library(ggplot2)
  # Reference genes, chrs
  chrs <- c(as.character(1:22), "X")
  meta_genes <- annotables::grch38 %>%
    dplyr::select(ensembl_gene_id=ensgene, ref_gene_symbol=symbol, chr) 
  meta_genes <- meta_genes[!duplicated(meta_genes$ensembl_gene_id), ]
  # dim(meta_genes)
  results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
  
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
  pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
  pair_df <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  rownames(pair_df) <- pair_df$desc
  
  save_dir <- input_dir
 
  pair_df <- pair_df %>% 
    dplyr::filter(comp_type=='treated_vs_untreated')%>% 
    dplyr::select(-result_fn)
  dim(pair_df)
  
  results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
  stat <- tibble::tibble()
  
  for(de in pair_df$desc){
    de_df <- data.table::fread(paste0(results_dir,pair_df[de,'datatag'],'/',de,'/signif_genes_FDR_0.01.csv')) %>% as.data.frame()
    
    de_df <- de_df %>%
      dplyr::filter(Gene_Type %in% c('In_cis_Increase_UpRegulated','In_cis_Decrease_DownRegulated')) %>%
      inner_join(meta_genes, by='ensembl_gene_id') %>%
      dplyr::filter(chr %in% chrs) %>%
      dplyr::group_by(chr, Gene_Type) %>%
      dplyr::summarise(nb_genes=n())%>%
      dplyr::ungroup()
    
    de_df$Gene_Type <- ifelse(de_df$Gene_Type=='In_cis_Increase_UpRegulated','IncisIU','IncisDD')
    de_df$desc <- de
    de_df$datatag <- pair_df[de,'datatag']
    stat <- dplyr::bind_rows(stat, de_df)
    # dim(de_df)
    # View(de_df)
    # de_df <- de_df %>% 
    #   tidyr::pivot_wider(names_from = 'Gene_Type', values_from = 'nb_genes', values_fill=0)
  }
  dim(stat)
  plt_desc_df <- data.table::fread(paste0(save_dir,'plt_desc_v2.csv')) %>% as.data.frame()
  dim(plt_desc_df)
  stat$plt_desc[1]
  stat$Gene_Type
  csv_dir <- paste0(results_dir, 'chr_distribution/')
  dir.create(csv_dir)
  # dtag <- 'SA535'
  stat <- stat %>%
    dplyr::inner_join(plt_desc_df, by='desc') 
  data.table::fwrite(stat, paste0(csv_dir, 'stat_chrs_distribution.csv'))
  stat <- data.table::fread(paste0(csv_dir, 'stat_chrs_distribution.csv'))

  for(dtag in unique(stat$datatag)){
    stat_pt <- stat %>%
      dplyr::filter(datatag==dtag) #%>%
    # dplyr::select(-Gene_Type)
    dim(stat_pt)
    if(dtag=='SA535'){
      stat_pt <- stat_pt %>% 
      dplyr::filter(plt_desc!="Res(X9:A) vs. Sen(X9:G)") 
    }
    stat_pt12 <- stat_pt %>%
      dplyr::group_by(plt_desc, chr) %>%
      dplyr::summarise(total_genes=sum(nb_genes)) %>%
      dplyr::mutate(plt_desc=paste0(plt_desc,': total genes')) %>%
      tidyr::pivot_wider(names_from = 'plt_desc', values_from = 'total_genes', values_fill=0) %>%
      tibble::column_to_rownames('chr')
    
    stat_pt12$avg_nb_genes <- round(rowSums(stat_pt12)/ncol(stat_pt12),1)
    stat_pt12$chr <- rownames(stat_pt12)
    
    stat_pt11 <- stat_pt %>%
      dplyr::mutate(plt_desc=paste0(plt_desc,': ',Gene_Type)) %>%
      dplyr::select(-desc, -datatag, -Gene_Type) %>%
      tidyr::pivot_wider(names_from = 'plt_desc', values_from = 'nb_genes', values_fill=0)%>%
      inner_join(stat_pt12, by='chr')
    
    stat_pt11$chr <- factor(stat_pt11$chr, levels = chrs)
    stat_pt11 <- stat_pt11[gtools::mixedorder(stat_pt11$chr),]
    data.table::fwrite(stat_pt11, paste0(csv_dir, dtag,'_chrs_distribution.csv'))
    
  }
  
  for(dtag in unique(stat$datatag)){
    print("------------------------")
    print(dtag)
    pt <- data.table::fread(paste0(csv_dir, dtag,'_chrs_distribution.csv')) %>% as.data.frame()
    pt <- pt[order(pt$avg_nb_genes, decreasing=T),]
    print('Top 5 chrs with highest number of positive influence genes (IU, DD): ')
    print(pt[1:5,'chr'])
    
    comps <- colnames(pt)
    compsIU <- grep('IncisIU',comps, value=T)
    pt$avg_IU <- 0
    for(c in compsIU){
      pt$avg_IU <- pt$avg_IU + pt[,c]
    }
    pt$avg_IU <- round(pt$avg_IU / length(compsIU), 1)
    pt <- pt[order(pt$avg_IU, decreasing=T),]
    print('Top 5 chrs with highest number of IU genes: ')
    print(pt[1:5,'chr'])
    
    compsDD <- grep('IncisDD',comps, value=T)
    pt$avg_DD <- 0
    for(c in compsDD){
      pt$avg_DD <- pt$avg_DD + pt[,c]
    }
    pt$avg_DD <- round(pt$avg_DD / length(compsDD), 1)
    pt <- pt[order(pt$avg_DD, decreasing=T),]
    print('Top 5 chrs with highest number of DD genes: ')
    print(pt[1:5,'chr'])
    print("------------------------")
  }  
  
 
  
  
}

results_dir
pw <- data.table::fread(paste0(results_dir,'pathways_IU_DD/pathways_summary.csv'))
View(pw)
pwDD <- pw %>%
  filter(grepl('DD',datatag)) #%>%
  # select(pathway, datatag)

pwIU <- pw %>%
  filter(grepl('IU',datatag))
dim(pwDD)
View(pwDD)
sum(pwIU$padj<0.05)
dim(pwIU)
View(pwIU)

# TO DO: a table of chrs, nb of cis genes. 

