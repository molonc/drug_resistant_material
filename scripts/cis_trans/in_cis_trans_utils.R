suppressPackageStartupMessages({
  library(gtools)
  library(parallel)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggrepel)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
})

# library(extrafont)
# font_import(prompt=F) # import all your fonts
# fonts()
# library(extrafont)
# font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
# fonts()

# stat_df: ks statistical tests, and metasample info
# de_df: df of logFC values
get_boxplot_ksStat <- function(stat_df, de_df, save_dir, datatag='allSeries', plottype='boxplot'){
  # library(dplyr)
  # library(ggplot2)
  yl <- c(-3.5, 3.5)
  de_df <- de_df %>%
    dplyr::mutate(logFC=case_when(
      logFC < yl[1] ~ yl[1],
      logFC > yl[2] ~ yl[2],
      TRUE ~ logFC
    ))
  de_df <- de_df %>%
    dplyr::mutate(classified_gene=case_when(
      classified_gene =='in_cis' ~ "In cis",
      classified_gene =='in_trans' ~ "In trans",
      TRUE ~ classified_gene
    ))
  GT <- c("chocolate", "blue2")
  names(GT) <- c("In cis", "In trans")
  stat_df$labels_detailed <- gsub(' vs. ','\nvs. ',stat_df$labels_detailed)
  meta_info <- stat_df %>%
    select(file_header,labels_detailed, ks_stat, series, order)
  de_df <- de_df %>%
    left_join(meta_info, by='file_header') %>%
    arrange(order)
  
  # de_df$labels_detailed <- gsub(' vs. ','\nvs. ',de_df$labels_detailed)
  de_df$labels_detailed <- factor(de_df$labels_detailed, levels=unique(de_df$labels_detailed))
  
  p <- ggplot(de_df, aes(x=labels_detailed, y=logFC, fill=classified_gene)) 
  if(plottype=='boxplot'){
    p <- p + geom_boxplot(position=position_dodge(1), outlier.size = 0.02) 
  }else{
    p <- p + geom_violin(draw_quantiles = c(0.5))
  }
  p <- p + scale_fill_manual(values = GT, name='Gene type') +
    facet_grid(. ~ series, scales="free", space='free') +
    theme(panel.spacing = unit(1, "lines"),
          axis.text.x = element_text(color="black", size=9, hjust = 0.5, angle = 90)) + 
    ylim(yl[1],yl[2]+0.3) + 
    labs(x=NULL,y='log2 FC',title='Degree of gene type expression')
  # p
  ## converting ks signf to star symbol
  signif_level <- stats::symnum(stat_df$ks_stat, corr = FALSE, na = FALSE, 
                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                symbols = c("***", "**", "*", ".", " "))
  stat_df$signif_level <- as.character(signif_level)
  stat_df$classified_gene <- 'in_cis'
  
  p1 <- p + geom_text(data=stat_df, aes(x = labels_detailed, y = yl[2]+0.2, label=signif_level,
  ), hjust=0.5, size=4) #, hjust=hjust, vjust=vjust
  # p1
  
  my_font <- "Helvetica"
  thesis_theme <- ggplot2::theme(
    text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
    axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),  
    axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(color="black",size=6, hjust = 0.5, family=my_font, angle = 90),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(color="black",size=7, hjust = 0.5, family=my_font),
    plot.title = element_text(color="black",size=10, face="bold", hjust=0, family=my_font),
    legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font), 
    legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
    strip.text.x = element_text(color="black",size=9, family=my_font),
    strip.text.y = element_text(color="black",size=9, family=my_font),
    legend.spacing.x = unit(0.1, 'mm'),
    legend.spacing.y = unit(0.1, 'mm'),
    legend.key.height=unit(1,"line"),
    legend.position = 'bottom',
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  p1 <- p1 + thesis_theme
  # png(paste0(save_dir,datatag,'_ks_test.png'), height = 2*500, width=2*1000,res = 2*72)
  # print(p1)
  # dev.off()
  
  ggsave(paste0(save_dir,datatag,'_ks_test.png'),
         plot = p1,
         height = 5,
         width = 8,
         # useDingbats=F,
         dpi=250)
  
  return(p1)
}


get_pathway_gprofiler <- function(pair_groups, base_dir, datatag){
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
  pair_groups <- data.table::fread(paste0(input_dir,'comparisons_drug_res.csv')) %>% as.data.frame()
  View(pair_groups)
  pair_groups <- pair_groups %>%
    dplyr::filter(comp_type!='untreated' & datatag!='SA1035')
  pathway_stat <- tibble::tibble()
  rownames(pair_groups) <- pair_groups$result_fn
  
  for(fn in pair_groups$result_fn){
    # fn <- gsub('_logfc_results.csv','',fn)
    de_genes <- data.table::fread(paste0(input_dir, fn)) %>% as.data.frame()
    de_genes <- data.table::fread(paste0(base_dir, fn,'/','signif_genes_FDR0.01.csv')) %>% as.data.frame()
    print(dim(de_genes))
    de_genes <- de_genes %>%
      dplyr::filter(abs(logFC)>0.25 & FDR<0.01 & PValue<0.05)
    print(dim(de_genes))
    stat <- get_gprofiler_pathways_obsgenes(de_genes$en_gene_id, save_dir, datatag)
    if(!is.null(gsea_out)){
      signf_ls <- dplyr::bind_rows(pathway_stat, stat)  
    }
  }
  
  
  pathway_stat <- as.data.frame(pathway_stat)
  
  
  pathway_stat <- pathway_stat %>%
    select(pathway, desc, padj, datatag, everything())
  
  data.table::fwrite(pathway_stat, paste0(input_dir,datatag,'_reference_sets.csv'))
  datatag <- ''
  data.table::fwrite(pathway_stat, paste0(input_dir,datatag,'_pathways.csv'))
  return(pathway_stat)
}

get_reference_set_stat_v2 <- function(pair_groups, base_dir, datatag){
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
  pair_groups <- data.table::fread(paste0(input_dir,'comparisons_drug_res.csv')) %>% as.data.frame()
  View(pair_groups)
  pair_groups <- pair_groups %>%
    dplyr::filter(comp_type!='untreated' & datatag!='SA1035')
  signf_ls <- tibble::tibble()
  fgs <- c("COSMIC","CISPLATIN_RESISTANCE","CORE_FITNESS","BS_ESSENTIAL_CANCER") #,"METASTASIS"
  reference_genes_set <- 'COSMIC'
  rownames(pair_groups) <- pair_groups$result_fn
  for(reference_genes_set in fgs){
    for(fn in pair_groups$result_fn){
      # fn <- gsub('_logfc_results.csv','',fn)
      de_genes <- data.table::fread(paste0(input_dir, fn)) %>% as.data.frame()
      de_genes <- de_genes %>%
        dplyr::filter(abs(logFC)>0.25 & FDR<0.01 & PValue<0.05)
      print(dim(de_genes))
      gsea_out <- get_custom_pathway_results(de_genes,                      # named vector of statistical significance 
                                             fn, base_name=pair_groups[fn,'datatag'],   
                                             pathway_name='custom_pathways', 
                                             groups_use=NULL,    # vector of 2 elements, 2 group name used for DE analysis
                                             save_dir=paste0(input_dir),
                                             reference_genes_set)
      if(!is.null(gsea_out)){
        signf_ls <- dplyr::bind_rows(signf_ls, gsea_out)  
      }
    }
  }
  
  signf_ls <- as.data.frame(signf_ls)
  
  
  signf_ls <- signf_ls %>%
    select(pathway, desc, padj, datatag, everything())
  
  data.table::fwrite(signf_ls, paste0(input_dir,datatag,'_reference_sets.csv'))
  datatag <- ''
  data.table::fwrite(signf_ls, paste0(input_dir,datatag,'_cisplatin_resistance_609_535.csv'))
  return(signf_ls)
}

get_reference_set_stat <- function(pair_groups, base_dir, datatag){
  signf_ls <- tibble::tibble()
  fgs <- c("COSMIC","CISPLATIN_RESISTANCE","CORE_FITNESS","BS_ESSENTIAL_CANCER") #,"METASTASIS"
  for(reference_genes_set in fgs){
    for(fn in pair_groups$result_fn){
      # fn <- gsub('_logfc_results.csv','',fn)
      de_genes <- data.table::fread(paste0(base_dir, fn,'/','signif_genes_FDR0.01.csv')) %>% as.data.frame()
      print(dim(de_genes))
      gsea_out <- get_custom_pathway_results(de_genes,                      # named vector of statistical significance 
                                             fn, base_name=datatag,   
                                             pathway_name='custom_pathways', 
                                             groups_use=NULL,    # vector of 2 elements, 2 group name used for DE analysis
                                             save_dir=paste0(base_dir, fn,'/'),
                                             reference_genes_set)
      if(!is.null(gsea_out)){
        signf_ls <- dplyr::bind_rows(signf_ls, gsea_out)  
      }
    }
  }
  
  signf_ls <- as.data.frame(signf_ls)
  signf_ls <- signf_ls %>%
    select(pathway, padj, desc, datatag, everything())
  data.table::fwrite(signf_ls, paste0(base_dir,datatag,'_reference_sets.csv'))
  return(signf_ls)
}
get_gene_type_edgeR <- function(pair_groups_fn, cnv_mat, 
                                datatag, subtag=NULL, input_dir, 
                                cancer_ref_genes_fn=NULL, outlier_FC_thrs=4.25,
                                FDR_cutoff=0.01, minLogFC=0.25, pValueThrs=0.05){  #FDR_cutoff=0.01, FDR_cutoff=0.05
  
  
  if(is.null(cancer_ref_genes_fn)){
    ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
    cancer_ref_genes_fn <- paste0(ref_dif,'Behan_CFgenes.csv')
  }
  cancer_ref_genes_df <- read.csv(cancer_ref_genes_fn, stringsAsFactors=F, check.names = F)
  print(dim(cancer_ref_genes_df))
  colnames(cancer_ref_genes_df)[which(colnames(cancer_ref_genes_df) == "ADAM PanCancer Core-Fitness genes")] <- "PanCancer_Fitness_genes"
  
  pair_groups <- read.csv(pair_groups_fn, header=T, check.names=F, stringsAsFactors=F)
  pair_groups <- pair_groups[pair_groups$datatag==datatag,]
  # View(pair_groups1)
  # if(!is.null(subtag)){
  #   pair_groups <- pair_groups[pair_groups$title==subtag,]
  #   # if(subtag=='SA535_cisplatin'){
  #   #   pair_groups$desc[6] <- paste0(pair_groups$desc[6],'_2')
  #   # }
  # }
  # print(dim(pair_groups))
  # print(pair_groups)
  
  
  rownames(pair_groups) <- pair_groups$desc
  frac_dlp <- list()
  genes_summary_stat <- list()
  genes_summary <- list()
  c <- 0
  for(de in pair_groups$desc){
    output_dir <- paste0(input_dir,de,"/")
    if (!file.exists(output_dir)){
      dir.create(output_dir, recursive=T)
    }
    # signif_genes <- read.csv(paste0(output_dir,'edgeR_significant_genes.csv'), check.names=F, stringsAsFactors=F)
    signif_genes <- read.csv(paste0(input_dir,pair_groups[de,'result_fn']), check.names=F, stringsAsFactors=F)
    
    print(dim(signif_genes))
    signif_genes <- signif_genes %>%
                    dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
    
    print(dim(signif_genes))
    # colnames(signif_genes)[which(names(signif_genes)=='gene_id')] <- 'ensembl_gene_id'
    # print(colnames(signif_genes))
    # signif_genes <- signif_genes[abs(signif_genes$logFC)>0.25,]
    signif_genes <- signif_genes %>% left_join(cancer_ref_genes_df, by = c("ensembl_gene_id"="ensemble_id"))
    s1 <- as.character(pair_groups[de,'clone1'])
    s2 <- as.character(pair_groups[de,'clone2'])
    if(sum(colnames(cnv_mat) %in% c(s1,s2))!=2){
      print(de)
      stop("DEBUG: do not exist 2 clones in cnv median CN profiles")
    }
    cnv_mat_tmp <- cnv_mat[,c(s1, s2)]
    dim(cnv_mat_tmp)
    var_genes <- apply(cnv_mat_tmp, 1, var)
    cnv_mat_tmp <- cnv_mat_tmp[var_genes > 0,]
    sum(is.na(cnv_mat_tmp$R))
    sum(is.na(cnv_mat_tmp$H))
    
    for(g in rownames(cnv_mat_tmp)){
      if(cnv_mat_tmp[g,s1]-cnv_mat_tmp[g,s2]>0){
        cnv_mat_tmp[g,'classified_gene_dlp'] <- 'Increase'
      }else if(cnv_mat_tmp[g,s1]-cnv_mat_tmp[g,s2]<0){
        cnv_mat_tmp[g,'classified_gene_dlp'] <- 'Decrease'
      }else{
        cnv_mat_tmp[g,'classified_gene_dlp'] <- 'No_variance'
      }
    }  
    print(dim(cnv_mat_tmp))
    print(summary(as.factor(cnv_mat_tmp$classified_gene_dlp)))
    cnv_mat_tmp <- cnv_mat_tmp[cnv_mat_tmp$classified_gene_dlp != 'No_variance',]
    print(dim(cnv_mat_tmp))
    print(paste0('Nb of variance dlp genes: ',nrow(cnv_mat_tmp)))  #675 genes with variance
    cnv_mat_tmp$ensembl_gene_id <- rownames(cnv_mat_tmp)
    # write.csv(cnv_mat_tmp, paste0(output_dir, "cnv_mat_total.csv"),quote=F, row.names=F)
    signif_genes <- signif_genes %>% left_join(cnv_mat_tmp, by = c("ensembl_gene_id"))
    in_cis <- intersect(signif_genes$ensembl_gene_id, cnv_mat_tmp$ensembl_gene_id)
    # print(paste0('Nb of variance in-cis genes: ',length(in_cis)))
    rownames(signif_genes) <- signif_genes$ensembl_gene_id
    signif_genes$classified_gene <- 'in_trans'
    signif_genes[in_cis,'classified_gene'] <- 'in_cis'
    # print(summary(signif_genes$logFC))
    
    signif_genes_outliers <- signif_genes[abs(signif_genes$logFC)>outlier_FC_thrs,]
    if(nrow(signif_genes_outliers)>0){
      # write.csv(signif_genes_outliers, paste0(output_dir,'signif_genes_outliers.csv'), quote = F, row.names = F)
    }
    # signif_genes <- signif_genes[abs(signif_genes$log2FoldChange) <= outlier_FC_thrs,]
    
    
    signif_genes$classified_gene_rna <- ifelse(signif_genes$logFC>0,'UpRegulated',
                                               ifelse(signif_genes$logFC<0,'DownRegulated','Invalid'))
    signif_genes$Gene_Type <- ifelse(!is.na(signif_genes$classified_gene_dlp),
                                     paste0('In_cis_',signif_genes$classified_gene_dlp, '_',signif_genes$classified_gene_rna),
                                     paste0('In_trans_',signif_genes$classified_gene_rna))
    
    signif_genes$Gene_Type <- as.factor(signif_genes$Gene_Type)
    
    signif_genes$is_fitness_gene <- as.factor(!is.na(signif_genes$PanCancer_Fitness_genes))
    
    # plot_DE(de, signif_genes, pair_groups, output_dir,
    #         minLogFC=0.25, pValueThrs=0.05, nbtopup=30, nbtopdown=30)
    
    # observed_signif_genes <- signif_genes[abs(signif_genes$logFC)>0.25 
    #                                       & signif_genes$FDR<0.01
    #                                       & signif_genes$PValue<0.05,]
    observed_signif_genes <- signif_genes
    print(paste0('Remove nb genes with logFC < 0.25: ',sum(signif_genes$logFC<=0.25)))
    
    print(paste0('Comparison: ',de,"  Number of fitness genes: "))
    print(summary(as.factor(observed_signif_genes$is_fitness_gene)))
    
    # Summary incis genes
    nb_incr_incis <- sum(observed_signif_genes$classified_gene_dlp=='Increase', na.rm = T)
    nb_decr_incis <- sum(observed_signif_genes$classified_gene_dlp=='Decrease', na.rm = T)
    frac_dlp[[as.character(de)]] <- c(pair_groups[de,'result_fn'],
                                      paste0(datatag,'-',s1,'_vs_',s2),
                                      nrow(cnv_mat_tmp),length(in_cis),
                                      nb_incr_incis,
                                      nb_decr_incis,
                                      round(length(in_cis)/nrow(cnv_mat_tmp) * 100,2))
   
    s <- as.data.frame(summary(as.factor(observed_signif_genes$Gene_Type)))
    # summary(-log10(signif_genes$pvalue))
    # summary(signif_genes$padj)
    colnames(s)[1] <- 'nb_genes'
    # View(s)
    s$pct_genes <- round(s$nb_genes/colSums(s), 3) * 100
    s$datatag <- pair_groups[de,'datatag']
    s$desc <- as.character(de)
    s$result_fn <- pair_groups[de,'result_fn']
    s$de_analysis <- paste0(datatag,'-',s1,'_vs_',s2)
    # colnames(s)[which(colnames(s) == "pct_genes")] 
    # colnames(s)[which(colnames(s) == "nb_genes")] <- paste0(datatag,'-',s1,'_vs_',s2,'_nbgenes')
    c <- c + 1
    s$gene_type <- rownames(s)
    genes_summary_stat[[c]] <- s
    
    write.csv(s, paste0(output_dir,de,'.csv'), quote = F, row.names = F)
    write.csv(signif_genes, paste0(output_dir, 'signif_genes.csv'), quote = F, row.names = F)
    genes_summary[[paste0(de,'_signif_genes')]] <- signif_genes
    genes_summary[[paste0(de,'_stat')]] <- s
    
  }
  frac_dlp_df <- as.data.frame(do.call(rbind, frac_dlp))
  colnames(frac_dlp_df) <- c('result_fn','de_analysis','total_var_dlp',
                             'incis','incis_incr','incis_decr','pct_incis')# 
  # frac_dlp_df$desc <- rownames(frac_dlp_df)
  # frac_dlp_df$pct_incis_incr <- round(as.numeric(frac_dlp_df$incis_incr)/as.numeric(frac_dlp_df$total_var_dlp) * 100, 2)
  # frac_dlp_df$pct_incis_decr <- round(as.numeric(frac_dlp_df$incis_decr)/as.numeric(frac_dlp_df$total_var_dlp) * 100, 2)
  # frac_dlp_df$pct_others <- 100 - frac_dlp_df$pct_incis_incr - frac_dlp_df$pct_incis_decr
  # View(frac_dlp_df)
  # frac_dlp_df$pct_incis_inc <- round(as.numeric(frac_dlp_df$incis_incr) / as.numeric(frac_dlp_df$total_var_dlp) * 100,2)
  # frac_dlp_df$pct_incis_dec <- round(as.numeric(frac_dlp_df$incis_decr) / as.numeric(frac_dlp_df$total_var_dlp) * 100,2)
  # frac_dlp_df$unaffected <- 100 - frac_dlp_df$pct_incis_inc - frac_dlp_df$pct_incis_dec
  # write.csv(frac_dlp_df, paste0(input_dir,datatag,'_',subtag,'_fraction_dlp.csv'), quote = F, row.names = F)
  
  stat <- do.call(rbind, genes_summary_stat)
  # View(head(stat))
  # stat$gene_type <- rownames(stat)
  # write.csv(stat, paste0(input_dir,datatag,'_',subtag,'_summary.csv'), quote = F, row.names = F)
  
  # saveRDS(genes_summary,paste0(input_dir, datatag,'_',subtag,'_summary_statistic.rds'))
}

get_gene_type_stat <- function(signif_genes_df, cnv_mat, datatag, output_dir, outlier_FC_thrs=3)
{
  signif_genes_df <- signif_genes_df %>%
    dplyr::filter(abs(logFC)>0.25 & PValue <0.05)
  print(dim(signif_genes_df))
  # dir.create(output_dir, showWarnings = F)
  output_dir <- paste0(output_dir,datatag,'_result/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  reference_genes_df <- read.csv(paste0(ref_dif,'Symbol_ensembl.csv'), 
                                 stringsAsFactors=F, check.names = F)
  
  cancer_ref_genes_df <- read.csv(paste0(ref_dif,'Behan_CFgenes.csv'), 
                                 stringsAsFactors=F, check.names = F)
  dim(cancer_ref_genes_df)
  colnames(cancer_ref_genes_df)[which(colnames(cancer_ref_genes_df) == "ADAM PanCancer Core-Fitness genes")] <- "PanCancer_Fitness_genes"
  
  
  signif_genes_df <- signif_genes_df %>% left_join(reference_genes_df, by = c("gene_id"="Ensembl"))
  signif_genes_df <- signif_genes_df %>% left_join(cancer_ref_genes_df, by = c("gene_id"="ensemble_id"))
  signif_genes_df$de_analysis <- paste0(signif_genes_df$sample1,'_vs_',signif_genes_df$sample2)
  
  de_pairs <- unique(signif_genes_df$de_analysis)
  genes_summary_stat <- list()
  genes_summary <- list()
  c <- 0
  for(de in de_pairs){
    signif_genes <- signif_genes_df[signif_genes_df$de_analysis==de,]
    print(dim(signif_genes))
    s1 <- as.character(signif_genes$sample1[1])
    s2 <-  as.character(signif_genes$sample2[1])
    print(paste0('clone ',s1, ' versus clone ',s2, ' nb significant genes: ',nrow(signif_genes)))
    print(sum(colnames(cnv_mat) %in% c(s1, s2))==2)
    cnv_mat_tmp <- cnv_mat[,c(s1, s2)]
    print(dim(cnv_mat_tmp))
    cnv_mat_tmp <- cnv_mat_tmp[(cnv_mat_tmp[,s1]-cnv_mat_tmp[,s2])!=0,]
    for(g in rownames(cnv_mat_tmp)){
      if(cnv_mat_tmp[g,s1]-cnv_mat_tmp[g,s2]>0){
        cnv_mat_tmp[g,'classified_gene_dlp'] <- 'Increase'
      }else if(cnv_mat_tmp[g,s1]-cnv_mat_tmp[g,s2]<0){
        cnv_mat_tmp[g,'classified_gene_dlp'] <- 'Decrease'
      }else{
        cnv_mat_tmp[g,'classified_gene_dlp'] <- 'No_variance'
      }
    }  
    print(summary(as.factor(cnv_mat_tmp$classified_gene_dlp)))
    cnv_mat_tmp <- cnv_mat_tmp[cnv_mat_tmp$classified_gene_dlp != 'No_variance',]
    print(dim(cnv_mat_tmp))
    print(paste0('Nb of variance dlp genes: ',nrow(cnv_mat_tmp)))
    cnv_mat_tmp$gene_id <- rownames(cnv_mat_tmp)
    # signif_genes <- signif_genes_df %>% as.data.frame()
    signif_genes <- signif_genes %>% left_join(cnv_mat_tmp, by = c("gene_id"))
    # View(head(signif_genes))
    # dim(signif_genes)
    in_cis <- intersect(signif_genes$gene_id, cnv_mat_tmp$gene_id)
    
    print(paste0('Nb of variance in-cis genes: ',length(in_cis)))
    
    rownames(signif_genes) <- signif_genes$gene_id
    signif_genes$classified_gene <- 'in_trans'
    signif_genes[in_cis,'classified_gene'] <- 'in_cis'
    print(summary(signif_genes$logFC))
    signif_genes_outliers <- signif_genes[abs(signif_genes$logFC)>outlier_FC_thrs,]
    write.csv(signif_genes_outliers, paste0(output_dir,datatag,'_',de,'_signif_genes_outliers.csv'), quote = F, row.names = F)
    # signif_genes <- signif_genes[abs(signif_genes$logFC) <= outlier_FC_thrs,]
    
    
    signif_genes$classified_gene_rna <- ifelse(signif_genes$logFC>0,'UpRegulated',
                                               ifelse(signif_genes$logFC<0,'DownRegulated','Invalid'))
    signif_genes$Gene_Type <- ifelse(!is.na(signif_genes$classified_gene_dlp),
                                     paste0('In_cis_',signif_genes$classified_gene_dlp, '_',signif_genes$classified_gene_rna),
                                     paste0('In_trans_',signif_genes$classified_gene_rna))
    
    signif_genes$Gene_Type <- as.factor(signif_genes$Gene_Type)
    p <- ggplot(signif_genes) +
      geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour=Gene_Type),size = 2) +  #aes(colour = factor(is_fitness_gene)), 
      scale_color_discrete(labels = paste0(levels(signif_genes$Gene_Type),': ',table(signif_genes$Gene_Type)))
    
    p <- p + labs(title=paste0(datatag,' ',de))
    p <- p + theme(plot.title = element_text(color="black", size=12, hjust = 0.5),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   axis.text.y = element_text(color="black", size=10),
                   axis.text.x = element_text(color="black", size=10),
                   axis.title = element_text(color="black", size=10),
                   legend.text = element_text(size=7))
    # p <- ggplot(signif_genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
    #   geom_point(aes(colour = Gene_Type), size = 2)
    # 
    png(paste0(output_dir,datatag,'_',de,'_signif_genes.png'), height = 2*400, width=2*500,res = 2*72)
    print(p)
    dev.off()
    
    signif_genes$is_fitness_gene <- as.factor(!is.na(signif_genes$PanCancer_Fitness_genes))
    print(paste0('Comparison: ',de,"  Number of fitness genes: "))
    print(summary(as.factor(signif_genes$is_fitness_gene)))
    
    p_ref <- ggplot(signif_genes) +
      geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour=is_fitness_gene),size = 2) +  #aes(colour = factor(is_fitness_gene)), 
      scale_color_discrete(labels = paste0(levels(signif_genes$is_fitness_gene),': ', table(signif_genes$is_fitness_gene)))
    p_ref <- p_ref + labs(title=paste0(datatag,' ',de))
    p_ref <- p_ref + theme(plot.title = element_text(color="black", size=12, hjust = 0.5),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                           axis.text.y = element_text(color="black", size=10),
                           axis.text.x = element_text(color="black", size=10),
                           axis.title = element_text(color="black", size=10),
                           legend.text = element_text(size=7))
   
    # TO DO: add number to legends
    png(paste0(output_dir,datatag,'_',de,'_signif_genes_cancer_ref.png'), height = 2*400, width=2*400,res = 2*72)
    print(p_ref)
    dev.off()
    
    s <- as.data.frame(summary(as.factor(signif_genes$Gene_Type)))
    # summary(-log10(signif_genes$pvalue))
    # summary(signif_genes$padj)
    colnames(s)[1] <- 'nb_genes'
    s$pct_genes <- round(s$nb_genes/colSums(s), 3) * 100
    colnames(s)[which(colnames(s) == "pct_genes")] <- paste0(datatag,'-',de)
    colnames(s)[which(colnames(s) == "nb_genes")] <- paste0(datatag,'-',de,'_nbgenes')
    c <- c + 1
    genes_summary_stat[[c]] <- s
    s$gene_type <- rownames(s)
    write.csv(s,paste0(output_dir,datatag,'_',de,'.csv'), quote = F, row.names = F)
    write.csv(signif_genes, paste0(output_dir,datatag,'_',de,'_signif_genes.csv'), quote = F, row.names = F)
    # write.csv(signif_genes, paste0(output_dir,datatag,'_',de,'_all_genes.csv'), quote = F, row.names = F)
    
    genes_summary[[paste0(datatag,'_',de,'_signif_genes')]] <- signif_genes
    genes_summary[[paste0(datatag,'_',de,'_stat')]] <- s
    
  }
  stat <- do.call(cbind, genes_summary_stat)
  # View(head(stat))
  stat$gene_type <- rownames(stat)
  write.csv(stat, paste0(output_dir,datatag,'_summary.csv'), quote = F, row.names = F)
  
  saveRDS(genes_summary,paste0(output_dir, datatag,'_',de,'_summary_statistic.rds'))
  
} 


get_gene_id <- function(genes, cores_use=8) {
  
  labels <- mclapply(strsplit(genes, "_"), function(x) {
    return(x[1])
  }, mc.cores = cores_use)
  return(as.character(labels))
}

filter_data <- function(gene_expression_data, gex_var_quantile=0.5){
  var_log_counts <- colVars(as.matrix(log2(gene_expression_data+1)))
  
  var_quantile <- quantile(var_log_counts, probs = as.numeric(gex_var_quantile),na.rm=T)
  # # rowData(sce)$total_counts > 0 & 
  print(var_quantile)
  genes <- colnames(gene_expression_data)
  genes <- genes[var_log_counts > var_quantile]
  print(paste0("Nb common genes after twice filtering: ",length(genes)))
  gene_expression_data <- gene_expression_data[,colnames(gene_expression_data) %in% genes]
  return(t(gene_expression_data))
}

compute_normalized_mean_gene <- function(normalizedCount, obs_clones, 
                                         colData, outputPrefix,
                                         cores_use=5){
  ls_val <- list()
  tmp <- rep(0, times = nrow(normalizedCount))
  for(c in obs_clones){
    ls_val[[c]] <- tmp
  }
  normalizedMean <- data.frame(matrix(unlist(ls_val), ncol=length(ls_val), byrow=F),
                               row.names = rownames(normalizedCount))
  colnames(normalizedMean) <- names(ls_val)
  dim(normalizedMean)
  
  cells_clones <- list()
  for(cn in obs_clones){
    cells_use <- rownames(colData[colData$clone==cn,])
    # cells_use <- rownames(colData[colData[,condition]==cn,])
    cells_clones[[cn]] <- cells_use
  }
  
  gene_mean <- NULL
  gene_mean <- mclapply(rownames(normalizedMean), function(g) {
    g_tmp <- normalizedMean[g,]
    for(cn in obs_clones){
      SUM = 0
      for(c in cells_clones[[cn]]){
        SUM = SUM + normalizedCount[g, c]
      }
      g_tmp[g,cn] <- SUM / length(cells_clones[[cn]])
    }
    g_tmp
  }, mc.cores = cores_use) %>% bind_rows()
  
  saveRDS(gene_mean, paste0(outputPrefix, "normalized_mean.rds"))
  return(gene_mean)
}




sce_cbind_func <- function(sce_list, genes_df, cut_off_overall = 0.05, exprs = c("counts", "logcounts","normcounts"), 
                           colData_names = NULL, meta_data=NULL) {
  n_batch <- length(sce_list)
  # method = "intersect"
  
  assay_list <- list()
  for (i in seq_len(length(exprs))) {
    assay_list[[i]] <- do.call(cbind, lapply(sce_list, 
                                             function(y) assay(y, exprs[i])))
  }
  names(assay_list) <- exprs
  colData_list <- do.call(DelayedArray::rbind, 
                          lapply(sce_list, function(y) colData(y)[, colData_names, drop = FALSE]))
  sce_combine <- SingleCellExperiment::SingleCellExperiment(assay = assay_list, 
                                                            colData = colData_list)
  
  rowData(sce_combine) <- genes_df
  # if (is.null(batches)) {
  #   batches <- paste("batch", seq_len(n_batch))
  # }
  # meta_data$batch_info <- gsub('^CHIP','C_',meta_data$batch_info)
  # print(meta_data$batch_info[1:4])
  
  # meta_data$pdxid <- gsub('^X0847_','',meta_data$pdxid)
  # print(meta_data$pdxid[1:4])
  
  # sce_combine$mouse_id <- meta_data$mouse_id
  # sce_combine$pdxid <- meta_data$pdxid
  # sce_combine$passage <- meta_data$passage
  # sce_combine$treatmentSt <- meta_data$treatmentSt
  # sce_combine$cell_cycle_phases <- meta_data$cphases
  # sce_combine$library_label <- meta_data$library_label
  # sce_combine$batch_info <- meta_data$batch_info
  # sce_combine$batch <- rep(batches, unlist(lapply(sce_list, ncol)))  #id
  # sce_combine$libid <- rep(libid, unlist(lapply(sce_list, ncol)))
  # sce_combine$treatment_status <- rep(treatment_status, unlist(lapply(sce_list, ncol)))
  
  if(cut_off_overall > 0){
    zero_cbind <- DelayedArray::rowMeans(assay(sce_combine, exprs[1]) == 0)
    sce_combine <- sce_combine[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
  }
  print(paste0("Dim sce combine: ",dim(sce_combine)))
  print(paste0("sce combine assay name: ",assayNames(sce_combine)))
  return(sce_combine)
}

load_clone_labels <- function(sce, input_dir, datatag, sample_fn, save_dir){
  print("Load sample metadata file: ")
  sample_df <- read.csv(sample_fn, check.names=F, stringsAsFactors=F)
  # View(sample_df)
  print(dim(sample_df))
  print(length(unique(sample_df$mouse_id)))
  
  clones_list <- list()
  c <- 0
  for(s in unique(sample_df$mouse_id)){
    cls_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/', s,'.csv')
    if(file.exists(cls_fn)){
      sclone_df <- read.csv(cls_fn, check.names = F, stringsAsFactors = F)
      print(s)
      print(dim(sclone_df))
      if(nrow(sclone_df)>0){
        c <- c + 1
        clones_list[[c]] <- sclone_df
      }
    }else{
      print(paste0('*** DEBUG: do not exist sample: ',s))
    }
  }  
  
  clone_df <- do.call(rbind, clones_list)
  print(dim(clone_df))

  if(sum(c('Sample','id') %in% colnames(clone_df))==2){
    clone_df$Sample <- gsub('(.cache/)','',clone_df$Sample)
    clone_df$Sample <- gsub('(/filtered_feature_bc_matrix)','',clone_df$Sample)
    colnames(clone_df)[which(names(clone_df) == "Sample")] <- "library_id"
    colnames(clone_df)[which(names(clone_df) == "id")] <- "mouse_id"
    print(head(clone_df))
    print(summary(as.factor(clone_df$library_id)))
  }
  
  write.csv(clone_df, paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_clones.csv'), row.names = F, quote = F)
  
  cells_use <- intersect(colnames(sce), clone_df$cell_id)
  print(length(cells_use))
  if(length(cells_use)>0){
    sce <- sce[,cells_use]
    rownames(clone_df) <- clone_df$cell_id
    # sce$library_id <- clone_df[cells_use,'library_id']
    sce$clone <- clone_df[cells_use,'clone']
    print(length(clone_df[cells_use,'library_id']))
  }
  
  data <- counts(sce)
  colData <- clone_df[clone_df$cell_id %in% colnames(data),]
  print(dim(data))
  print(dim(colData))
  write.csv(colData, paste0(save_dir,'colData.csv'), row.names = F, quote = F)
  saveRDS(data, file = paste0(save_dir,'data.rds'))
  print(dim(sce))
  saveRDS(sce, paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce_clones.rds'))
  
}

load_sce_data <- function(input_dir, save_dir, datatag, sample_fn, cutoff_zero_thrs=0.05){
  sample_df <- read.csv(sample_fn, check.names=F, stringsAsFactors=F)
  # View(sample_df)
  print(dim(sample_df))
  print(length(unique(sample_df$mouse_id)))
  sce_list <- list()
  sce_fls <- list()
  c <- 0
  
  for(s in unique(sample_df$mouse_id)){
    sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/', s,'.rdata')
    if(file.exists(sce_fn)){
      sce_fls[[s]] <- sce_fn
    }else{
      stop(paste0('Double check input files, do not exist sample: ',s))
    }
  }
  for(f in sce_fls){
    sce <- readRDS(f)
    print(f)
    print(dim(sce))
    if(nrow(sce)>0){
      c <- c + 1
      sce_list[c] <- sce
    }
  }
  
  # From 33k to >10k genes, remove cells with zeroes expression values in 95% cells
  sce_combine <- sce_cbind_func(sce_list, rowData(sce), 
                                cut_off_overall = cutoff_zero_thrs, 
                                exprs = c("counts", "logcounts"), 
                                colData_names = colnames(colData(sce)), meta_data=NULL)
  
  print(dim(rowData(sce_combine)))
  print(dim(sce_combine))
  print(dim(colData(sce_combine)))
  # View(head(rowData(sce_combine)))
  print(rowData(sce_combine)$Symbol[1:3])
  library(stringr)
  mito_genes <- str_detect(rowData(sce_combine)$Symbol, "^MT\\-")
  print(paste0('Nb mito genes: ',sum(mito_genes==TRUE)))
  
  ribo_genes <- str_detect(rowData(sce_combine)$Symbol, "^RP(L|S)")  # or ^RP[L|S]
  print(paste0('Nb ribo genes: ',sum(ribo_genes==TRUE)))
  
  sce_combine <- sce_combine[(!mito_genes) & (!ribo_genes), ]
  print(dim(sce_combine))
  if(!is.null(rowData(sce_combine)$ID)){
    rownames(sce_combine) <- rowData(sce_combine)$ID
  }else{
    print('Double check input ens gene id in ID column')
  }
  
  
  # View(head(colData(sce_combine)$listData))     
  
  # colnames(sce_combine)[1:3]
  # dim(sce)
  # sce$library_ids[1:3]
  # sce_combine$Sample[1:3]
  # sce$series[1:3]
  # sce_combine$sample[1:3]
  # sce$timepoint[1:3]
  # sce_combine$Barcode[1:3]
  
  # Process library_ids
  if(!is.null(sce_combine$Sample)){
    library_ids <- sce_combine$Sample
    library_ids <- str_replace_all(library_ids, "(.cache/)", "")
    library_ids <- str_replace_all(library_ids, "(/filtered_feature_bc_matrix)", "")
    print(summary(as.factor(library_ids)))
    sce_combine$library_ids <- library_ids
  }

  
  colnames(sce_combine) <- paste0(sce_combine$sample,'_',sce_combine$Barcode)
  print(colnames(sce_combine)[1:3])
  
  # library(parallel)
  ts <- mclapply(strsplit(sce_combine$series, "-"), function(x) {
    return(x[3])
  }, mc.cores = 2)
  sce_combine$treatmentSt <- as.character(ts)
  saveRDS(sce_combine, paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce.rds'))
  
  cell_ids = paste0(sce_combine$library_ids,'_',sce_combine$Barcode)
  meta_df <- data.frame(cell_id=cell_ids, cell_sid=colnames(sce_combine), 
                        library_id=library_ids,
                        sample_id=sce_combine$sample,
                        treatmentSt=sce_combine$treatmentSt,
                        time_point=sce_combine$timepoint)
  
  
  write.csv(meta_df, file = paste0(save_dir,'filtered_cells_v62.csv'), row.names = F, quote = F)
  print(dim(meta_df))
}

plot_stack_barplot <- function(cls_df, colorcode, xlabel, ylabel, plottitle, output_dir, tag='', #axis_x=F, 
                               fa='gene_type', xa='PDX', ya='Fitness_genes_pct', yl=c(0,37), lg=F){
  # df_cumsum <- ddply(cls_df, "timepoint",
  #                    transform, 
  #                    label_ypos=cumsum(Freq) - 0.5*Freq + 1)
  # # 
  if(!lg){
    # xa <- element_text(size=7, hjust = 0.5, family=my_font, angle = 90)
    lg_pos <- "none"
  }else{
    # xa <- element_blank()
    lg_pos <- "bottom"
  }
  # 
  cls_df$pct <- round(cls_df[,ya],1)
  p <- ggplot2::ggplot(cls_df, ggplot2::aes_string(fill=fa, y=ya, x=xa)) + 
    ggplot2::geom_bar(position="dodge", stat="identity", width = 0.8) + 
    ggplot2::geom_text(ggplot2::aes(label=pct), size=2.5, color="black", position=position_dodge(width=0.9), vjust=-0.1)+ #vjust=-0.3, 
    ggplot2::facet_grid(. ~ PDX, scales="free", space='free') +  #, scales="free_x"
    # geom_text(aes(y=label_ypos, label=Freq), vjust=1.6, 
    #           color="white", size=3.5) +
    ggplot2::scale_fill_manual(values = colorcode)
  p <- p + ggplot2::labs(x=xlabel, y=ylabel, title=plottitle)
  # p <- p + facet_wrap(~timepoint)
  my_font <- "Helvetica"
  thesis_theme <- ggplot2::theme(
    text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
    axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),  
    axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(color="black",size=6, hjust = 0.5, family=my_font, angle = 90),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(color="black",size=7, hjust = 0.5, family=my_font),
    plot.title = element_text(color="black",size=10, face="bold", hjust=0, family=my_font),
    legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font), 
    legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
    strip.text.x = element_text(color="black",size=9, family=my_font),
    strip.text.y = element_text(color="black",size=9, family=my_font),
    legend.spacing.x = unit(0.1, 'mm'),
    legend.spacing.y = unit(0.1, 'mm'),
    legend.key.height=unit(1,"line"),
    legend.position = lg_pos,
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  
  p <- p + thesis_theme
  # p <- p + ggplot2::theme(legend.title = ggplot2::element_text(color="black", size=11),
  #                legend.text = ggplot2::element_text(color="black", size=10),
  #                plot.title = ggplot2::element_text(color="black", size=14, hjust = 0.5, face='bold'),
  #                # legend.position = "none", 
  #                # axis.line = element_blank(), 
  #                panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
  #                panel.border = ggplot2::element_blank(),
  #                # axis.text.x =ggplot2::element_blank(),
  #                axis.text.x = ggplot2::element_text(color="black", size=10, hjust=0.5, angle = 90),
  #                axis.text.y = ggplot2::element_text(color="black", size=8, hjust = 0.5),
  #                axis.title.y = ggplot2::element_text(color="black", size=10, hjust = 0.5),
  #                axis.title.x = ggplot2::element_text(color="black", size=12, hjust = 0.5, face='bold')
  #                # axis.title.x = ggplot2::element_blank()
  #                # axis.ticks = element_blank()
  # )
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  
  saveRDS(p, paste0(output_dir,tag,"_stat.rds"))
  png(paste0(output_dir,tag,"_stat.png"), height = 2*380, width=2*720, res = 2*72)
  print(p)
  dev.off()
  
  ggsave(file=paste0(output_dir,tag,"_stat.pdf"), plot=p, height=3.8, width=6.5, useDingbats=F)  #, dpi=2*72
  return(p)
  
}



get_gene_type_stat_testing <- function(signif_genes_df, datatag, 
                                       output_dir, outlier_FC_thrs=3)
{
  # dir.create(output_dir, showWarnings = F)
  output_dir <- paste0(output_dir,datatag,'_result/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  reference_genes_df <- read.csv(paste0(ref_dif,'Symbol_ensembl.csv'), 
                                 stringsAsFactors=F, check.names = F)
  
  cancer_ref_genes_df <- read.csv(paste0(ref_dif,'Behan_CFgenes.csv'), 
                                  stringsAsFactors=F, check.names = F)
  dim(cancer_ref_genes_df)
  colnames(cancer_ref_genes_df)[which(colnames(cancer_ref_genes_df) == "ADAM PanCancer Core-Fitness genes")] <- "PanCancer_Fitness_genes"
  
  
  signif_genes_df <- signif_genes_df %>% left_join(reference_genes_df, by = c("gene_id"="Ensembl"))
  signif_genes_df <- signif_genes_df %>% left_join(cancer_ref_genes_df, by = c("gene_id"="ensemble_id"))
  signif_genes_df$de_analysis <- paste0(signif_genes_df$sample1,'_vs_',signif_genes_df$sample2)
  
  de_pairs <- unique(signif_genes_df$de_analysis)
  genes_summary_stat <- list()
  genes_summary <- list()
  c <- 0
  for(de in de_pairs){
    signif_genes <- signif_genes_df[signif_genes_df$de_analysis==de,]
    print(dim(signif_genes))
    s1 <- as.character(signif_genes$sample1[1])
    s2 <-  as.character(signif_genes$sample2[1])
    print(paste0('clone ',s1, ' versus clone ',s2, ' nb significant genes: ',nrow(signif_genes)))
    # print(sum(colnames(cnv_mat) %in% c(s1, s2))==2)
    # cnv_mat_tmp <- cnv_mat[,c(s1, s2)]
    # print(dim(cnv_mat_tmp))
    # 
    # for(g in rownames(cnv_mat_tmp)){
    #   if(cnv_mat_tmp[g,s1]-cnv_mat_tmp[g,s2]>0){
    #     cnv_mat_tmp[g,'classified_gene_dlp'] <- 'Increase'
    #   }else if(cnv_mat_tmp[g,s1]-cnv_mat_tmp[g,s2]<0){
    #     cnv_mat_tmp[g,'classified_gene_dlp'] <- 'Decrease'
    #   }else{
    #     cnv_mat_tmp[g,'classified_gene_dlp'] <- 'No_variance'
    #   }
    # }  
    # print(summary(as.factor(cnv_mat_tmp$classified_gene_dlp)))
    # cnv_mat_tmp <- cnv_mat_tmp[cnv_mat_tmp$classified_gene_dlp != 'No_variance',]
    # print(dim(cnv_mat_tmp))
    # print(paste0('Nb of variance dlp genes: ',nrow(cnv_mat_tmp)))
    # cnv_mat_tmp$gene_id <- rownames(cnv_mat_tmp)
    # signif_genes <- signif_genes %>% left_join(cnv_mat_tmp, by = c("gene_id"))
    # in_cis <- intersect(signif_genes$gene_id, cnv_mat_tmp$gene_id)
    # print(paste0('Nb of variance in-cis genes: ',length(in_cis)))
    
    rownames(signif_genes) <- signif_genes$gene_id
    signif_genes$classified_gene <- 'in_trans'
    # signif_genes[in_cis,'classified_gene'] <- 'in_cis'
    print(summary(signif_genes$log2FoldChange))
    signif_genes_outliers <- signif_genes[abs(signif_genes$log2FoldChange)>outlier_FC_thrs,]
    write.csv(signif_genes_outliers, paste0(output_dir,datatag,'_',de,'_signif_genes_outliers.csv'), quote = F, row.names = F)
    signif_genes <- signif_genes[abs(signif_genes$log2FoldChange) <= outlier_FC_thrs,]
    
    
    signif_genes$classified_gene_rna <- ifelse(signif_genes$log2FoldChange>0,'UpRegulated',
                                               ifelse(signif_genes$log2FoldChange<0,'DownRegulated','Invalid'))
    # signif_genes$Gene_Type <- ifelse(!is.na(signif_genes$classified_gene_dlp),
    #                                  paste0('In_cis_',signif_genes$classified_gene_dlp, '_',signif_genes$classified_gene_rna),
    #                                  paste0('In_trans_',signif_genes$classified_gene_rna))
    
    signif_genes$Gene_Type <- signif_genes$classified_gene_rna
    signif_genes$Gene_Type <- as.factor(signif_genes$Gene_Type)
    p <- ggplot(signif_genes) +
      geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour=Gene_Type),size = 1.5) +  #aes(colour = factor(is_fitness_gene)), 
      scale_color_discrete(labels = paste0(levels(signif_genes$Gene_Type),': ',table(signif_genes$Gene_Type)))
    
    p <- p + labs(title=paste0(datatag,' ',de))
    p <- p + theme(plot.title = element_text(color="black", size=12, hjust = 0.5),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   axis.text.y = element_text(color="black", size=10),
                   axis.text.x = element_text(color="black", size=10),
                   axis.title = element_text(color="black", size=10),
                   legend.text = element_text(size=7))
    # p <- ggplot(signif_genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
    #   geom_point(aes(colour = Gene_Type), size = 2)
    # 
    png(paste0(output_dir,datatag,'_',de,'_signif_genes.png'), height = 2*400, width=2*500,res = 2*72)
    print(p)
    dev.off()
    
    signif_genes$is_fitness_gene <- as.factor(!is.na(signif_genes$PanCancer_Fitness_genes))
    print(paste0('Comparison: ',de,"  Number of fitness genes: "))
    print(summary(as.factor(signif_genes$is_fitness_gene)))
    
    p_ref <- ggplot(signif_genes) +
      geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour=is_fitness_gene),size = 1.5) +  #aes(colour = factor(is_fitness_gene)), 
      scale_color_discrete(labels = paste0(levels(signif_genes$is_fitness_gene),': ', table(signif_genes$is_fitness_gene)))
    p_ref <- p_ref + labs(title=paste0(datatag,' ',de))
    p_ref <- p_ref + theme(plot.title = element_text(color="black", size=12, hjust = 0.5),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                           axis.text.y = element_text(color="black", size=10),
                           axis.text.x = element_text(color="black", size=10),
                           axis.title = element_text(color="black", size=10),
                           legend.text = element_text(size=7))
    
    # TO DO: add number to legends
    png(paste0(output_dir,datatag,'_',de,'_signif_genes_cancer_ref.png'), height = 2*400, width=2*400,res = 2*72)
    print(p_ref)
    dev.off()
    
    s <- as.data.frame(summary(as.factor(signif_genes$Gene_Type)))
    # summary(-log10(signif_genes$pvalue))
    # summary(signif_genes$padj)
    colnames(s)[1] <- 'nb_genes'
    s$pct_genes <- round(s$nb_genes/colSums(s), 3) * 100
    colnames(s)[which(colnames(s) == "pct_genes")] <- paste0(datatag,'-',de)
    colnames(s)[which(colnames(s) == "nb_genes")] <- paste0(datatag,'-',de,'_nbgenes')
    c <- c + 1
    genes_summary_stat[[c]] <- s
    s$gene_type <- rownames(s)
    write.csv(s,paste0(output_dir,datatag,'_',de,'.csv'), quote = F, row.names = F)
    write.csv(signif_genes, paste0(output_dir,datatag,'_',de,'_signif_genes.csv'), quote = F, row.names = F)
    genes_summary[[paste0(datatag,'_',de,'_signif_genes')]] <- signif_genes
    genes_summary[[paste0(datatag,'_',de,'_stat')]] <- s
    
  }
  stat <- do.call(cbind, genes_summary_stat)
  # View(head(stat))
  stat$gene_type <- rownames(stat)
  write.csv(stat, paste0(output_dir,datatag,'_summary.csv'), quote = F, row.names = F)
  
  saveRDS(genes_summary,paste0(output_dir, datatag,'_',de,'_summary_statistic.rds'))
  
} 

get_stat_logFC <- function(pair_groups, input_dir, save_dir=NULL){
  ## To evaluate cis, trans genes, using this threshold
  FDR_cutoff <- 0.01
  minLogFC <- 0.5
  pValueThrs <- 0.05
  # Load DE comparisons by file header
  # Get gene type, and statistics
  if(!dir.exists(save_dir)){
    dir.create(save_dir, recursive=T)
  }
  rownames(pair_groups) <- pair_groups$file_header
  pair_groups$ks_stat <- 1 # initial stat value
  stat <- tibble::tibble()
  combined_DE <- tibble::tibble()
  for(idx in pair_groups$file_header){
    fn <- paste0(input_dir, idx,'/signif_genes.csv')
    if(file.exists(fn)){
      de <- data.table::fread(fn) %>% as.data.frame()
      de <- de %>%
        dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs & classified_gene!='unmapped')
      # print(dim(de))
      de_cis <- de %>%
        dplyr::filter(classified_gene=='in_cis')
      de_trans <- de %>%
        dplyr::filter(classified_gene=='in_trans')
      # "Two-sample Kolmogorov-Smirnov test"
      res <- ks.test(abs(de_cis$logFC), abs(de_trans$logFC), alternative="l")
      pair_groups[idx,'ks_stat'] <- res$p.value
      de <- de %>%
        dplyr::select(ensembl_gene_id, logFC, classified_gene, Gene_Type)
      de$file_header <- idx
      combined_DE <- dplyr::bind_rows(combined_DE, de)
    }
  }
  # View(pair_groups)
  # dim(pair_groups)
  pair_groups$ks_signf <- ifelse(pair_groups$ks_stat<0.05,'Signif','No')
  print(sum(pair_groups$ks_stat<0.05))
  data.table::fwrite(pair_groups, paste0(save_dir, 'pair_groups_ks_test.csv.gz'))  
  data.table::fwrite(combined_DE, paste0(save_dir, 'combined_DE_logFC.csv.gz'))  
  for(tag in unique(pair_groups$datatag)){
    out_tmp <- pair_groups %>%
      dplyr::filter(datatag==tag) %>%
      dplyr::arrange(order) %>%
      dplyr::select(datatag, labels_detailed, ks_signf)
    data.table::fwrite(out_tmp, paste0(save_dir, tag, '_ks_test.csv.gz'))  
  }
  
  return(list(pair_groups=pair_groups, combined_DE=as.data.frame(combined_DE)))  
}

## plottitle: 'UnRx', 'Rx'
# df <- df1
# pc1_label <- viz_labels(df, lb_use='clone1',plottitle='UnRx', color_scheme='predefined_clone_col',tag='SA609')
# pc2_label <- viz_labels(df, lb_use='clone2',plottitle='Rx', color_scheme='predefined_clone_col',tag='SA609')
# p_passage <- viz_labels(df, lb_use='passage',plottitle='Passage', color_scheme='gradient')

viz_labels <- function(df, lb_use='clone1',plottitle='Rx',color_scheme=NULL,tag='SA'){
  # library(ggplot2)
  if(tag=='SA609'){
    df$clone1 <- ifelse(df$clone1=='R','A',df$clone1)
    df$clone2 <- ifelse(df$clone2=='R','A',df$clone2)
  }
  df <- df %>%
    dplyr::group_by(file_header, order, !!sym(lb_use)) %>%
    dplyr::summarise(nb_gt = sum(nb_genes)) ## just to extract unique clone 1 and file header info
  # dim(df)
  df$xst <- 1
  my_font <- "Helvetica"
  # df <- data.frame(passage=c('X2','X3','X4','X5','X6'), xst=c(rep(1,5)))
  # df
  df <- df %>%
    dplyr::arrange(-order)
  
  
  ps <- as.character(unlist(unique(df[,lb_use])))
  # ps <- c('X2','X4','X5','X3','X6')
  ps <- gtools::mixedsort(ps)
  if(color_scheme=='gradient'){
    colfunc <- colorRampPalette(c("#767676","black"))  
    cols_use <- colfunc(length(ps))  
    names(cols_use) <- ps
  }else if(color_scheme=='predefined_clone_col'){
    # loading from file, predefined color scheme for each clone
    base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/'
    color_df <- data.table::fread(paste0(base_dir,'umap_figs/colorcode_total_v4.csv'))
    
    tmp <- color_df %>%
      dplyr::filter(datatag==tag)
    cols_use <- tmp$colour
    names(cols_use) <- tmp$clone_id
    # head(color_df)
  }else{
    colfunc <- colorRampPalette(c("black","#767676"))  ## to do for untreated color here
    cols_use <- colfunc(length(ps))  
    names(cols_use) <- ps
  }
  # print(cols_use)
  
  df$file_header <- factor(df$file_header, levels = df$file_header)
  p <- ggplot(df, aes_string(x='xst', y='file_header', label = lb_use)) + 
    geom_label(aes_string(fill = lb_use), colour = "white", fontface = "bold") + 
    scale_fill_manual(values = cols_use) +
    # scale_fill_gradientn(colours = cols_use) + 
    theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          # panel.border = element_blank(),
          # axis.title = element_blank(),
          # axis.text = element_blank(),
          # axis.ticks = element_blank(),
          legend.position = "none",
          plot.title = element_text(size=12, hjust = 0.5, family=my_font),
          # panel.background = element_rect(fill = 'white', colour = 'white')
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "null"),
          panel.spacing = unit(c(0, 0, 0, 0), "null"), #= unit(2, "lines")
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          axis.ticks.length = unit(0, "null"),
          # axis.ticks.margin = unit(0, "null"),
          # legend.margin = unit(0, "null")
          )+
    labs(title=plottitle, x=NULL, y=NULL)#, x=' '
  # p
  return(p)
}



get_gene_type_stat_v3 <- function(pair_groups, input_dir, save_dir)
{
  ## To evaluate cis, trans genes, using this threshold
  FDR_cutoff <- 0.01
  minLogFC <- 0.5
  pValueThrs <- 0.05
  # Load DE comparisons by file header
  # Get gene type, and statistics
  if(!dir.exists(save_dir)){
    dir.create(save_dir, recursive=T)
  }
  stat <- tibble::tibble()
  for(idx in pair_groups$file_header){
    fn <- paste0(input_dir, idx,'/signif_genes.csv')
    if(file.exists(fn)){
      de <- data.table::fread(fn) %>% as.data.frame()
      de <- de %>%
        dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
      print(dim(de))
      # nb_genes,pct_genes,datatag,desc,result_fn,de_analysis,gene_type
      stat_de <- de %>%
        dplyr::group_by(Gene_Type) %>%
        dplyr::summarise(nb_genes=n(),
                         pct_genes=round(n()/dim(de)[1]*100,2)#,
                         # avgLogFC=round(mean(abs(logFC)),2),
                         # sdLogFC=round(sd(logFC),2)
                         ) %>%
        dplyr::rename(gene_type=Gene_Type) %>%
        dplyr::mutate(file_header=idx)
      
      stat <- dplyr::bind_rows(stat, stat_de)
      
    }else{
      print('Do not exist the comparison: ')
      print(idx)
    }
  }
  
  print(dim(stat))
  stat <- stat %>% left_join(pair_groups, by='file_header')  
  data.table::fwrite(stat, paste0(input_dir, 'wholedata_summary.csv'))
  return(stat)
}  

get_gene_type_stat_manuscript <- function(pair_groups, input_dir, save_dir)
{
  ## To evaluate cis, trans genes, using this threshold
  FDR_cutoff <- 0.01
  minLogFC <- 0.5
  pValueThrs <- 0.05
  # Load DE comparisons by file header
  # Get gene type, and statistics
  if(!dir.exists(save_dir)){
    dir.create(save_dir, recursive=T)
  }
  pos_gt <- c('In_cis_Decrease_DownRegulated','In_cis_Increase_UpRegulated')
  neg_gt <- c('In_cis_Decrease_UpRegulated','In_cis_Increase_DownRegulated')
  trans_gt <- c('In_trans_UpRegulated','In_trans_DownRegulated')
  stat <- tibble::tibble()
  for(idx in pair_groups$file_header){
    fn <- paste0(input_dir, idx,'/signif_genes.csv')
    if(file.exists(fn)){
      de <- data.table::fread(fn) %>% as.data.frame()
      de <- de %>%
        dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
      print(dim(de))
      ## nb_genes,pct_genes,datatag,desc,result_fn,de_analysis,gene_type
      # stat_de <- de %>%
      #   mutate(gt_stat = case_when(Gene_Type %in% pos_gt ~ 'cis_pos',
      #                              Gene_Type %in% neg_gt ~ 'cis_neg',
      #                              Gene_Type %in% trans_gt ~ 'trans',
      #                                           TRUE ~ 'unmapped'))
      stat_de <- de
      stat_de$gt_stat <- stat_de$Gene_Type
      # stat_de <- de %>%
      #   mutate(gt_stat = case_when(Gene_Type %in% c(pos_gt,neg_gt) ~ 'cis',
      #                              Gene_Type %in% trans_gt ~ 'trans',
      #                              TRUE ~ 'unmapped'))
      stat_de <- stat_de %>%
        dplyr::group_by(gt_stat) %>%
        dplyr::summarise(nb_genes=n(),
                         pct_genes=round(n()/dim(de)[1]*100,2),
                         avgLogFC=round(mean(abs(logFC)),2),
                         sdLogFC=round(sd(logFC),2)) %>%
        dplyr::rename(gene_type=gt_stat) %>%
        dplyr::mutate(file_header=idx)
      
      stat <- dplyr::bind_rows(stat, stat_de)
      
    }else{
      print('Do not exist the comparison: ')
      print(idx)
    }
  }
  
  print(dim(stat))
  stat <- stat %>% left_join(pair_groups, by='file_header')  
  # data.table::fwrite(stat, paste0(input_dir, 'wholedata_summary_manuscript.csv.gz'))
  return(stat)
}  

get_corrected_gene_symbol <- function(ref_df, reference_set_name='', obs_genes=NULL){
  ref_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/cancer_reference_genes/'
  database_genes <- data.table::fread(paste0(ref_dir, 'custom_ens_genes_symbols.txt'), sep='\t') %>% as.data.frame()
  # dim(database_genes)
  database_genes$corrected_gene_symb <- database_genes$`Approved symbol`
  database_genes$previous_symbols <- database_genes$`Previous symbols`
  database_genes$alias <- database_genes$`Alias symbols`
  
   meta_genes <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/Symbol_ensembl.csv') %>% as.data.frame()
  
  if(is.null(obs_genes)){
    obs_genes <- ref_df %>%
      dplyr::filter(!gene_symb %in% meta_genes$Symbol) %>%
      dplyr::pull(gene_symb)
    obs_genes <- unique(obs_genes)
  }
  # t <- annotables::grch38 %>% 
  #   dplyr::select(gene_symb = symbol,) #%>% 
  # # inner_join(cf, by='gene_symb')
  # dim(t)
  # t1 <- t$gene_symb[grepl('RAB8A',t$gene_symb)]
  # head(t1)
  anno <- tibble::tibble()
  for(t in obs_genes){
    # print(t)
    re <- grepl(toupper(t),database_genes$alias)
    if(sum(re==TRUE)>0){
      tmp <- database_genes[re,c('corrected_gene_symb','previous_symbols','alias')]
      tmp$gene_symb <- t
      anno <- dplyr::bind_rows(anno, tmp)
    }
    re1 <- grepl(toupper(t),database_genes$previous_symbols)
    if(sum(re1==TRUE)>0){
      tmp1 <- database_genes[re1,c('corrected_gene_symb','previous_symbols','alias')]
      tmp1$gene_symb <- t
      anno <- dplyr::bind_rows(anno, tmp1)
    }
  }
  out <- tibble::tibble()
  for(i in seq(1:dim(anno)[1])){
    # gs <- unlist(strsplit(anno$previous_symbols[i],', '))
    gs <- unlist(strsplit(anno$alias[i],', '))
    if(length(gs)>0){
      tmp <- tibble::tibble(corrected_gene_symb=rep(anno$corrected_gene_symb[i],length(gs)),gene_symb=gs)
      out <- dplyr::bind_rows(out, tmp)  
    }
    
  }
  for(i in seq(1:dim(anno)[1])){
    # gs <- unlist(strsplit(anno$previous_symbols[i],', '))
    gs <- unlist(strsplit(anno$previous_symbols[i],', '))
    if(length(gs)>0){
      tmp <- tibble::tibble(corrected_gene_symb=rep(anno$corrected_gene_symb[i],length(gs)),gene_symb=gs)
      out <- dplyr::bind_rows(out, tmp)  
    }
    
  }
  print(dim(out))
  out <- out %>%
    dplyr::filter(gene_symb %in% toupper(obs_genes))
  # correct_mapping <- summary(as.factor(out$gene_symb))
  # correct_mapping <- correct_mapping[correct_mapping==1] # only one mapping between alias and correct gene name
  # out1 <- out %>%
  #   dplyr::filter(gene_symb %in% names(correct_mapping))
  # out2 <- out %>%
  #   dplyr::filter(!gene_symb %in% names(correct_mapping))
  # out1
  ref_df <- ref_df %>%
    dplyr::filter(!gene_symb %in% out$gene_symb)
  # dim(ref_df)
  out <- out %>%
    dplyr::select(corrected_gene_symb) %>%
    dplyr::rename(gene_symb=corrected_gene_symb)
  ref_df <- dplyr::bind_rows(ref_df, out)
  dim(ref_df)
  if(reference_set_name!=''){
    data.table::fwrite(ref_df, paste0(ref_dir, reference_set_name,'.csv'))  
  }
  
  return(ref_df)
  
}


get_proportion_cistrans_reference_set <- function(pair_groups, input_dir, save_dir)
{
  ## Load reference data first
  ref_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/cancer_reference_genes/'
  reference_sets=c('CoreFitness_corrected','cisplatin_resistance_corrected') #'BroadSanger','cosmic'
    
  cf <- data.table::fread(paste0(ref_dir, reference_sets[1],'.csv')) %>% as.data.frame()
  
  cf <- cf[!duplicated(cf$gene_symb),, drop=F]
  cr <- data.table::fread(paste0(ref_dir, reference_sets[2],'.csv')) %>% as.data.frame()
  
  cr <- cr[!duplicated(cr$gene_symb),, drop=F]
  meta_genes <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/Symbol_ensembl.csv') %>% as.data.frame()
  
  cf <- cf %>%
    dplyr::filter(gene_symb %in% meta_genes$Symbol)
  
  cr <- cr %>%
    dplyr::filter(gene_symb %in% meta_genes$Symbol)
  
  print(dim(cf))
  print(dim(cr))
  # Correct genes symbols
  # cf <- get_corrected_gene_symbol(cf, 'CoreFitness_corrected')
  # cr <- get_corrected_gene_symbol(cr, 'cisplatin_resistance_corrected')
  
  ## To evaluate cis, trans genes, using this threshold
  FDR_cutoff <- 0.01
  minLogFC <- 0.5
  pValueThrs <- 0.05
  # Load DE comparisons by file header
  # Get gene type, and statistics
  if(!dir.exists(save_dir)){
    dir.create(save_dir, recursive=T)
  }
  stat_cf <- tibble::tibble()
  stat_cr <- tibble::tibble()
  
  for(idx in pair_groups$file_header){
    fn <- paste0(input_dir, idx,'/signif_genes.csv')
    if(file.exists(fn)){
      de <- data.table::fread(fn) %>% as.data.frame()
      de <- de %>%
        dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & 
                      PValue<pValueThrs & Gene_Type!='unmapped')
      
      print(dim(de))
      de <- de %>% 
        dplyr::mutate(Gene_Type = case_when(grepl('In_cis',Gene_Type) ~ 'In_cis',
                                     grepl('In_trans',Gene_Type) ~ 'In_trans'))
      
      de_cf <- de %>%
        dplyr::filter(gene_symbol %in% cf$gene_symb) %>%
        dplyr::group_by(Gene_Type) %>%
        dplyr::summarise(nb_genes=n(),
                         pct_genes=round(n()*100/dim(cf)[1],2)) %>%
        dplyr::rename(gene_type=Gene_Type) %>%
        dplyr::mutate(file_header=idx, ref_set='CF')
      stat_cf <- dplyr::bind_rows(stat_cf, de_cf)
      
      de_cr <- de %>%
        dplyr::filter(gene_symbol %in% cr$gene_symb) %>%
        dplyr::group_by(Gene_Type) %>%
        dplyr::summarise(nb_genes=n(),
                         pct_genes=round(n()*100/dim(cr)[1],2)) %>%
        dplyr::rename(gene_type=Gene_Type) %>%
        dplyr::mutate(file_header=idx, ref_set='CR')
      stat_cr <- dplyr::bind_rows(stat_cr, de_cr)
      
    }else{
      print('Do not exist the comparison: ')
      print(idx)
    }
  }
  
  # dim(stat)
  stat_cf <- stat_cf %>% left_join(pair_groups, by='file_header')  
  stat_cr <- stat_cr %>% left_join(pair_groups, by='file_header')  
  data.table::fwrite(stat_cf, paste0(input_dir, 'wholedata_summary_CF.csv'))
  data.table::fwrite(stat_cr, paste0(input_dir, 'wholedata_summary_CR.csv'))
  return(list(stat_cf=stat_cf, stat_cr=stat_cr))
}  

get_pathways_reference_set <- function(pair_groups, input_dir, save_dir)
{
  ## Load reference data first
  ref_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/cancer_reference_genes/'
  ref_genes_ls <- load_reference_genes_set(ref_dir, ref_fn=NULL)  
  print(names(ref_genes_ls))
  length(ref_genes_ls$CISPLATIN_RESISTANCE)
  length(ref_genes_ls$CORE_FITNESS)
  
  ## To compute pathway, using this threshold
  FDR_cutoff <- 0.01
  minLogFC <- 0.25
  pValueThrs <- 0.05
  # Load DE comparisons by file header
  # Get gene type, and statistics
  if(!dir.exists(save_dir)){
    dir.create(save_dir, recursive=T)
  }
  pathway_stat <- tibble::tibble()
  for(idx in pair_groups$file_header){
    fn <- paste0(input_dir, idx,'/signif_genes.csv')
    if(file.exists(fn)){
      de <- data.table::fread(fn) %>% as.data.frame()
      de <- de %>%
        dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
      
      print(dim(de))
      if(dim(de)[1]>0){
        # Deal with overlapping and NA gene symbols
        de <- de %>%
          dplyr::select(gene_symbol, logFC) %>%
          na.omit() %>%
          distinct() %>%
          group_by(gene_symbol) %>%
          summarize(logFC=mean(logFC))
        
        deg_stat <- de$logFC
        names(deg_stat) <- de$gene_symbol
        gsea_out <- fgsea::fgsea(pathways=ref_genes_ls, stats=deg_stat, nPermSimple = 10000) #, nperm=10000
        if(dim(gsea_out)[1]>0){
          gsea_out$important_genes <- ''
          for(i in rep(1:length(gsea_out$leadingEdge),1)){
            gsea_out$important_genes[i] <- as.character(paste(unlist(gsea_out$leadingEdge[i]), collapse = ','))
          }
          gsea_out <- gsea_out %>%
            dplyr::select(-leadingEdge)%>%
            as.data.frame()%>%
            dplyr::mutate(file_header=idx)
          
          pathway_stat <- dplyr::bind_rows(pathway_stat, gsea_out)  
          
        }
        
      }
    }else{
      print('Do not exist the comparison: ')
      print(idx)
    }
  }
  
  pathway_stat <- as.data.frame(pathway_stat)
  data.table::fwrite(pathway_stat, paste0(input_dir, 'wholedata_pathways_CF_CR_sets.csv'))
  return(pathway_stat)
}  


get_metadata <- function(){
  # library(stringr)
  # library(dplyr)
  id609 <- paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/')
  id1035 <- paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/')
  id535 <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/')
  res_fn_df <- list()
  c <- 0
  for(id in list(id609,id1035,id535)){
    logFC_files <- list.files(path=id, pattern = "_logfc_results.csv$", recursive = F)
    
    desc <- gsub('_logfc_results.csv','',logFC_files)
    desc <- lapply(strsplit(desc, "_"), function(x) {
      to_str <- ''
      for(i in rep(2:(length(x)-1), 1)){
        to_str <- paste0(to_str, x[i],'_')
      }
      to_str <- paste0(to_str, x[length(x)])
      return(to_str)
      # return(paste(x[2],x[3],x[4],x[5],x[6],sep='_'))
    })
    # desc[1] <- 'comp535c2-SA535X9XB03617_T_SA535X9XB03776_J'
    comp_id = lapply(strsplit(logFC_files, "_"), function(x) {
      return(x[1])
    })
    
    
    res_ls <- data.frame(result_fn=logFC_files,desc=as.character(desc),comp_id=as.character(comp_id))
    c <- c + 1
    res_fn_df[[c]] <- res_ls
  }
  fns <- do.call(rbind,res_fn_df)
  # View(fns)
  pair_groups <- read.csv(pair_groups_fn, check.names=F, stringsAsFactors=F)
  pair_groups$comp_id <- paste0(pair_groups$comp,rownames(pair_groups))
  pair_groups <- pair_groups %>% left_join(fns, by=c('comp_id'))
  pair_groups$datatag <- pair_groups$series
  # View(pair_groups)
  pair_groups <- pair_groups %>%
    dplyr::select(-comp)
  
  # View(pair_groups)
  write.csv(pair_groups,file = paste0(id609,'comparisons_drug_res.csv'), row.names = F, quote = F)
  write.csv(pair_groups,file = paste0(id1035,'comparisons_drug_res.csv'), row.names = F, quote = F)
  write.csv(pair_groups,file = paste0(id535,'comparisons_drug_res.csv'), row.names = F, quote = F)
  
}

update_metadata <- function(input_fn){
  
  cis_df <- read.csv(input_fn, check.names = F, stringsAsFactors = F)
  dim(cis_df)
  colnames(cis_df)
  
  ## Broad Sanger df do not have right format, need to change it here
  # cis_df <- cis_df %>%
  #       dplyr::select(pct_incis, pct_intrans, desc)
  # cis_df <- cis_df %>% pivot_longer(!desc, names_to = "gene_type", values_to = "gene_pct")
  # cis_df$gene_type <- ifelse(grepl("*incis",cis_df$gene_type),"in cis",'in trans')
  
  cis_df$gene_type <- ifelse(grepl("*incis",cis_df$desc),"in cis",'in trans')
  # viz_pancancer_genes(cancer_pct_df, output_dir, tag="totaldata")
  summary(as.factor(cis_df$gene_type))
  
  # cis_df$de_pair <- cis_df$desc
  # cancer_pct_df$de_pair <- gsub('_incis', '', grep("*incis",cancer_pct_df$de_pair, value=T))
  # cancer_pct_df$de_pair <- gsub('_intrans', '', grep("*intrans",cancer_pct_df$de_pair, value=T))
  cis_df$desc <- ifelse(grepl("*incis",cis_df$desc), gsub('_incis', '', cis_df$desc),
                           ifelse(grepl("*intrans",cis_df$desc), gsub('_intrans', '', cis_df$desc),cis_df$desc))
  
  pair_df <- read.csv('/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/incis_intrans/comparisons_drug_res.csv', stringsAsFactors=F, check.names=F)
  # head(pair_df)

  pair_df <- pair_df %>% 
    dplyr::select(title, DE_desc1, desc)
  
  ## Broad Sanger case
  # pair_df <- pair_df %>% 
  #   dplyr::rename(PDX=title)
  # colnames(pair_df)
  # rownames(pair_df) <- pair_df$desc
  
  cis_df <- cis_df %>% left_join(pair_df, by=c('desc'))
  unique(cis_df$PDX)
  cis_df$PDX <- ifelse(cis_df$PDX=='SA535_CISPLATIN','SA535:Cisplatin',
                      ifelse(cis_df$PDX=='SA535_CX5461','SA535:CX5461', cis_df$PDX))
  cis_df$PDX <- ifelse(cis_df$PDX=='SA535_cisplatin','SA535:Cisplatin',
                       ifelse(cis_df$PDX=='SA535_CX5461','SA535:CX5461', cis_df$PDX))
  cis_df$PDX <- gsub('_CISPLATIN','',cis_df$PDX)
  cis_df$PDX <- factor(cis_df$PDX, levels = unique(cis_df$PDX))
  cis_df$DE_desc1 <- gsub('M_N_O','MNO',cis_df$DE_desc1)
  cis_df$DE_desc1 <- gsub('S_T','ST',cis_df$DE_desc1)
  cis_df$de_pair <- cis_df$DE_desc1
  # cis_df$de_pair <- gsub('SA609_UTTT_R','SA609_UTTT_A',cis_df$de_pair)
  # cis_df$de_pair <- gsub('SA609_UTTTT_R','SA609_UTTTT_A',cis_df$de_pair)
  # cis_df$de_pair <- gsub('^(SA609_)|(SA1035_)|(SA535_)','',cis_df$de_pair)
  # level_pdx <- c('SA609_CISPLATIN','SA1035_CISPLATIN','SA535_CISPLATIN','SA535_CX5461')
  # cis_df$PDX <- factor(cis_df$PDX, levels=level_pdx)
  return(cis_df)
}




get_xlabel <- function(datatag, de=NULL, meta=NULL){
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  de_files <- c('SA609_summary.csv','SA1035_summary.csv','SA535_cisplatin_summary.csv','SA535_cx_summary.csv')
  names(de_files) <- c('SA609','SA1035','SA535_cisplatin','SA535_CX')
  
  meta_files <- c('SA609_10x.csv','SA1035_10x.csv','SA535_10x.csv','SA535_10x.csv')
  names(meta_files) <- c('SA609','SA1035','SA535_cisplatin','SA535_CX')
  ls_meta <- list()
  for(m in meta_files[1:3]){
    print(m)
    meta <- read.csv(paste0(base_dir,'rnaseq_v6/incis_intrans/',m), stringsAsFactors=F, check.names=F)
    ls_meta[[m]] <- meta
  }
  meta_df <- dplyr::bind_rows(ls_meta)
  class(meta_df)  
  # View(meta_df)
  colnames(meta_df)
  meta_df <- meta_df %>%
    dplyr::select(PDX, passage, mouse_id,library_id, treatmentSt, batch_info)
  write.csv(meta_df, paste0(base_dir,'rnaseq_v6/incis_intrans/total_metadata.csv'), quote=F, row.names = F)
  meta_df <- meta_df[!duplicated(meta_df$mouse_id),]
  meta_df <- meta_df %>%
    dplyr::select(passage, mouse_id)
  rownames(meta_df) <- meta_df$mouse_id
  
  pair_df <- read.csv(paste0(base_dir,'rnaseq_v6/incis_intrans/comparisons_drug_res.csv'), stringsAsFactors=F, check.names=F)
  # if(is.null(de)){
  #   de <- read.csv(paste0(base_dir,'rnaseq_v6/incis_intrans/',de_files[[datatag]]), stringsAsFactors=F, check.names=F)
  # }
  # if(is.null(meta)){
  #   meta <- read.csv(paste0(base_dir,'rnaseq_v6/incis_intrans/',meta_files[[datatag]]), stringsAsFactors=F, check.names=F)
  # }
  meta_tmp <- meta_df %>%
    dplyr::rename(passage1=passage, sample1=mouse_id)
  pair_df <- pair_df %>% left_join(meta_tmp, by=c('sample1'))
  meta_tmp <- meta_df %>%
    dplyr::rename(passage2=passage, sample2=mouse_id)
  pair_df[pair_df$clone1=='R' & pair_df$datatag=='SA609','clone1'] <- 'A'
  pair_df <- pair_df %>% left_join(meta_tmp, by=c('sample2'))
  pair_df$DE_desc <- paste0(pair_df$datatag,':Res(',pair_df$passage1,':',pair_df$clone1,') vs Sen(',pair_df$passage2,':',pair_df$clone2,')')
  
  pair_df$DE_desc1 <- paste0('Res(',pair_df$passage1,':',pair_df$clone1,') vs Sen(',pair_df$passage2,':',pair_df$clone2,')')
  pair_df$DE_desc2 <- paste0(pair_df$datatag,': Rx ',pair_df$passage1,':clone',pair_df$clone1,' vs UnRx ',pair_df$passage2,':clone',pair_df$clone2)
  
  write.csv(pair_df, paste0(base_dir,'rnaseq_v6/incis_intrans/comparisons_drug_res.csv'), quote=F, row.names = F)
  
  # SA609:Res(X7:A) vs Res(X7:A)
  # SA535:CX5461:Rx X7:cloneA
  # SA609:Rx X7:cloneA
  dim(de)
  # View(pair_df)
  dim(meta)
  de <- de[!duplicated(de$desc),]
  de <- de %>%
    dplyr::select(datatag, desc)
  de$de_desc <- gsub(paste0(datatag,'_'),'',de$desc)
  # View( head(de))
}

load_data <- function(pair_groups_fn, datatag, cnv_fn){
  pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  pair_groups <- pair_groups[pair_groups$datatag==datatag,] 
  
  # View(pair_groups)
  # pair_groups <- pair_groups %>%
  #   dplyr::filter(series!="")
  # pair_groups <- pair_groups %>%
  #   dplyr::filter(series=="")
  # pair_groups <- pair_groups %>%
  #     dplyr::filter(order<4 & order!=0)
  # data.table::fwrite(pair_groups, paste0(save_dir, 'tested_SA535_early.csv'))
  
  
  cnv_mat <- data.table::fread(cnv_fn) %>% as.data.frame()
  # head(cnv_mat)
  rownames(cnv_mat) <- cnv_mat$ensembl_gene_id
  cnv_mat$ensembl_gene_id <- NULL
  # cnv_mat <- cnv_mat %>% 
  #   tibble::column_to_rownames(var="ensembl_gene_id")
  # dim(cnv_mat)
  # colnames(cnv_mat)
  # pair_groups$desc <- paste0(pair_groups$series,'_',pair_groups$labels_detailed)
  pair_groups$desc <- pair_groups$file_header
  # pair_groups <- pair_groups %>% 
  #   dplyr::filter(clone1!='' & clone2 !='')
  pair_groups <- pair_groups %>%
    dplyr::filter(clone1!=clone2)
  
  return(list(pair_groups=pair_groups,cnv_mat=cnv_mat))
}

get_gene_type_edgeR_v2 <- function(pair_groups, cnv_mat_df, 
                                  datatag, input_dir, save_dir,
                                  FDR_cutoff=0.01, minLogFC=0.25, pValueThrs=0.05){  #FDR_cutoff=0.01, FDR_cutoff=0.05
  
  
  if (!file.exists(save_dir)){
    dir.create(save_dir, recursive=T)
  }
  
  rownames(pair_groups) <- pair_groups$desc
  frac_dlp <- list()
  genes_summary_stat <- list()
  genes_summary <- list()
  c <- 0
  for(de in pair_groups$desc){
    print('-------------------------------')
    print(de)
    cnv_mat <- cnv_mat_df
    output_dir <- paste0(save_dir,pair_groups[de,'file_header'],"/")
    if (!file.exists(output_dir)){
      dir.create(output_dir, recursive=T)
    }
    
    signif_genes <- data.table::fread(paste0(input_dir,pair_groups[de,'result_fn'])) %>% as.data.frame()
    # print(dim(signif_genes))
    # colnames(signif_genes)
    signif_genes <- signif_genes %>%
      dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
    # print('After genes filtering: ')
    # print(dim(signif_genes))
    # colnames(signif_genes)[which(names(signif_genes)=='gene_id')] <- 'ensembl_gene_id'
    # print(colnames(signif_genes))
    # signif_genes <- signif_genes[abs(signif_genes$logFC)>0.25,]
    # signif_genes <- signif_genes %>% left_join(cancer_ref_genes_df, by = c("ensembl_gene_id"="ensemble_id"))
    s1 <- as.character(pair_groups[de,'clone1'])
    s2 <- as.character(pair_groups[de,'clone2'])
    print(s1)
    print(s2)
    # print(colnames(cnv_mat))
    if(sum(colnames(cnv_mat) %in% c(s1,s2))!=2){
      print(de)
      stop("DEBUG: do not exist 2 clones in cnv median CN profiles")
    }
    
    cnv_mat <- cnv_mat[rownames(cnv_mat) %in% signif_genes$ensembl_gene_id,]
    # dim(cnv_mat)
    # dim(signif_genes)
    # var_genes <- apply(cnv_mat, 1, var)
    cnv_mat <- cnv_mat[,c(s1, s2)]
    # cnv_mat <- cnv_mat[rowSums(is.na(cnv_mat))==0 & var_genes > 0,]
    cnv_mat <- cnv_mat[rowSums(is.na(cnv_mat))==0,]
    
    # print(paste0('Nb bins with CN change: ',dim(cnv_mat)[1]))
    # t <- cnv_mat[,s1]-cnv_mat[,s2]
    cnv_mat$subtr <- cnv_mat[,s1]-cnv_mat[,s2]
    cnv_mat <- cnv_mat %>% 
      mutate(classified_gene_dlp = case_when(subtr>0 ~ 'Increase',
                                             subtr<0 ~ 'Decrease',
                                             subtr==0 ~ 'No_variance'))
    cnv_mat$subtr <- NULL
    # print(summary(as.factor(cnv_mat$classified_gene_dlp)))
    
    cnv_mat$ensembl_gene_id <- rownames(cnv_mat)
    # write.csv(cnv_mat_tmp, paste0(output_dir, "cnv_mat_total.csv"),quote=F, row.names=F)
    signif_genes <- signif_genes %>% left_join(cnv_mat, by = c("ensembl_gene_id"))
    rownames(signif_genes) <- signif_genes$ensembl_gene_id
    # summary(as.factor(signif_genes$classified_gene_dlp))
    
    signif_genes <- signif_genes %>% 
      mutate(classified_gene = case_when(classified_gene_dlp!='No_variance' & !is.na(classified_gene_dlp) ~ "in_cis",
                                         classified_gene_dlp=='No_variance' & !is.na(classified_gene_dlp) ~ "in_trans",
                                         is.na(classified_gene_dlp) ~ "unmapped"))
    
    # print(summary(as.factor(signif_genes$classified_gene)))
    
    
    signif_genes <- signif_genes %>% 
      mutate(classified_gene_rna = case_when(logFC>0 ~ 'UpRegulated',
                                             logFC<0 ~ 'DownRegulated',
                                             is.na(logFC) ~ 'Invalid'))
    
    signif_genes <- signif_genes %>% 
      mutate(Gene_Type = case_when(signif_genes$classified_gene=='in_cis' ~ paste0('In_cis_',classified_gene_dlp, '_',classified_gene_rna),
                                   signif_genes$classified_gene=='in_trans' ~ paste0('In_trans_',classified_gene_rna),
                                   signif_genes$classified_gene== 'unmapped' ~ "unmapped"))
    
    # signif_genes$Gene_Type <- as.factor(signif_genes$Gene_Type)
    # plot_DE(de, signif_genes, pair_groups, output_dir,
    #         minLogFC=0.25, pValueThrs=0.05, nbtopup=30, nbtopdown=30)
    print('SUMMARY cis positive percentage:')
    # print(summary(as.factor(signif_genes$Gene_Type)))
    st <- summary(as.factor(signif_genes$Gene_Type))
    cis_pos <- as.numeric(st['In_cis_Decrease_DownRegulated']) + as.numeric(st['In_cis_Increase_UpRegulated'])
    cis_neg <- as.numeric(st['In_cis_Decrease_UpRegulated'])+as.numeric(st['In_cis_Increase_DownRegulated'])
    pct <- round(100*cis_pos/(cis_pos + cis_neg),2)
    print(pct)
    data.table::fwrite(signif_genes, paste0(output_dir, 'signif_genes.csv'))
  }
  
  
  # frac_dlp_df <- as.data.frame(do.call(rbind, frac_dlp))
  # colnames(frac_dlp_df) <- c('result_fn','de_analysis','total_var_dlp',
  #                            'incis','incis_incr','incis_decr','pct_incis')# 
  # frac_dlp_df$desc <- rownames(frac_dlp_df)
  # frac_dlp_df$pct_incis_incr <- round(as.numeric(frac_dlp_df$incis_incr)/as.numeric(frac_dlp_df$total_var_dlp) * 100, 2)
  # frac_dlp_df$pct_incis_decr <- round(as.numeric(frac_dlp_df$incis_decr)/as.numeric(frac_dlp_df$total_var_dlp) * 100, 2)
  # frac_dlp_df$pct_others <- 100 - frac_dlp_df$pct_incis_incr - frac_dlp_df$pct_incis_decr
  # View(frac_dlp_df)
  # frac_dlp_df$pct_incis_inc <- round(as.numeric(frac_dlp_df$incis_incr) / as.numeric(frac_dlp_df$total_var_dlp) * 100,2)
  # frac_dlp_df$pct_incis_dec <- round(as.numeric(frac_dlp_df$incis_decr) / as.numeric(frac_dlp_df$total_var_dlp) * 100,2)
  # frac_dlp_df$unaffected <- 100 - frac_dlp_df$pct_incis_inc - frac_dlp_df$pct_incis_dec
  # write.csv(frac_dlp_df, paste0(input_dir,datatag,'_',subtag,'_fraction_dlp.csv'), quote = F, row.names = F)
  
  # stat <- do.call(rbind, genes_summary_stat)
  # View(head(stat))
  # stat$gene_type <- rownames(stat)
  # data.table::fwrite(stat, paste0(save_dir, datatag,'_cistrans_summary.csv'), quote = F)
  # saveRDS(genes_summary,paste0(save_dir, datatag,'_summary_statistic.rds'))
  # return(stat)
  return(signif_genes)
}



get_top_genes <- function(de_genes, minLogFC=0.25, nbtopup=25, nbtopdown=25){
  markers_ls_upreg <- de_genes[de_genes$logFC>minLogFC,]
  markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$logFC,decreasing = T),] 
  markers_ls_upreg <- markers_ls_upreg[1:nbtopup,]
  markers_ls_downreg <- de_genes[de_genes$logFC<(-minLogFC),]
  markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$logFC,decreasing = F),] 
  markers_ls_downreg <- markers_ls_downreg[1:nbtopdown,]
  topGenes <- c(as.character(markers_ls_upreg$gene_symbol),as.character(markers_ls_downreg$gene_symbol))
  return(topGenes)
}

# colData contain clone, cell_id columns, and row names are cell ids
# data: counts dataframe with colnames are cell ids
extract_genes_regression_multi_clones <- function(data, clones_df, pair_clones, colData,
                                                  outputPrefix, condition='clone',
                                                  computeMeanVal=FALSE, 
                                                  subsample_data=TRUE, maxExtCells=1000){
  library(DESeq2)
  rownames(clones_df) <- clones_df$obs_clones
  # sensitive_treatment_st <- c('UT','U','UU','UUU','UUUU','UUUUU')
  # resistant_treatment_st <- c('UT','UTT','UTTT','UTTTT','UTTTTT')
  sensitive_treatment_st <- c("UUUU")
  resistant_treatment_st <- c('UTTT')
  
  signif_threshold <- 0.05
  colData = colData[colData$clone %in% clones_df$obs_clones, ]
  print(summary(as.factor(colData$treatment_st)))
  
  ext_cells <- c()
  if(subsample_data){
    clones <- unique(colData$clone)
    for(c in clones){
      if(clones_df[c,'treatment_effect']=='sensitive'){
        obs_ts <- sensitive_treatment_st
      }else{
        obs_ts <-resistant_treatment_st
      } 
      print(paste0("Clone is: ",c, ' obs treatment status are: '))
      print(obs_ts)
      # tmp <- colData[colData$clone==c & colData$treatment_st %in% obs_ts,]
      tmp <- colData[colData$clone==c,]
      if(!is.null(tmp) & nrow(tmp) < maxExtCells){
        cells_to_sample <- tmp$cell_id
      } else if(!is.null(tmp) & nrow(tmp) >= maxExtCells){
        cells_to_sample <- sample(tmp$cell_id, maxExtCells, replace = FALSE)
      } else{
        cells_to_sample <- NULL
        print(paste0('DEBUG: check input cells name for clone', c))
        # obs_clones <- obs_clones[obs_clones!=c]
        clones_df <- clones_df[clones_df$obs_clones != c,]
      }
      
      if(length(cells_to_sample)>0){
        print(paste0('Extract ',length(cells_to_sample),' from clone', c))
        ext_cells <- c(ext_cells, cells_to_sample)
      }
      
    }
    ext_cells <- intersect(colnames(data), ext_cells)
  }else{
    ext_cells <- intersect(colnames(data),rownames(colData))
  }
  print(paste0('Number of obs cells: ',length(ext_cells)))
  
  
  colData <- colData[ext_cells,]
  colData$clone <- as.factor(colData$clone)
  print(dim(colData))
  data <- data[,ext_cells]
  print(dim(data))
  
  # print('Creating all pair of comparison...')
  # pair_clones <- combinations(length(obs_clones),2, obs_clones)
  # print(paste0('Nb of pair comparisons: ',nrow(pair_clones)))
  # print(pair_clones)
  
  print('Creating DESeq object...')
  dds_fn <- paste0(outputPrefix, "dds.rds")
  if(file.exists(dds_fn)){
    dds <- readRDS(dds_fn)
  }else{
    dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ clone)
    # #----- estimate dispersions -----#
    dds <- estimateSizeFactors(dds, type = "poscounts")
    # dds <- estimateDispersions(dds, fitType='local')
    saveRDS(dds, dds_fn)
  }
  
  print('Multi-conditions testing...')
  cores_use <- 3
  nbcores <- detectCores()
  if(cores_use > nbcores){
    cores_use <- nbcores
  }
  p <- MulticoreParam(cores_use)
  dds_deseq_fn <- paste0(outputPrefix, "dds_deseq.rds")
  if(file.exists(dds_deseq_fn)){
    dds_rg <- readRDS(dds_deseq_fn)
  }else{
    dds_rg <- DESeq(dds, fitType = "local", parallel = TRUE, BPPARAM = p)  # TO DO: parallel computation
    saveRDS(dds_rg, dds_deseq_fn)
  }  
  
  DE_results <- list()
  for(i in rep(1:nrow(pair_clones),1)){
    lb <- paste0(pair_clones[i,1],pair_clones[i,2])
    print(lb)
    res <- results(dds_rg, contrast = c("clone", pair_clones[i,1], pair_clones[i,2]),
                   parallel = TRUE, BPPARAM = p)
    rownames(res) <- rownames(data)
    
    res_trans <- data.frame(transform(res, sample1 = pair_clones[i,1], sample2 = pair_clones[i,2]),
                            row.names = paste0(rownames(res), "_", lb))
    DE_results[[i]] <- res_trans
  }
  length(DE_results)
  dim(res_trans)
  all <- do.call(rbind, DE_results)
  print(dim(all))
  
  print('FDR test, BH adjusted test...')
  q <- p.adjust(all$pval, method = "BH")
  result <- transform(all, qval = q)
  print(summary(result$qval))
  signif_genes_df <- subset(result, qval < signif_threshold)
  signif_genes_df$gene_id <- get_gene_id(rownames(signif_genes_df), cores_use) 
  print(paste0('Nb significant genes: ',length(unique(signif_genes_df$gene_id))))
  # saveRDS(result, paste0(outputPrefix, "results.rds"))
  
  print('Mean exp for each clone...')
  normalizedCount <- as.data.frame(counts(dds_rg, normalized = TRUE))
  print(dim(normalizedCount))
  
  # saveRDS(normalizedCount, paste0(outputPrefix, "normalized_count.rds"))
  total_data_results <- list(normalizedCount=normalizedCount,
                             colData=colData,
                             data=data,
                             obs_clones=clones_df,
                             result=result,
                             signif_genes_df=signif_genes_df,
                             outputPrefix=outputPrefix)
  
  saveRDS(total_data_results, paste0(outputPrefix, "total_data_results.rds"))
  
  if(computeMeanVal){
    gene_mean <- compute_normalized_mean_gene(normalizedCount, clones_df$obs_clones, 
                                              colData, outputPrefix, cores_use)
    saveRDS(gene_mean, paste0(outputPrefix, "gene_mean.rds"))
    write.table(gene_mean, paste(outputPrefix, "normalized_mean", sep = ""),
                col.names = F, row.names = T, append = F, quote = F, sep = "\t")
  }
  
  # data.table::fwrite(normalizedCount, file=paste(outputPrefix, "normalized_count", sep = ""),
  #                    col.names = T, row.names = T, append = F, quote=F, sep = "\t")
  # 
  # write.table(normalizedCount, paste(outputPrefix, "normalized_count", sep = ""), 
  #             col.names = F, row.names = T, append = F, quote = F, sep = "\t")
  
  write.table(result, paste(outputPrefix, "results", sep = ""),
              col.names = T, row.names = T, quote = F, sep = "\t")
  if(!is.null(signif_genes_df) && nrow(signif_genes_df)>0){
    print(paste0('Number of significant genes: ',nrow(signif_genes_df)))
    write.table(signif_genes_df, paste(outputPrefix, "significant", sep = ""), col.names = T, row.names = T, quote = F, sep = "\t")
  }
  
  
}


load_reference_genes_set <- function(ref_dif, ref_fn=NULL){
  pathway_name=c('cisplatin_resistance','core_fitness') #'cosmic' or 'cisplatin_resistance', or 'metastasis'
  # ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  if(is.null(ref_dif)){
    ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'  
  }
  if(is.null(ref_fn)){
    ref_fn <- paste0(ref_dif, 'reference_cancer_genes_set_CF_CR.rds')  
  }
  if(file.exists(ref_fn)){
    ref_genes_ls <- readRDS(ref_fn)
    # print(names(ref_genes_ls))
    return(ref_genes_ls)
  }else{
    ref_genes_ls <- list()
    if('cosmic' %in% pathway_name){
      ref_genes <- read.table(paste0(ref_dif, 'oncogene_cosmic.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
      print(length(ref_genes$Gene_Symbol))  
      ref_genes_ls[['COSMIC']] <- toupper(ref_genes$Gene_Symbol)
    }
    if('cisplatin_resistance' %in% pathway_name){
      ref_genes <- data.table::fread(paste0(ref_dif, 'cisplatin_resistance_corrected.csv')) %>% as.data.frame()
      ref_genes <- ref_genes[!duplicated(ref_genes$gene_symb),, drop=F]
      ref_genes_ls[['CISPLATIN_RESISTANCE']] <- toupper(ref_genes$gene_symb)
    }
    if('broadsanger' %in% pathway_name){
      ref_genes <- data.table::fread(paste0(ref_dif, 'integrated_broad_sanger_essential_983_genes.csv')) %>% as.data.frame()
      print(length(ref_genes$gene_symbol))  
      ref_genes_ls[['BS_ESSENTIAL_CANCER']] <- toupper(ref_genes$gene_symbol)
    }
    if('core_fitness' %in% pathway_name){
      ref_genes <- data.table::fread(paste0(ref_dif, 'CoreFitness_corrected.csv')) %>% as.data.frame()
      ref_genes <- ref_genes[!duplicated(ref_genes$gene_symb),, drop=F]
      ref_genes_ls[['CORE_FITNESS']] <- toupper(ref_genes$gene_symb)
    }
    print(names(ref_genes_ls))
    saveRDS(ref_genes_ls, ref_fn)
    return(ref_genes_ls)
  }
  
}

