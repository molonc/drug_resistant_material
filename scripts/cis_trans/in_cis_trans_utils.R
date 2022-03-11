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
library(extrafont)
font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
fonts()
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


viz_pancancer_genes <- function(cancer_pct_df, output_dir, tag=""){
  p <- ggplot(data=cancer_pct_df, aes(x=PDX, y=Fitness_genes_pct, fill=series)) +
    geom_bar(stat="identity", width=0.7)+ #, fill="white"
    geom_text(aes(label=Fitness_genes_pct), vjust=-0.3, size=3.5, color="black")+
    theme_minimal()
  p <- p + scale_fill_manual(values=c("#30F50A", "#8DD87F", "#1A6A0B", "#560C56"))
  
  p <- p + labs(x='PDX drug treatment', y="% Fitness_genes", title='Percentage genes in ADAM PanCancer Core-Fitness genes')
  p <- p + theme(legend.title = element_text(size=10), 
                 legend.text = element_text(size=9),
                 plot.title = element_text(color="black", size=13, hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 # panel.border = element_blank(),
                 axis.text.x = element_text(size=9, angle = 90, color="black"),
                 #axis.title = element_blank(),
                 axis.ticks.x = element_blank()
  )
  # p
  png(paste0(output_dir,tag,"_fitness_genes_stat.png"), height = 2*400, width=2*620, res = 2*72)
  print(p)
  dev.off()
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

viz_heatmap <- function(input_dir, output_dir){
  library(ComplexHeatmap)
  SA609_summary <- read.csv(paste0(input_dir,'SA609/SA609_result/SA609_summary.csv'), 
                            stringsAsFactors=F, check.names = F)
  # dim(SA609_summary)
  # View(SA609_summary)
  SA1035_summary <- read.csv(paste0(input_dir,'SA1035/SA1035_result/SA1035_summary.csv'), 
                             stringsAsFactors=F, check.names = F)
  
  SA535_cis_summary <- read.csv(paste0(input_dir,'SA535/SA535_cis_result/SA535_cis_summary.csv'), 
                                stringsAsFactors=F, check.names = F)
  SA535_cx_summary <- read.csv(paste0(input_dir,'SA535/SA535_cx_result/SA535_cx_summary.csv'), 
                               stringsAsFactors=F, check.names = F)
  
  
  stat_df <- list(SA609_summary, SA1035_summary, SA535_cis_summary, SA535_cx_summary)
  stat_df2 <- list()
  c <- 0
  for(df in stat_df){
    rownames(df) <- df$gene_type
    df <- df[,colnames(df) !='gene_type']
    exclude_cols <- grep('_nbgenes',colnames(df), value = T)
    cols_use <- colnames(df)[!colnames(df) %in% exclude_cols]
    df <- df[,cols_use]
    c <- c + 1
    stat_df2[[c]] <- t(df)
  }
  # View(SA535_cx_summary)
  stat <- do.call(rbind, stat_df2)
  # View(stat)
  # rownames(stat)
  coln <- colnames(stat)
  coln <- gsub('In_cis','InCis',coln)
  coln <- gsub('In_trans','InTrans',coln)
  colnames(stat) <- coln
  gt <- lapply(strsplit(colnames(stat), "_"), function(x) {
    return(x[1])
  })
  
  annotation_col = data.frame(
    Gene_Type = factor(as.character(gt))#, 
    #GeneType = 1:6
  )
  # dim(annotation_col)
  rownames(annotation_col) = colnames(stat)
  
  series <- lapply(strsplit(rownames(stat), "-"), function(x) {
    return(x[1])
  })
  annotation_row = data.frame(
    PDX = factor(as.character(series))
  )
  # dim(annotation_row)
  rownames(annotation_row) = rownames(stat)
  
  top_anno = ComplexHeatmap::HeatmapAnnotation(Gene_Type = factor(as.character(gt)),
                                               col = list(Gene_Type = c(InCis = "#FF3232", InTrans = "#40A0E0")))
  left_anno = ComplexHeatmap::rowAnnotation(PDX = factor(as.character(series)),
                                            col = list(PDX=c(SA609="#1A6A0B",SA1035 = "#8DD87F", SA535_cis = "#30F50A", SA535_cx = "#560C56"))  #SA609="#1A6A0B"
  )
  
  
  cell_func = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", stat[i, j]), x, y, gp = gpar(fontsize = 13))
  }
  
  # ha = HeatmapAnnotation(summary = anno_summary(height = unit(4, "cm")))
  p <- ComplexHeatmap::Heatmap(stat, na_col = "white",
                               show_column_names=T,
                               show_row_names = T,
                               cluster_rows=F,cluster_columns=F,
                               name = "% Genes", 
                               row_order = rownames(stat), ## TO DO: check rows order
                               row_split= annotation_row,
                               row_title_rot = 0,
                               row_gap = unit(2, "mm"),
                               column_split = annotation_col, 
                               column_title = "Genes Copy Number Driven",
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 10),
                               row_names_gp = grid::gpar(fontsize = 10),
                               show_heatmap_legend = T,
                               top_annotation=top_anno,
                               left_annotation = left_anno,
                               cell_fun = cell_func,
                               row_dend_reorder=F)
  
  
  png(paste0(output_dir,'summary_incis_intrans_genes.png'), height = 2*700, width=2*1000,res = 2*72)
  print(p)
  dev.off()
  
  # View(stat)
  write.csv(stat, paste0(output_dir,'summary_incis_intrans_genes.csv'), quote = F, row.names = T)
  
}

viz_prevalence <- function(input_dir, output_dir){
  stat <- read.csv(paste0(output_dir,'summary_incis_intrans_genes.csv'), check.names = F, stringsAsFactors=F, row.names = 1)
  # dim(stat)
  # View(stat)
  
  # stat$InCis_Increase <- stat$InCis_Increase_DownRegulated + stat$InCis_Increase_UpRegulated
  # stat$InCis_Decrease <- stat$InCis_Decrease_DownRegulated + stat$InCis_Decrease_UpRegulated
  # stat$InTrans <- stat$InTrans_DownRegulated + stat$InTrans_UpRegulated
  # in-cis up (proportion of CN change up)
  # in-cis down (proportion of CN change down)
  # in-trans (proportion of in-trans)
  # stat <- stat[,c('InCis_Increase','InCis_Decrease','InTrans')]
  
  stat1 <- stat
  # stat1$desc <- as.character(rownames(stat1))
  
  stat1 <- stat1 %>%
    pivot_longer(!desc, names_to = "gene_type", values_to = "percentage")
  
  series <- mclapply(strsplit(stat1$desc, "-"), function(x) {
    return(x[1])
  }, mc.cores = 2)
  
  desc <- mclapply(strsplit(stat1$desc, "-"), function(x) {
    return(x[2])
  }, mc.cores = 2)
  
  stat1$desc <- as.character(desc)
  stat1$series <- as.character(series)
  stat1$series <- factor(stat1$series, levels = unique(stat1$series))
  p <- ggplot(stat1, aes(fill=gene_type, y=percentage, x=desc)) + 
    geom_bar(position="fill", stat="identity")+
    facet_grid(. ~ series, scales="free", space='free')
  p <- p + labs(x='Resistant vs Sensitive cells comparisons', y="Percentage gene type (%)", title='Genes - Copy Number Driven Percentage')
  p <- p + theme(legend.title = element_text(size=10), 
                 legend.text = element_text(size=7),
                 plot.title = element_text(color="black", size=12, hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.x = element_text(size=10, angle = 90, color="black"),
                 axis.text.y = element_text(size=10, color="black")
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  
  png(paste(output_dir,"incis_intrans_genes_prevalence_summary.png",sep="/"), height = 2*380, width=2*700, res = 2*72)
  print(p)
  dev.off()
}


plot_DE_genes_edgeR <- function(df, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                plttitle="A versus B", save_dir="",legendVisible=F,
                                iscaption=TRUE, legend_verbose='none', save_plot=TRUE){  
  #legend_verbose is none or 'right', 'left','bottom'
  # df <- de_genes
  # library(EnhancedVolcano)
  
  # colnames(df)[which(names(df) == "avg_logFC")] <- "log2FoldChange"
  # colnames(df)[which(names(df) == "p_val_adj")] <- "padj"
  # capstr <- paste0("FDR cutoff, ",FDRcutoff,"; logFC cutoff, ",logFCcutoff, "; nb genes signif, ",nbgenessig)
  # summary(as.factor(df$gene_type))
  Gene_Type=c('In_cis_Decrease_DownRegulated','In_cis_Decrease_UpRegulated',
              'In_cis_Increase_DownRegulated','In_cis_Increase_UpRegulated',
              'In_trans_DownRegulated','In_trans_UpRegulated'
  )
  gt <- data.frame(Gene_Type=Gene_Type,
                   gene_type=c('In-cis_DD','In-cis_DU',
                               'In-cis_ID','In-cis_IU',
                               'In-trans_D','In-trans_U'),
                   col_gt=c('#E61367','#D35B53',
                            '#EF816E','#E82503',
                            '#7BEADF','#093CB4'))
  
  df <- df %>% left_join(gt, by='Gene_Type')
  # unique(df$col_gt)
  keyvals_colour <- as.character(df$col_gt)
  names(keyvals_colour) <- df$gene_type
  # keyvals_colour <- factor(keyvals_colour, levels = unique(keyvals_colour))
  # unique(keyvals_colour)
  # names(keyvals_colour[1:3])
  df <- df[abs(df$logFC)>logFCcutoff & df$FDR<FDRcutoff  & df$PValue<pValuecutoff,]
  df$logFC <- sapply(df$logFC, function(x) replace(x, x > 3, 3))
  df$logFC <- sapply(df$logFC, function(x) replace(x, x < (-3), -3))
  
  st <- summary(df$gene_type)
  keyvals_shape <- ifelse(df$is_fitness_gene==T, 3, 16)
  names(keyvals_shape) <- ifelse(df$is_fitness_gene==T,'Fitness genes','Others')
  
  # df$gene_type <- factor(df$gene_type, levels = unique(df$gene_type))
  if(capstr=='' & iscaption){
    # capstr <- paste0(capstr,'With abs(logFC)>0.25, FDR<0.01, pValue<0.05 \n')
    capstr <- paste0(capstr,names(st[1]),':',as.numeric(st[1]), ', ')
    capstr <- paste0(capstr,names(st[2]),':',as.numeric(st[2]), ', ')
    capstr <- paste0(capstr,names(st[3]),':',as.numeric(st[3]), ', ')
    capstr <- paste0(capstr,names(st[4]),':',as.numeric(st[4]), ' \n')
    capstr <- paste0(capstr,names(st[5]),':',as.numeric(st[5]), ', ')
    capstr <- paste0(capstr,names(st[6]),':',as.numeric(st[6]), ' ')
    # for(i in rep(1:length(st),1)){
    #   # print(st[i])
    #   capstr <- paste0(capstr,names(st[i]),':',as.numeric(st[i]), ' ')
    # }
    
  }
  df$mlog10FDR <- -log10(df$FDR)
  
  df$mlog10FDR <- sapply(df$mlog10FDR, function(x) replace(x, is.infinite(x), 300))
  # legend_verbose <- 'top'
  p <- EnhancedVolcano::EnhancedVolcano(df,
                                        lab = df$gene_symbol,   #df$symbol or rownames(df)
                                        x = 'logFC',
                                        y = 'FDR',
                                        # selectLab = as.character(topGenes),
                                        selectLab = topGenes,
                                        xlim=c(-3,3),
                                        ylim = c(0, max(df$mlog10FDR, na.rm = TRUE)+10),
                                        xlab = bquote(~Log[2] ~ ' fold change'),
                                        ylab = bquote(~-Log[10]~italic(FDR)),
                                        pCutoff = FDRcutoff,
                                        FCcutoff = logFCcutoff,
                                        # pointSize = 3.2,
                                        pointSize = c(ifelse(df$is_fitness_gene==T, 4, 3)),
                                        labSize = 3,
                                        colAlpha = 0.8,
                                        gridlines.major = FALSE,
                                        gridlines.minor = FALSE,
                                        colCustom = keyvals_colour,
                                        shapeCustom = keyvals_shape,
                                        cutoffLineType = 'twodash',
                                        cutoffLineCol = 'grey0',
                                        cutoffLineWidth = 0.5,
                                        # hline = c(0.01),
                                        # hlineCol = c('grey0'),
                                        # hlineType = 'longdash',
                                        legend=NULL,
                                        # legend=c('NS','Log2 FC','Adjusted p-value',
                                        #          'Adjusted p-value & Log2 FC'),
                                        legendPosition = legend_verbose,
                                        legendLabSize = 8,
                                        legendIconSize = 3.0,
                                        legendVisible = legendVisible,
                                        # col=c('black', 'black', 'black', 'red3'),
                                        col=unique(keyvals_colour),
                                        title = NULL,
                                        subtitle = plttitle,
                                        caption = capstr,
                                        drawConnectors = TRUE,
                                        widthConnectors = 0.15,
                                        colConnectors = 'grey30')
  # p
  # p <- p + theme(legend.position="none")
  
  tag <- ifelse(legend_verbose=='none','','_with_legend')
  png(paste0(save_dir,"DE_",plttitle,tag,".png"), height = 2*550, width=2*650,res = 2*72)
  print(p)
  dev.off()
  if(save_plot){
    plttitle <- gsub(':','_',plttitle)
    plttitle <- gsub(' ','_',plttitle)
    saveRDS(p, file=paste0(save_dir,"DE_",plttitle,".rds"))
  }
  
  return(p)
}


plot_DE <- function(de, de_genes, pair_groups, output_dir,
                    minLogFC=0.25, pValueThrs=0.05,
                    nbtopup=30, nbtopdown=30){
  FDR_cutoff <- 0.01
  nbtop_extract <- 20
  de_genes <- de_genes[abs(de_genes$logFC)>minLogFC & de_genes$FDR<FDR_cutoff  & de_genes$PValue<pValueThrs,]
  # Plot DE genes
  # plttitle <- paste0(pair_groups[de,'datatag'],":  ",pair_groups[de,'clone1']," versus ",pair_groups[de,'clone2'])
  plttitle <- de
  
  markers_ls_upreg <- de_genes[de_genes$logFC>minLogFC,]
  markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$logFC,decreasing = T),] 
  dim(markers_ls_upreg)
  # pair_groups[de,'resistant_genes'] <- nrow(markers_ls_upreg)
  
  # Extract top incis, intrans up-regulated
  markers_ls_upreg_incis <- markers_ls_upreg[!is.na(markers_ls_upreg$classified_gene_dlp),]
  dim(markers_ls_upreg_incis)
  markers_ls_upreg_incis <- markers_ls_upreg_incis[order(markers_ls_upreg_incis$logFC,decreasing = T),] 
  if(nrow(markers_ls_upreg_incis) > nbtop_extract){
    markers_ls_upreg_incis <- markers_ls_upreg_incis[1:nbtop_extract,]
  }
  write.csv(markers_ls_upreg_incis, file=paste0(output_dir,'topgenes_upregulated_incis.csv'), quote = F, row.names = F)
  
  markers_ls_upreg_intrans <- markers_ls_upreg[is.na(markers_ls_upreg$classified_gene_dlp),]
  dim(markers_ls_upreg_intrans)
  markers_ls_upreg_intrans <- markers_ls_upreg_intrans[order(markers_ls_upreg_intrans$logFC,decreasing = T),] 
  if(nrow(markers_ls_upreg_intrans) > nbtop_extract){
    markers_ls_upreg_intrans <- markers_ls_upreg_intrans[1:nbtop_extract,]
  }
  write.csv(markers_ls_upreg_intrans, file=paste0(output_dir,'topgenes_upregulated_intrans.csv'), quote = F, row.names = F)
  
  
  # dim(markers_ls_upreg)
  if(nrow(markers_ls_upreg) < nbtopup){
    nbtopup <- nrow(markers_ls_upreg)
  }
  markers_ls_upreg <- markers_ls_upreg[1:nbtopup,]
  
  # Get top down-regulated genes
  markers_ls_downreg <- de_genes[de_genes$logFC<(-minLogFC),]
  markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$logFC,decreasing = F),] 
  
  # Extract top incis, intrans down-regulated
  markers_ls_downreg_incis <- markers_ls_downreg[!is.na(markers_ls_downreg$classified_gene_dlp),]
  dim(markers_ls_downreg_incis)
  markers_ls_downreg_incis <- markers_ls_downreg_incis[order(markers_ls_downreg_incis$logFC,decreasing = F),] 
  if(nrow(markers_ls_downreg_incis) > nbtop_extract){
    markers_ls_downreg_incis <- markers_ls_downreg_incis[1:nbtop_extract,]
  }
  write.csv(markers_ls_downreg_incis, file=paste0(output_dir,'topgenes_downregulated_incis.csv'), quote = F, row.names = F)
  
  markers_ls_downreg_intrans <- markers_ls_downreg[is.na(markers_ls_downreg$classified_gene_dlp),]
  dim(markers_ls_downreg_intrans)
  markers_ls_downreg_intrans <- markers_ls_downreg_intrans[order(markers_ls_downreg_intrans$logFC,decreasing = F),] 
  if(nrow(markers_ls_downreg_intrans) > nbtop_extract){
    markers_ls_downreg_intrans <- markers_ls_downreg_intrans[1:nbtop_extract,]
  }
  write.csv(markers_ls_downreg_intrans, file=paste0(output_dir,'topgenes_downregulated_intrans.csv'), quote = F, row.names = F)
  
  
  if(nrow(markers_ls_downreg) < nbtopdown){
    nbtopdown <- nrow(markers_ls_downreg)
  }
  markers_ls_downreg <- markers_ls_downreg[1:nbtopdown,]
  topGenes <- c(as.character(markers_ls_upreg$gene_symbol),as.character(markers_ls_downreg$gene_symbol))
  genes_df <- data.frame(desc=rep(de, length(topGenes)), topGenes=topGenes)
  write.csv(genes_df, file=paste0(output_dir,'topgenes.csv'), quote = F, row.names = F)
  
  # rownames(markers_ls_tmp) <- markers_ls_tmp$gene_symb
  
  plot_DE_genes_edgeR(de_genes, topGenes, capstr='', 
                      FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                      plttitle, output_dir, legendVisible=F,
                      iscaption=TRUE, legend_verbose='none', save_plot=TRUE)
  plot_DE_genes_edgeR(de_genes, topGenes, capstr='', 
                      FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                      plttitle, output_dir, legendVisible=T,
                      iscaption=TRUE, legend_verbose='top', save_plot=FALSE)
  
    
}

plot_dlp_prevalence <- function(output_dir=NULL){
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/SA535-v6/'
  base_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/')
  dlp609 <- read.csv(paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609__fraction_dlp_corrected.csv'), check.names=F, stringsAsFactors=F)
  dlp1035 <- read.csv(paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/SA1035__fraction_dlp.csv'), check.names=F, stringsAsFactors=F)
  dlpSA535_cis <- read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_SA535_cisplatin_fraction_dlp.csv'), check.names=F, stringsAsFactors=F)
  dlpSA535_cx <- read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_SA535_CX5461_fraction_dlp.csv'), check.names=F, stringsAsFactors=F)
  
  # View(dlp609)
  # dlp609$desc_2 <- gsub('SA609_','',dlp609$desc)
  # dlp1035$desc_2 <- gsub('SA1035_','',dlp1035$desc)
  # dlpSA535_cis$desc_2 <- gsub('SA535_','',dlpSA535_cis$desc)
  # dlpSA535_cx$desc_2 <- gsub('SA535_','',dlpSA535_cx$desc)
  stat_dlp <- do.call(rbind,list(dlp609, dlp1035, dlpSA535_cis, dlpSA535_cx))
  # View(stat_dlp)
  stat_dlp$pct_incis_incr <- NULL
  stat_dlp$pct_incis_decr <- NULL
  stat_dlp$pct_others <- NULL 
  stat_dlp$pct_incis_incr <- round(as.numeric(stat_dlp$incis_incr)/as.numeric(stat_dlp$total_var_dlp) * 100, 2)
  stat_dlp$pct_incis_decr <- round(as.numeric(stat_dlp$incis_decr)/as.numeric(stat_dlp$total_var_dlp) * 100, 2)
  stat_dlp$pct_others <- 100 - stat_dlp$pct_incis_incr - stat_dlp$pct_incis_decr
  
  
  write.csv(stat_dlp, paste0(input_dir,'summary_dlp.csv'), row.names = F, quote = F)
  colnames(stat_dlp)
  stat_dlp1 <- stat_dlp[,c("desc","pct_incis_incr","pct_incis_decr","pct_others")]
  colnames(stat_dlp1) <- c("desc","Incis_incr","Incis_decr","Others_CNchange")
 
  stat_dlp1 <- stat_dlp1 %>%
    tidyr::pivot_longer(!desc, names_to = "gene_type", values_to = "percentage")
  
  dim(stat_dlp1)
  PDX <- lapply(strsplit(stat_dlp1$desc, "_"), function(x) {
    return(x[1])
  })
  unique(stat_dlp1$gene_type)
  stat_dlp1$PDX <- as.character(PDX)
  stat_dlp1$PDX <- ifelse(grepl('X',stat_dlp1$desc),'SA535_CX5461',
                          ifelse(!grepl('X',stat_dlp1$desc),paste0(stat_dlp1$PDX,'_CISPLATIN'),stat_dlp1$PDX))
  
  cols_dlp <- c("#E82503","#E61367","#00B300")
  stat_dlp1$gene_type <- factor(stat_dlp1$gene_type, levels = unique(stat_dlp1$gene_type))
  stat_dlp1$PDX <- factor(stat_dlp1$PDX, levels = unique(stat_dlp1$PDX))
  # stat_dlp1$desc <- gsub('SA609_','',stat_dlp1$desc)
  stat_dlp1$desc <- str_replace_all(stat_dlp1$desc, "(SA609_)|(SA535_)|(SA1035_)", "")
  p <- ggplot(stat_dlp1, aes(fill=gene_type, y=percentage, x=desc)) + 
    geom_bar(position="fill", stat="identity", width = 0.6) +
    scale_y_continuous(labels = scales::percent_format()) +
    facet_grid(. ~ PDX, scales="free", space='free') + 
    scale_fill_manual(values = cols_dlp) + 
    labs(x='DE Analysis: Resistant vs Sensitive cells', y="(%) genes in total CN change ", title='Proportion in-cis genes in total copy number change') +
    theme(legend.title = element_text(size=8), 
          legend.text = element_text(size=7),
          plot.title = element_text(color="black", size=13, hjust = 0.5),
          # legend.position = "none", 
          # axis.line = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          axis.text.x = element_text(size=8, angle = 90, color="black"),
          axis.text.y = element_text(size=8, color="black")
          #axis.title = element_blank(),
          # axis.ticks = element_blank()
    )
  
  png(paste0(output_dir,"cn_change_proportion.png"), height = 2*420, width=2*720, res = 2*72)
  print(p)
  dev.off()
  write.csv(stat_dlp,paste0(output_dir,'cn_change_proportion.csv'), row.names = F, quote = F)
  
  saveRDS(p, paste0(output_dir,"cn_change_proportion.rds"))
  
}

plot_intrans_incis_prevalence_v2 <- function(){
  
  # base_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/')
  # datatag <- 'SA609'
  # stat609 <- data.table::fread(paste0(base_dir,datatag,'_rna/cis_trans/',datatag,'_summary.csv')) %>% as.data.frame()
  # dim(stat609)
  # save_dir <- paste0(base_dir,datatag,'_rna/cis_trans/')
  # # View(stat609)
  # datatag <- 'SA1035'
  # stat1035 <- data.table::fread(paste0(base_dir,datatag,'_rna/cis_trans/',datatag,'_summary.csv')) %>% as.data.frame()
  # dim(stat1035)
  # # View(stat1035)
  # datatag <- 'SA535'  
  # stat535 <- data.table::fread(paste0(base_dir,datatag,'_rna/cis_trans/',datatag,'_summary.csv')) %>% as.data.frame()
  # dim(stat535)
  # stat609$PDX <- 'SA609'
  # stat1035$PDX <- 'SA1035'
  # stat535$PDX <- 'SA535'
  # stat <- dplyr::bind_rows(list(stat609, stat1035, stat535))
  # stat$PDX <- c(rep('SA609',nrow(stat609)),rep('SA1035',nrow(stat1035)),rep('SA535',nrow(stat535)))
  gene_type=c('In_cis_Decrease_DownRegulated','In_cis_Decrease_UpRegulated',
              'In_cis_Increase_DownRegulated','In_cis_Increase_UpRegulated',
              'In_trans_DownRegulated','In_trans_UpRegulated'
  )
  gt <- data.frame(gene_type=gene_type,
                   Gene_Type=c('In cis Loss-Down','In cis Loss-Up',
                               'In cis Gain-Down','In cis Gain-Up',
                               'In trans Down','In trans Up'),
                   col_gt=c("#C3D7A4","#52854C",
                            "#FFDB6D","#D16103",
                            "#4E84C4","#293352"))
  
  # gt <- data.frame(gene_type=gene_type,
  #                  Gene_Type=c('Incis_DD','Incis_DU',
  #                              'Incis_ID','Incis_IU',
  #                              'Intrans_D','Intrans_U'),
  #                  col_gt=c('#E61367','#D35B53',
  #                           '#EF816E','#E82503',
  #                           '#7BEADF','#093CB4'))
  keyvals_colour <- c("In cis Loss-Down"="#C3D7A4","In cis Loss-Up"="#52854C",
                      "In cis Gain-Down"="#FFDB6D","In cis Gain-Up"="#D16103",
                      "In trans Down"="#4E84C4","In trans Up"="#293352")
  
  # keyvals_colour <- c('#E61367','#D35B53',
  #           '#EF816E','#E82503',
  #           '#7BEADF','#093CB4')
  stat <- stat %>% left_join(gt, by='gene_type')
  unique(stat$col_gt)
  # stat$DE_analysis <- gsub('(SA609_)|(SA1035_)|(SA535)','',stat$desc)
  # stat$DE_analysis <- gsub('^_','',stat$DE_analysis)
  write.csv(stat,paste0(save_dir,'summary_incis_intrans_genes_all_series.csv'), row.names = F, quote = F)
  stat <- data.table::fread(paste0(save_dir,'summary_incis_intrans_genes_all_series.csv')) %>% as.data.frame()
  # input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/cis_trans_landscape_treated_and_untreated/differential_expression/'
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
  pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
  pair_df <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  # pair_df <- read.csv(paste0(base_dir,'rnaseq_v6/incis_intrans/comparisons_drug_res.csv'), stringsAsFactors=F, check.names=F)
  # pair_df <- pair_df %>% 
  #   dplyr::select(DE_desc1, desc)
  # colnames(pair_df)
  rownames(pair_df) <- pair_df$desc
  pair_df <- pair_df %>%
    dplyr::select(-result_fn)
  # Corrected R to A
  # stat1 <- read.csv(paste0(input_dir,'summary_incis_intrans_genes_corrected.csv'), check.names = F, stringsAsFactors = F)
  # View(stat1)
  stat1 <- stat
  # View(head(stat1))
  # stat1$pct_genes <- stat1$pct_genes * 100
  stat1 <- stat1 %>% left_join(pair_df, by=c('desc','datatag'))
  
  pdx_level = c('SA501','SA530','SA604','SA609','SA535','SA1035')
  stat1$PDX <- stat1$datatag
  stat1$PDX <- factor(stat1$PDX, levels = pdx_level)
  stat1$Gene_Type <- factor(stat1$Gene_Type, levels = gt$Gene_Type)
  # stat1$DE_desc1 <- gsub('M_N_O','MNO',stat1$DE_desc1)
  # stat1$DE_desc1 <- gsub('S_T','ST',stat1$DE_desc1)
  # summary(stat1$plt_desc)
  
  # t <- data.frame(plt_desc=unique(stat1$plt_desc))
  # data.table::fwrite(t, paste0(save_dir,'plt_desc.csv'), quote = F)
  # 
  my_font <- "Helvetica"
  # stat1$desc
  # stat1$plt_desc <- NULL
  save_dir <- input_dir
  data.table::fwrite(stat1, paste0(save_dir,'summary_incis_intrans_genes_plots.csv'), quote = F)
  
  
  # Plotting
  stat1 <- data.table::fread(paste0(save_dir,'summary_incis_intrans_genes_plots.csv')) %>% as.data.frame()
  plt_desc_df <- data.table::fread(paste0(save_dir,'plt_desc_v2.csv')) %>% as.data.frame()
  dim(plt_desc_df)
  dim(stat1)
  # View(plt_desc_df)
  
  plt_desc_df$ord <- rep(1:dim(plt_desc_df)[1],1)
  # stat1$plt_desc <- NULL
  # stat1$ord <- NULL
  # stat1 <- stat1 %>% inner_join(plt_desc_df, by=c('desc'))
  
  
  table(stat1$Gene_Type,stat1$PDX)
  rownames(stat1) <- paste0(stat1$desc, stat1$Gene_Type)
  # stat2 <- stat1 %>%
  #   dplyr::filter(PDX=='SA1035')
  # table(stat2$Gene_Type, stat2$desc)
  # for(d in unique(stat1$desc)){
  #   gt <- unique(stat1[stat1$desc==d,'Gene_Type'])
  #   gs <- gt$Gene_Type[!gt$Gene_Type %in% gt]
  #   for(g in gs){
  #     stat1[stat1$desc==d,'Gene_Type']
  #   }
  #   if(is.null(stat1[r,'pct_genes'])){
  #     print('error')
  #   }
  # }
  # colnames(stat1)
  # plt_desc_df$plt_desc
  
  p <- ggplot(stat1, aes(fill=Gene_Type, y=pct_genes, x=plt_desc)) + 
       geom_bar(stat="identity",position = "dodge",width = 0.8) +  #identity position="fill" , 
       scale_x_discrete(drop = FALSE) +
       # geom_text(aes(label=pct_genes), size=3.5) + #, color="white", vjust=1.6, 
       # geom_text(aes(label=ifelse(pct_genes < 5, NA, pct_genes)) , position = position_stack(vjust = 0.5), color="white", size=2.6) + #position=position_dodge(width=0.9), 
       scale_y_continuous(labels = scales::percent_format(scale = 1)) +
       facet_grid(. ~ PDX, scales="free", space='free') + 
       scale_fill_manual(values = keyvals_colour, name='Gene Type')
  # p <- p + labs(x='DE Analysis: Resistant vs Sensitive cells', y="(%) Gene Type ", title='Differentially Expressed Genes - Copy Number Driven Proportion')
  p <- p + labs(x=NULL, y="(%) Gene Type ", title='Differentially Expressed Genes - Copy Number Driven Proportion')
  # p <- p + theme(text=element_text(family=my_font),
  #                legend.title = element_text(color="black", size=13, hjust = 0.5, family=my_font), 
  #                legend.text = element_text(color="black", size=11, hjust = 0.5, family=my_font),
  #                plot.title = element_text(color="black", size=14, hjust = 0.5, family=my_font, face = "bold"),
  #                legend.position = "bottom",
  #                # axis.line = element_blank(), 
  #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  #                panel.border = element_blank(),
  #                axis.text.x = element_text(size=12, angle = 90, color="black", family=my_font),
  #                axis.text.y = element_text(size=12, color="black", family=my_font),
  #                axis.title = element_text(size=12, color="black", family=my_font)
  #                # axis.ticks = element_blank()
  # )
  thesis_theme <- ggplot2::theme(
    text = element_text(size = 8, hjust = 0.5, family=my_font),
    axis.title.x = element_text(size=8, hjust = 0.5, family=my_font),  
    axis.title.y = element_text(size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(size=7, hjust = 0.5, family=my_font, angle = 90),  
    axis.text.y = element_text(size=7, hjust = 0.5, family=my_font),
    plot.title = element_text(size=10, face="bold", hjust=0, family=my_font),
    strip.text.x = element_text(size=9),
    strip.text.y = element_text(size=9),
    legend.title=element_text(size=7, hjust = 0.5, family=my_font), 
    legend.text=element_text(size=7, hjust = 0.5, family=my_font),
    legend.spacing.x = unit(0.1, 'mm'),
    legend.spacing.y = unit(0.1, 'mm'),
    legend.key.height=unit(1,"line"),
    legend.position = "bottom",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  )
  
  
  p <- p + thesis_theme
  output_dir <- save_dir
  png(paste0(output_dir,"incis_intrans_genes_prevalence_summary.png"), height = 2*600, width=2*650, res = 2*72)
  print(p)
  dev.off()
  ggsave(paste0(output_dir,"incis_intrans_genes_prevalence_summary.pdf"),
         plot = p,
         height = 5.5,
         width = 8,
         useDingbats=F)
  saveRDS(p, paste0(output_dir,"incis_intrans_genes_prevalence_summary.rds"))
  
  
  
  # dotplot here 
  save_dir <- input_dir
  stat1 <- data.table::fread(paste0(save_dir,'summary_incis_intrans_genes_plots.csv')) %>% as.data.frame()
  dim(stat1)
  # stat1 <- stat1 %>% left_join(gt, by='gene_type')
  # stat1$pct_genes
  # stat_backup <- stat1
  # pd <- 'SA535'
  unique(stat1$comp_type)
  meta_ptx <- data.frame(datatag=c("SA501","SA530","SA604","SA609","SA535","SA1035"), 
                       pt=paste0("Pt",rep(1:6,1)))
  stat1 <- stat1 %>% inner_join(meta_ptx, by='datatag')
  
  stat1$plt_desc <- paste0(stat1$pt,':\n',stat1$plt_desc)
  # stat1$plt_desc[1]
  # chr_lg <- max(nchar(stat1$plt_desc))
  # dc <- stat1$plt_desc
  # length(dc)
  # 
  # for(i in rep(1:length(dc),1)){
  #   if(nchar(dc[i])<chr_lg){
  #     stat1$plt_desc[i] <- paste0(paste(rep(' ',chr_lg-nchar(dc[i])), collapse=''),dc[i])
  #   }
  # }
  pdxs = c('SA609','SA535','SA1035')
  pdxs_untreated = c('SA501','SA530','SA604','SA609','SA535','SA1035')
  my_font <- "Helvetica"
  
  dim(stat1)
  # First get untreated comparison
  df_untreated <- stat1 %>% 
    dplyr::filter(PDX %in% pdxs_untreated & comp_type=='untreated')
  dim(df_untreated)
  # View(df$ord)
  unique(df_untreated$desc)
  
  df_untreated$ord <- rep(1:dim(df_untreated)[1],1)
  df_untreated <- df_untreated[order(df_untreated$ord, decreasing = T),]
  # df$datatag <- factor(df$datatag, levels = pdxs_untreated)
  df_untreated$plt_desc <- factor(df_untreated$plt_desc, levels = unique(df_untreated$plt_desc))
  # View(df)
  # unique(df$plt_desc)
  
  
  plot_ls <- list()
  df_untreated <- df_untreated %>%
    arrange(ord)
  
  # Untreated series
  plot_ls[['untreated']] <- dotplot_prevalence(df_untreated)
  res_untreated_plt <- plot_cis_pos_neg(df_untreated)
  plot_ls[['untreated_proportion']] <- res_untreated_plt$p
  res_untreated_plt$plg
  unique(stat1$comp_type)
  # Time series
  for(pd in pdxs){
    df <- stat1 %>% 
      dplyr::filter(PDX==pd & comp_type=='treated_vs_untreated')
    dim(df)
    # df$ord <- rep(1:dim(df)[1],1)
    # df <- df[order(df$ord, decreasing = T),]
    plot_ls[[pd]] <- dotplot_prevalence(df)
    res_plt <- plot_cis_pos_neg(df)
    plot_ls[[paste0(pd,'_proportion')]] <- res_plt$p
    
  }
  
  p <- ggplot(df, aes(x = pct_genes, y = plt_desc, color = Gene_Type)) + 
    geom_point(size=5, alpha=1) + #, size=2*log10(pct_genes)
    scale_color_manual(name = NULL, values = keyvals_colour) + 
    theme(legend.position ='bottom',
          legend.text = element_text(size = 9, hjust = 0.5, family=my_font),
          legend.title = element_text(size = 10, hjust = 0.5, family=my_font, angle = 90))
  p <- p + guides(color = guide_legend(override.aes = list(size=8, ncol = 2)))#shape = 0, , ncol = 2
  lg <- cowplot::get_legend(p)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  plg

  
  
  plot_ls[['lg']] <- plg
  plot_ls[['lg_cis']] <- res_untreated_plt$plg
  
  
  # untreated_row <- cowplot::plot_grid(plot_ls$untreated, plot_ls$untreated_proportion, rel_widths = c(4,1), ncol=2)
  # sa609_row <- cowplot::plot_grid(plot_ls$SA609, plot_ls$SA609_proportion, rel_widths = c(4,1), ncol=2)
  # sa535_row <- cowplot::plot_grid(plot_ls$SA535, plot_ls$SA535_proportion, rel_widths = c(4,1), ncol=2)
  # sa1035_row <- cowplot::plot_grid(plot_ls$SA1035, plot_ls$SA1035_proportion, rel_widths = c(4,1), ncol=2)
  # lgt_row <- cowplot::plot_grid(plot_ls$lg, plot_ls$lg_cis, rel_widths = c(4,1.3), ncol=2)
  # p_total <- cowplot::plot_grid(untreated_row,sa609_row,sa535_row,sa1035_row,lgt_row,
  #                               ncol = 1,align='vh', rel_heights = c(6,4,3,2,1))
  
  lgs <- cowplot::plot_grid(prop_plt$plg, pathway_plt$plg, ncol=1)
  lgs2 <- cowplot::plot_grid(plot_ls$lg, plot_ls$lg_cis, ncol=1)
  
  lgs_total <- cowplot::plot_grid(lgs2, lgs,nrow = 2)
  untreated_row <- cowplot::plot_grid(plot_ls$untreated, plot_ls$untreated_proportion, 
                                      lgs_total,
                                      rel_widths = c(2.8,1,2.5), ncol=3)
  sa609_row <- cowplot::plot_grid(plot_ls$SA609, plot_ls$SA609_proportion, 
                                  prop_plt_ls$SA609, pathway_plt_ls$SA609,
                                  rel_widths = c(2.8,1,2,0.5), ncol=4)
  sa535_row <- cowplot::plot_grid(plot_ls$SA535, plot_ls$SA535_proportion, 
                                  prop_plt_ls$SA535, pathway_plt_ls$SA535,
                                  rel_widths = c(2.8,1,2,0.5), ncol=4)
  sa1035_row <- cowplot::plot_grid(plot_ls$SA1035, plot_ls$SA1035_proportion, 
                                   prop_plt_ls$SA1035, pathway_plt_ls$SA1035,
                                   rel_widths = c(2.8,1,2,0.5), ncol=4)
  # lgt_row <- cowplot::plot_grid(plot_ls$lg, plot_ls$lg_cis, rel_widths = c(4,1,2,1), ncol=4)
  # p_total <- cowplot::plot_grid(untreated_row,sa609_row,sa535_row,sa1035_row,lgt_row,
  #                               ncol = 1,align='vh', rel_heights = c(6,4,3,2,1))
  p_total <- cowplot::plot_grid(NULL, untreated_row,NULL,sa609_row,NULL,sa535_row,NULL,sa1035_row,
                                ncol = 1,align='v', rel_heights = c(0.5,6,0.5,4,0.5,3,0.5,2),
                                labels = c('a -Rx PDX tumors','',
                                           'b -Rx/+Rx time series Pt4','',
                                           'c -Rx/+Rx time series Pt5','',
                                           'd -Rx/+Rx time series Pt6',''))
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/Fig4_prevalence_cistrans/'
  # dir.create(save_fig_dir)
  saveRDS(untreated_row, paste0(save_fig_dir,'untreated_row.rds'))
  saveRDS(sa609_row, paste0(save_fig_dir,'sa609_row.rds'))
  saveRDS(sa535_row, paste0(save_fig_dir,'sa535_row.rds'))
  saveRDS(sa1035_row, paste0(save_fig_dir,'sa1035_row.rds'))
  saveRDS(lgt_row, paste0(save_fig_dir,'lgt_row.rds'))
  saveRDS(p_total, paste0(save_fig_dir,'p_total.rds'))
  
  # p_total <- cowplot::plot_grid(plotlist=plot_ls, ncol = 1, 
  #                               align='vh', rel_heights = c(6,4,3,2,1))
  # labels = c('UnRx vs. UnRx','SA609: Res Rx vs. Sen UnRx',
  #            'SA535: Res Rx vs. Sen UnRx','SA1035: Res Rx vs. Sen UnRx')
                                # labels = c('UnRx vs. UnRx','Rx vs. UnRx','Rx vs. UnRx','Rx vs. UnRx'),
                                # label_size = 10,
                                # label_fontfamily = my_font,
                                # label_x = 0, label_y = 0
                                # )#hjust = 0, vjust = 0
  # p_total <- p_total + labs(x=NULL, y="(%) Gene Type ", title='Differentially Expressed Genes - Copy Number Driven Proportion')
  png(paste0(save_fig_dir,"Fig3_incis_intrans_genes_prevalence_summary.png"), 
      height = 2*800, width=2*750, res = 2*72)
  print(p_total)
  dev.off()
  
  ggsave(paste0(save_fig_dir,"Fig3_incis_intrans_genes_prevalence_summary.pdf"),
                  plot = p_total,
                  height = 9,
                  width = 9.8,
                  useDingbats=F,
                  dpi=300)
  
  # saveRDS(p_total, paste0(save_dir,"incis_intrans_genes_prevalence_dotplot.rds"))
  
}
dotplot_prevalence <- function(df){
  keyvals_colour <- c("In cis Loss-Down"="#C3D7A4","In cis Loss-Up"="#52854C",
                      "In cis Gain-Down"="#FFDB6D","In cis Gain-Up"="#D16103",
                      "In trans Down"="#4E84C4","In trans Up"="#293352")
  df$plt_desc <- factor(df$plt_desc, levels = unique(df$plt_desc))
  # View(df)
  # unique(df$plt_desc)
  # xmax <- max(df$pct_genes) + 5
  # if(xmax>100){
  #   xmax <- 100
  # }
  df$sz <- ifelse(df$pct_genes>0, 2*log10(df$pct_genes), 0)
  p <- ggplot(df, aes(x = pct_genes, y = plt_desc)) + 
    geom_point(aes(color = Gene_Type, size=sz), alpha=0.9) + #, size=2*log10(pct_genes), size=4.5
    scale_color_manual(values = keyvals_colour) + 
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    scale_x_continuous(trans = scales::log2_trans(),
                       breaks = scales::trans_breaks("log2", function(x) 2^x)) + 
    # scale_x_continuous(breaks = scales::log2_trans()) + 
    # facet_grid(rows = vars(PDX)) + 
    # geom_hline(yintercept = 1:6, col = "#e2e2e2") +
    # geom_text(aes(label=celltype_desc))+
    # annotate('text', x = df$success, y = df$index, label = df$success, size=3, colour = col[df$gender])+
    
    # scale_color_manual(values = col) +
    theme_bw(base_size = 10) + 
    theme(legend.position = "none",
          # panel.grid.major.x = element_blank(),
          # panel.grid.minor.x = element_blank(),
          # axis.ticks.y = element_blank(),
          # axis.text  = element_blank(),
          text = element_text(size = 8, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=7, hjust = 0.5, family=my_font),  #, angle = 90
          axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=10, face="bold", hjust=0, family=my_font),
          axis.title = element_blank()) #+
    # xlim(0,xmax)
  # geom_text_repel( data = df, aes(label = pct_genes), nudge_x=0, nudge_y=0, max.overlaps = Inf,
  #                  min.segment.length = 0,segment.size=0.1,
  #                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
  return(p)
}
plot_cis_pos_neg <- function(df){
  my_font <- "Helvetica"
  stat2 <- df %>% 
    dplyr::filter(grepl('In_cis',gene_type))%>% 
    dplyr::select(gene_type,pct_genes,desc)%>% 
    tidyr::pivot_wider(names_from = 'gene_type', values_from = 'pct_genes')
  sum(is.na(stat2))
  stat2[is.na(stat2)] <- 0
  stat2 <- stat2 %>%
    dplyr::mutate(incis_positive=In_cis_Decrease_DownRegulated+In_cis_Increase_UpRegulated,
                  incis_negative=In_cis_Decrease_UpRegulated+In_cis_Increase_DownRegulated)
  stat3 <- stat2 %>%
    dplyr::select(incis_positive,incis_negative, desc)
  # tibble::column_to_rownames(desc)%>% 
  # stat3 <- stat3 %>% tibble::column_to_rownames('desc')
  stat3 <- stat3 %>% tidyr::pivot_longer(!desc, names_to = "gene_type", values_to = "pct_genes")
  
  cis_cols <- c('lightseagreen','red')
  names(cis_cols) <- c("incis_positive","incis_negative")
  p <- ggplot(stat3, aes(fill=gene_type, y=pct_genes, x=desc)) + 
    geom_bar(position="fill", stat="identity") + 
    coord_flip() + 
    scale_fill_manual(values = cis_cols, name='',labels = c("In cis same directions", 
                                                            "In cis opposite directions")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + 
    theme(text=element_text(family=my_font),
          axis.text.x = element_text(size=9, color="black", family=my_font),
          axis.text.y = element_blank(),
          # legend.position ='bottom',
          plot.title = element_text(size=10, hjust=0.5, family=my_font),
          legend.text = element_text(size=9, color="black", family=my_font),
          legend.title = element_blank())
  p <- p + guides(fill = guide_legend(override.aes = list(shape = 0, size=6, nrow = 2)))
  lg <- cowplot::get_legend(p)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  
  p1 <- ggplot(stat3, aes(fill=gene_type, y=pct_genes, x=desc)) + 
    geom_bar(position="fill", stat="identity") + 
    coord_flip() + 
    scale_fill_manual(values = cis_cols, name='Gene Type',labels = c("In Cis Positive", "In Cis Negative")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + 
    theme(plot.title = element_text(size=10, hjust=0.5, family=my_font),
          text=element_text(family=my_font),
          legend.position = 'none',
          axis.text.x = element_text(size=9, hjust = 0.5, family=my_font),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          # axis.ticks.length=unit(0.25, "cm"),
          # axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = 'white', colour = 'white'))
  p1 <- p1 + labs(title='In cis direction')
  return(list(plg=plg, p=p1))
}

plot_all_volcano <- function(base_dir, input_dir){
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results','/SA609_rna/deg_analysis/SA609-v6/')
  p609 <- readRDS(paste0(input_dir,'SA609_UTTT_R_UUUU_C/DE_SA609_UTTT_R_UUUU_C.rds'))
  input_dir <- paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/')
  p1035 <- readRDS(paste0(input_dir,'SA1035_UTTTT_G_UUUUU_E/DE_SA1035_UTTTT_G_UUUUU_E.rds'))
  input_dir <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/')
  pSA535_cis <- readRDS(paste0(input_dir,'SA535_UUTTT_S_T_UUUUU_Q/DE_SA535_UUTTT_S_T_UUUUU_Q.rds'))
  pSA535_cx <- readRDS(paste0(input_dir,'SA535_UXXXX_U_UUUUU_Q/DE_SA535_UXXXX_U_UUUUU_Q.rds'))
  save_dir <- input_dir
  main_plot <- cowplot::plot_grid(p609, p1035, pSA535_cis, pSA535_cx,
    ncol = 2,
    nrow=2,
    # rel_heights = c(2.5,1),
    align = 'hv'
  )
  png(paste0(save_dir,"volcano_3series.png"), height = 2*920, width=2*1050, res = 2*72)
  print(main_plot)
  dev.off()
}

plot_fitness_genes <- function(){
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  input_dir <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/')
  pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
  id609 <- paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/')
  
  id1035 <- paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/')
  id535 <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/')
  pair_groups <- read.csv(pair_groups_fn, header=T, check.names=F, stringsAsFactors=F)
  View(pair_groups)
  dim(pair_groups)
  length(unique(pair_groups$desc))
  unique(pair_groups$datatag)
  pair_groups$result_dir <- ifelse(pair_groups$datatag=='SA609',paste0(id609,pair_groups$desc,'/'),
                                   ifelse(pair_groups$datatag=="SA1035",paste0(id1035,pair_groups$desc,'/'),
                                          paste0(id535,pair_groups$desc,'/')))
  
  if(is.null(cancer_ref_genes_fn)){
    ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
    cancer_ref_genes_fn <- paste0(ref_dif,'Behan_CFgenes.csv')
  }
  cancer_ref_genes_df <- read.csv(cancer_ref_genes_fn, stringsAsFactors=F, check.names = F)
  print(dim(cancer_ref_genes_df))
  colnames(cancer_ref_genes_df)[which(colnames(cancer_ref_genes_df) == "ADAM PanCancer Core-Fitness genes")] <- "PanCancer_Fitness_genes"
  
  
  cancer_genes_ls <- list()
  rownames(pair_groups) <- pair_groups$desc
  minLogFC <- 0.25
  FDR_cutoff <- 0.01
  pValueThrs <- 0.05
  for(de in pair_groups$desc){
    signif_genes <- read.csv(paste0(pair_groups[de,'result_dir'],'signif_genes.csv'), check.names = F, stringsAsFactors=F)
    signif_genes <- signif_genes[abs(signif_genes$logFC)>minLogFC & signif_genes$FDR<FDR_cutoff  & signif_genes$PValue<pValueThrs,]
    signif_genes <- signif_genes[signif_genes$is_fitness_gene,]
    print(dim(signif_genes))
    # fitness_genes <- summary(as.factor(signif_genes$classified_gene_dlp))
    nb_incis <- sum(!is.na(signif_genes$classified_gene_dlp)==TRUE)
    nb_intrans <- sum(is.na(signif_genes$classified_gene_dlp)==TRUE)
    # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
    cancer_genes_ls[[paste0(de,'_incis')]] <- round(nb_incis*100/nrow(cancer_ref_genes_df),2)
    cancer_genes_ls[[paste0(de,'_intrans')]] <- round(nb_intrans*100/nrow(cancer_ref_genes_df),2)
  }
  summary(as.factor(pair_groups$title))
  unique(pair_groups$title)
  cancer_pct_df <- do.call(rbind, cancer_genes_ls)
  cancer_pct_df <- as.data.frame(cancer_pct_df)
  dim(cancer_pct_df)
  View(head(cancer_pct_df))
  colnames(cancer_pct_df)[which(colnames(cancer_pct_df)=='V1')] <- 'Fitness_genes_pct'
  series = lapply(strsplit(rownames(cancer_pct_df), "_"), function(x) {
    return(x[1])
  })
  cancer_pct_df$PDX <- as.character(series)
  
  level_pdx <- c('SA609_CISPLATIN','SA1035_CISPLATIN','SA535_CISPLATIN','SA535_CX5461')
  cancer_pct_df$PDX <- ifelse(cancer_pct_df$PDX=='SA609','SA609_CISPLATIN',
                              ifelse(cancer_pct_df$PDX=='SA1035','SA1035_CISPLATIN',cancer_pct_df$PDX))
  cancer_pct_df$PDX <- ifelse(grepl('SA535',rownames(cancer_pct_df)) & grepl('T',rownames(cancer_pct_df)),'SA535_CISPLATIN',
                              ifelse(grepl('SA535',rownames(cancer_pct_df)) & grepl('X',rownames(cancer_pct_df)),'SA535_CX5461',cancer_pct_df$PDX))
  summary(as.factor(cancer_pct_df$PDX))
  View(rownames(cancer_pct_df))
  # cancer_pct_df$PDX <- c(rep(level_pdx, c(12, 8, 14, 12)))
  cancer_pct_df$PDX <- factor(cancer_pct_df$PDX, levels=level_pdx)
  
  
  rn <- as.character(rownames(cancer_pct_df))
  cancer_pct_df$desc <- rn #factor(rn, levels = rn)
  output_dir <- input_dir
  write.csv(cancer_pct_df, paste0(output_dir,'total_cancer_pct_df.csv'), quote=F, row.names = F)
  cancer_pct_df <- read.csv(paste0(output_dir,'total_cancer_pct_df.csv'), check.names = F, stringsAsFactors = F)
  View(cancer_pct_df)
  cancer_pct_df$gene_type <- ifelse(grepl("*incis",cancer_pct_df$desc),"incis",'intrans')
  # viz_pancancer_genes(cancer_pct_df, output_dir, tag="totaldata")
  unique(cancer_pct_df$gene_type)
  # desc <- mclapply(strsplit(cancer_pct_df$desc, "_"), function(x) {
  #   if(length(x)==7){
  #     return(paste(x[2],x[3],x[4],x[5], sep='_'))
  #   }else{
  #     return(paste(x[2],x[3],x[4],x[5],x[6], sep='_'))
  #   }
  #   
  # }, mc.cores = 2)
  cancer_pct_df$de_pair <- cancer_pct_df$desc
  # cancer_pct_df$de_pair <- gsub('_incis', '', grep("*incis",cancer_pct_df$de_pair, value=T))
  # cancer_pct_df$de_pair <- gsub('_intrans', '', grep("*intrans",cancer_pct_df$de_pair, value=T))
  cancer_pct_df$de_pair <- ifelse(grepl("*incis",cancer_pct_df$de_pair), gsub('_incis', '', cancer_pct_df$de_pair),
                                  ifelse(grepl("*intrans",cancer_pct_df$de_pair), gsub('_intrans', '', cancer_pct_df$de_pair),cancer_pct_df$de_pair))
  # View(cancer_pct_df$de_pair)
  cancer_pct_df$de_pair <- gsub('^(SA609_)|(SA1035_)|(SA535_)','',cancer_pct_df$de_pair)
  colorcode <- c('#FF3232','#40A0E0')
  ylabel <- '(%) Core-Fitness genes '
  xlabel <- 'DE Analysis: Resistant versus Sensitive cells'
  plottitle <- 'Percentage genes in ADAM PanCancer Core-Fitness genes'
  
  # cancer_pct_df$de_pair <- factor(cancer_pct_df$de_pair, levels = unique(cancer_pct_df$de_pair))
  
  # cancer_pct_df$series <- factor(cancer_pct_df$series, levels = unique(cancer_pct_df$series))
  plot_stack_barplot(cancer_pct_df, colorcode, xlabel, ylabel, plottitle, output_dir,'totaldata',
                     fa='gene_type', xa='de_pair', ya='Fitness_genes_pct')
  
 
  
  # Genes that appear in 20 cisplatin genes in 3 series
  rownames(pair_groups) <- pair_groups$desc
  minLogFC <- 0.25
  FDR_cutoff <- 0.01
  pValueThrs <- 0.05
  de_genes_ls <- list()
  for(de in pair_groups$desc){
    signif_genes <- read.csv(paste0(pair_groups[de,'result_dir'],'signif_genes.csv'), check.names = F, stringsAsFactors=F)
    signif_genes <- signif_genes[abs(signif_genes$logFC)>minLogFC & signif_genes$FDR<FDR_cutoff  & signif_genes$PValue<pValueThrs,]
    # signif_genes <- signif_genes[signif_genes$is_fitness_gene,]
    print(dim(signif_genes))
    colnames(signif_genes)
    signif_genes <- signif_genes[,c("gene_symbol","ensembl_gene_id")]
    signif_genes$datatag <- pair_groups[de,'datatag']
    dim(signif_genes)
    de_genes_ls[[paste0(de)]] <- signif_genes
  }
  summary(as.factor(pair_groups$title))
  unique(pair_groups$title)
  de_genes_df <- do.call(rbind, de_genes_ls)
  de_genes_df <- as.data.frame(de_genes_df)
  View(head(de_genes_df))
  dim(de_genes_df)
  data.table::fwrite(de_genes_df,file = paste0(input_dir,'total_DE_genes_3series.csv'), sep=',', quote=F)
  de_genes_df$description <- paste(de_genes_df$gene_symbol,de_genes_df$datatag, sep='_')
  de_genes_df1 <- de_genes_df[!duplicated(de_genes_df$description),]
  dim(de_genes_df1)
  View(head(de_genes_df1))
  colnames(de_genes_df1)
  de_genes_df2 <- de_genes_df1[,c("gene_symbol","datatag")]
  # de_genes_df2 <- de_genes_df2 %>%
  #   pivot_wider(names_from = datatag, values_from = gene_symbol)
  

  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  cis_20 <- read.csv(paste0(base_dir,'biodatabase/cisplatin_20_genes.csv'), check.names = F, stringsAsFactors=F)

  dim(cis_20)
  de_genes_20 <- de_genes_df2[de_genes_df2$gene_symbol %in% cis_20$gene,]
  dim(de_genes_20)

  rownames(cis_20) <- cis_20$gene
  cis_20$DE_expression_series <- NA
  for(g in cis_20$gene){
    c=de_genes_20[de_genes_20$gene_symbol==g,'datatag']
    print(paste(c,collapse = '_'))
    if(!is.null(c)){
      cis_20[g,'DE_expression_series'] <- paste(c,collapse = '_')
    }
  }

  write.csv(cis_20, file = paste0(input_dir,'cisplatin_20genes_3series.csv'), quote=F)
  colnames(de_genes_20)
  # rownames(de_genes_20) <- de_genes_20
  de_genes_20$gene_symbol <- as.character(de_genes_20$gene_symbol)
  de_genes_20 <- de_genes_20 %>%
      pivot_wider(names_from = gene_symbol, values_from = datatag, values_fn = list)
  View(de_genes_20)
  de_genes_20 <- t(de_genes_20)
  de_genes_20 <- as.data.frame(de_genes_20)
  colnames(de_genes_20) <- 'DE_expression_series'
  de_genes_20$gene_symbol <- rownames(de_genes_20)
  de_genes_20$DE_expression_series <- paste(de_genes_20$DE_expression_series, collapse = '_')
  cis_20
  dim(de_genes_20)
  head(cis_20)
  library(dplyr)
  cis_20 <- cis_20 %>% left_join(de_genes_20,by=c("gene"="gene_symbol"))
  View(cis_20)
  
  
  
  cosmic_genes <- read.table(paste0(base_dir, 'biodatabase/oncogene_cosmic.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
  ref_cis_genes <- read.table(paste0(base_dir, 'biodatabase/cisplatin_resistance_genes.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
  # View(head(cosmic_genes))
  # View(head(ref_cis_genes))
  dim(cosmic_genes)
  dim(ref_cis_genes)
  cancer_genes_ls <- list()
  cosmic_genes_ls <- list()
  # rownames(pair_groups) <- pair_groups$desc
  minLogFC <- 0.25
  FDR_cutoff <- 0.01
  pValueThrs <- 0.05
  for(de in pair_groups$desc){
    signif_genes <- read.csv(paste0(pair_groups[de,'result_dir'],'signif_genes.csv'), check.names = F, stringsAsFactors=F)
    signif_genes <- signif_genes[abs(signif_genes$logFC)>minLogFC & signif_genes$FDR<FDR_cutoff  & signif_genes$PValue<pValueThrs,]
    signif_genes_cosmic <- signif_genes[signif_genes$gene_symbol %in% cosmic_genes$Gene_Symbol,]
    signif_genes_cis <- signif_genes[signif_genes$gene_symbol %in% ref_cis_genes$gene_symbol,]
    print(dim(signif_genes_cosmic))
    # fitness_genes <- summary(as.factor(signif_genes$classified_gene_dlp))
    nb_incis_cosmic <- sum(!is.na(signif_genes_cosmic$classified_gene_dlp))
    nb_intrans_cosmic <- sum(is.na(signif_genes_cosmic$classified_gene_dlp))
    # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
    cosmic_genes_ls[[paste0(de,'_incis')]] <- round(nb_incis_cosmic*100/nrow(cosmic_genes),2)
    cosmic_genes_ls[[paste0(de,'_intrans')]] <- round(nb_intrans_cosmic*100/nrow(cosmic_genes),2)
    
    nb_incis_cis <- sum(!is.na(signif_genes_cis$classified_gene_dlp))
    nb_intrans_cis <- sum(is.na(signif_genes_cis$classified_gene_dlp))
    # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
    cancer_genes_ls[[paste0(de,'_incis')]] <- round(nb_incis_cis*100/nrow(ref_cis_genes),2)
    cancer_genes_ls[[paste0(de,'_intrans')]] <- round(nb_intrans_cis*100/nrow(ref_cis_genes),2)
  }
  # summary(as.factor(pair_groups$title))
  # unique(pair_groups$title)
  plot_cisplatin_resistance_genes_pct(cancer_genes_ls)
  plot_cosmic_genes_pct(cosmic_genes_ls)

  
}

plot_cisplatin_resistance_genes_pct <- function(cancer_genes_ls){
  cis_df <- do.call(rbind, cancer_genes_ls)
  cis_df <- as.data.frame(cis_df)
  colnames(cis_df)[which(colnames(cis_df)=='V1')] <- 'cisplatin_resistance_gene_pct'
  dim(cis_df)
  View(head(cis_df))
  # colnames(cis_df) <- 'cosmic_genes_pct'
  
  # View(head(cis_df))
  level_pdx <- c('SA609_CISPLATIN','SA1035_CISPLATIN','SA535_CISPLATIN','SA535_CX5461')
  pdx <- lapply(strsplit(rownames(cis_df), "_"), function(x) {
    return(x[1])
  })
  cis_df$PDX <- as.character(pdx)
  unique(cis_df$PDX)
  cis_df$PDX <- ifelse(grepl('X',rownames(cis_df)),'SA535_CX5461',
                       ifelse(!grepl('X',rownames(cis_df)),paste0(cis_df$PDX,'_CISPLATIN'),cis_df$PDX))
  summary(as.factor(cis_df$PDX))
  # cis_df$PDX <- c(rep(level_pdx, c(12, 8, 14, 12)))
  cis_df$PDX <- factor(cis_df$PDX, levels=level_pdx)
  
  colnames(cancer_pct_df)[which(colnames(cancer_pct_df)=='V1')] <- 'Fitness_genes_pct'
  cis_df$desc <- as.character(rownames(cis_df)) #factor(rn, levels = rn)
  output_dir <- input_dir
  # View(head(cis_df))
  write.csv(cis_df, paste0(output_dir,'reference_cisplatin_resistance_gene_pct.csv'), quote=F, row.names = F)
 
  cis_df$gene_type <- ifelse(grepl("*incis",cis_df$desc),"incis",'intrans')
  # viz_pancancer_genes(cancer_pct_df, output_dir, tag="totaldata")
  summary(as.factor(cis_df$gene_type))

  cis_df$de_pair <- cis_df$desc
  # cancer_pct_df$de_pair <- gsub('_incis', '', grep("*incis",cancer_pct_df$de_pair, value=T))
  # cancer_pct_df$de_pair <- gsub('_intrans', '', grep("*intrans",cancer_pct_df$de_pair, value=T))
  cis_df$de_pair <- ifelse(grepl("*incis",cis_df$de_pair), gsub('_incis', '', cis_df$de_pair),
                           ifelse(grepl("*intrans",cis_df$de_pair), gsub('_intrans', '', cis_df$de_pair),cis_df$de_pair))
  
  cis_df$de_pair <- gsub('^(SA609_)|(SA1035_)|(SA535_)','',cis_df$de_pair)
  colorcode <- c('#FF3232','#40A0E0')
  ylabel <- ' (%) Ref cisplatin resistance genes'
  xlabel <- 'DE Analysis: Resistant versus Sensitive cells'
  plottitle <- 'Percentage DE genes in reference cisplatin resistance genes'
  View(head(cis_df))
  # cancer_pct_df$de_pair <- factor(cancer_pct_df$de_pair, levels = unique(cancer_pct_df$de_pair))
  cis_df$cisplatin_genes_pct
  # cancer_pct_df$series <- factor(cancer_pct_df$series, levels = unique(cancer_pct_df$series))
  plot_stack_barplot(cis_df, colorcode, xlabel, ylabel, plottitle, output_dir,'ref_cisplatin_resistance_genes_3series',
                     fa='gene_type', xa='de_pair', ya='cisplatin_resistance_gene_pct')
  
}  
plot_cosmic_genes_pct <- function(cosmic_genes_ls){
  cosmic_df <- do.call(rbind, cosmic_genes_ls)
  cosmic_df <- as.data.frame(cosmic_df)
  dim(cosmic_df)
  colnames(cosmic_df)[which(colnames(cosmic_df)=='V1')] <- 'cosmic_genes_pct'
  level_pdx <- c('SA609_CISPLATIN','SA1035_CISPLATIN','SA535_CISPLATIN','SA535_CX5461')
  pdx <- lapply(strsplit(rownames(cosmic_df), "_"), function(x) {
    return(x[1])
  })
  cosmic_df$PDX <- as.character(pdx)
  unique(cosmic_df$PDX)
  cosmic_df$PDX <- ifelse(grepl('X',rownames(cosmic_df)),'SA535_CX5461',
                          ifelse(!grepl('X',rownames(cosmic_df)),paste0(cosmic_df$PDX,'_CISPLATIN'),cosmic_df$PDX))
  summary(as.factor(cosmic_df$PDX))
  # cis_df$PDX <- c(rep(level_pdx, c(12, 8, 14, 12)))
  cosmic_df$PDX <- factor(cosmic_df$PDX, levels=level_pdx)
  
  
  cosmic_df$desc <- as.character(rownames(cosmic_df)) #factor(rn, levels = rn)
  output_dir <- input_dir
  View(head(cosmic_df))
  write.csv(cosmic_df, paste0(output_dir,'reference_cosmic_gene_pct.csv'), quote=F, row.names = F)
  cosmic_df$de_pair <- cosmic_df$desc
  # cancer_pct_df$de_pair <- gsub('_incis', '', grep("*incis",cancer_pct_df$de_pair, value=T))
  # cancer_pct_df$de_pair <- gsub('_intrans', '', grep("*intrans",cancer_pct_df$de_pair, value=T))
  cosmic_df$de_pair <- ifelse(grepl("*incis",cosmic_df$de_pair), gsub('_incis', '', cosmic_df$de_pair),
                              ifelse(grepl("*intrans",cosmic_df$de_pair), gsub('_intrans', '', cosmic_df$de_pair),cosmic_df$de_pair))
  cosmic_df$gene_type <- ifelse(grepl("incis",cosmic_df$desc),"incis",'intrans')
  cosmic_df$de_pair <- gsub('^(SA609_)|(SA1035_)|(SA535_)','',cosmic_df$de_pair)
  colorcode <- c('#FF3232','#40A0E0')
  ylabel <- '(%) Reference cosmic genes '
  xlabel <- 'DE Analysis: Resistant versus Sensitive cells'
  plottitle <- 'Percentage DE genes in COSMIC reference genes'
  
  # cancer_pct_df$de_pair <- factor(cancer_pct_df$de_pair, levels = unique(cancer_pct_df$de_pair))
  # cancer_pct_df$series <- factor(cancer_pct_df$series, levels = unique(cancer_pct_df$series))
  plot_stack_barplot(cosmic_df, colorcode, xlabel, ylabel, plottitle, output_dir,'ref_cosmic_genes_in_3series',
                     fa='gene_type', xa='de_pair', ya='cosmic_genes_pct')
  
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

plot_dlp_prevalence_v2 <- function(output_dir=NULL){
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/SA535-v6/'
  base_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/')
  cnv_535_fn <- paste0(base_dir, 'SA535_total_rna_v2/clonealign/whole_data/combined_clones/SA535_total_cnv_mat.csv')
  cnv_535 <- read.csv(cnv_535_fn, check.names = F, stringsAsFactors = F)
  head(cnv_535)
  cnv_1035_fn <- paste0(base_dir, 'SA1035_rna/clonealign/whole_data/SA1035_cnv_mat.csv')
  cnv_1035 <- read.csv(cnv_1035_fn, check.names = F, stringsAsFactors = F)
  
  cnv_609_fn <- paste0(base_dir, 'SA609_rna/added_segments/clonealign/whole_data/SA609_cnv_mat.csv')
  cnv_609 <- read.csv(cnv_609_fn, check.names = F, stringsAsFactors = F)
  
  
  de_df <- read.csv(paste0(base_dir, 'SA535_total_rna_v2/SA535-v6/comparisons_drug_res.csv'), check.names = F, stringsAsFactors = F)
  dim(de_df)  
  rownames(de_df) <- de_df$desc
  
  # de <- de_df$desc[8]
  cn_change <- list()
  for(de in de_df$desc){
    datatag <- de_df[de,'datatag']
    if(datatag=='SA609'){
      df_cnv <- cnv_609
    }else if(datatag=='SA535'){
      df_cnv <- cnv_535
    }else{
      df_cnv <- cnv_1035
    }
    df_cnv <- df_cnv %>% 
      dplyr::select(de_df[de,'clone1'],de_df[de,'clone2'])
    
    print(dim(df_cnv))
    df_cnv_var <- df_cnv[abs(df_cnv[,de_df[de,'clone1']]-df_cnv[,de_df[de,'clone2']])>=1,]
    print(df_cnv_var)
    # df_cnv_var <- df_cnv %>%
    #   dplyr::filter(apply(df_cnv, 1, var)>0)
    # dim(df_cnv)
    # dim(df_cnv_var)
    cn_change[[as.character(de)]] <- c(nrow(df_cnv_var),nrow(df_cnv))
  }
  
  cn_change_df <- as.data.frame(t(dplyr::bind_rows(cn_change)))
  colnames(cn_change_df) <- c('CN_change','total_CNs')
  cn_change_df$desc <- rownames(cn_change_df)
  write.csv(cn_change_df, paste0(output_dir,'cn_change.csv'), row.names = F, quote = F)
  View(cn_change_df)
  pdx <- lapply(strsplit(rownames(cn_change_df), "_"), function(x) {
    return(x[1])
  })
  cn_change_df$PDX <- as.character(pdx)
  cn_change_df$pct_change <- round(cn_change_df$CN_change/cn_change_df$total_CNs * 100,2)
  cn_change_df$PDX <- factor(cn_change_df$PDX, levels=c('SA609', 'SA1035','SA535'))
  p<-ggplot(data=cn_change_df, aes(x=desc, y=pct_change)) +
    geom_bar(stat="identity", fill="steelblue")+
    ggplot2::facet_grid(. ~ PDX, scales="free", space='free') +
    ggplot2::theme(legend.title = ggplot2::element_text(color="black", size=11),
                            legend.text = ggplot2::element_text(color="black", size=10),
                            plot.title = ggplot2::element_text(color="black", size=14, hjust = 0.5, face='bold'),
                            # legend.position = "none", 
                            # axis.line = element_blank(), 
                            panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
                            panel.border = ggplot2::element_blank(),
                            # axis.text.x =ggplot2::element_blank(),
                            axis.text.x = ggplot2::element_text(color="black", size=10, hjust=0.5, angle = 90),
                            axis.text.y = ggplot2::element_text(color="black", size=8, hjust = 0.5),
                            axis.title.y = ggplot2::element_text(color="black", size=10, hjust = 0.5),
                            axis.title.x = ggplot2::element_text(color="black", size=12, hjust = 0.5, face='bold')
                            # axis.title.x = ggplot2::element_blank()
                            # axis.ticks = element_blank()
    )
  p
  
  
  
}  


plot_DE_genes_ggplot <- function(df, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                 plttitle="A versus B", save_dir="",legendVisible=F,
                                 iscaption=TRUE, legend_verbose='none', save_plot=TRUE, 
                                 xl=NULL, yl=NULL){  
  # library(ggplot2)
  # library(ggrepel)
  #legend_verbose is none or 'right', 'left','bottom'
  # df <- de_genes
  # library(EnhancedVolcano)
  
  # colnames(df)[which(names(df) == "avg_logFC")] <- "log2FoldChange"
  # colnames(df)[which(names(df) == "p_val_adj")] <- "padj"
  # capstr <- paste0("FDR cutoff, ",FDRcutoff,"; logFC cutoff, ",logFCcutoff, "; nb genes signif, ",nbgenessig)
  # summary(as.factor(df$gene_type))
  Gene_Type=c('In_cis_Decrease_DownRegulated','In_cis_Decrease_UpRegulated',
              'In_cis_Increase_DownRegulated','In_cis_Increase_UpRegulated',
              'In_trans_DownRegulated','In_trans_UpRegulated'
  )
  gt <- data.frame(Gene_Type=Gene_Type,
                   gene_type=c('In cis Loss-Down','In cis Loss-Up',
                               'In cis Gain-Down','In cis Gain-Up',
                               'In trans Down','In trans Up'),
                   gene_type_classify=c('In cis positive tendency','In cis negative tendency',
                                        'In cis negative tendency','In cis positive tendency',
                                        'In trans','In trans'),
                   col_gt=c("#C3D7A4","#52854C",
                            "#FFDB6D","#D16103",
                            "#4E84C4","#293352"))
  
  df <- df %>% left_join(gt, by='Gene_Type')
  # unique(df$col_gt)
  # keyvals_colour <- as.character(df$col_gt)
  # names(keyvals_colour) <- df$gene_type
  # keyvals_colour <- gt$col_gt
  # names(keyvals_colour) <- gt$gene_type
  
  keyvals_colour <- c("In cis Loss-Down"="#C3D7A4","In cis Loss-Up"="#52854C",
                      "In cis Gain-Down"="#FFDB6D","In cis Gain-Up"="#D16103",
                      "In trans Down"="#4E84C4","In trans Up"="#293352")
  keyvals_colour_classify <- c("In cis positive tendency"="#6a0dad",
                               "In cis negative tendency"="#A8A8A8",
                               "In trans"="#000000")
  # keyvals_colour <- factor(keyvals_colour, levels = unique(keyvals_colour))
  # unique(keyvals_colour)
  # names(keyvals_colour[1:3])
  df <- df[abs(df$logFC)>logFCcutoff & df$FDR<FDRcutoff  & df$PValue<pValuecutoff,]
  maxLogFC <- 3.5
  # df$logFC <- ifelse(df$logFC>maxLogFC,maxLogFC,
  #                    ifelse(df$logFC<-maxLogFC,-maxLogFC,df$logFC))
  df$logFC <- sapply(df$logFC, function(x) replace(x, x > maxLogFC, maxLogFC))
  df$logFC <- sapply(df$logFC, function(x) replace(x, x < (-maxLogFC), -maxLogFC))
  
  st <- summary(df$gene_type)
  # keyvals_shape <- ifelse(df$is_fitness_gene==T, 3, 16)
  # names(keyvals_shape) <- ifelse(df$is_fitness_gene==T,'Fitness genes','Others')
  
  # df$gene_type <- factor(df$gene_type, levels = unique(df$gene_type))
  # if(capstr=='' & iscaption){
  #   # capstr <- paste0(capstr,'With abs(logFC)>0.25, FDR<0.01, pValue<0.05 \n')
  #   capstr <- paste0(capstr,names(st[1]),':',as.numeric(st[1]), ', ')
  #   capstr <- paste0(capstr,names(st[2]),':',as.numeric(st[2]), ' \n')
  #   capstr <- paste0(capstr,names(st[3]),':',as.numeric(st[3]), ', ')
  #   capstr <- paste0(capstr,names(st[4]),':',as.numeric(st[4]), ' \n')
  #   capstr <- paste0(capstr,names(st[5]),':',as.numeric(st[5]), ', ')
  #   capstr <- paste0(capstr,names(st[6]),':',as.numeric(st[6]), ' ')
  #   # for(i in rep(1:length(st),1)){
  #   #   # print(st[i])
  #   #   capstr <- paste0(capstr,names(st[i]),':',as.numeric(st[i]), ' ')
  #   # }
  # 
  # }
  df$mlog10FDR <- -log10(df$FDR)
  
  df$mlog10FDR <- sapply(df$mlog10FDR, function(x) replace(x, is.infinite(x), 300))
  # df$gt_alpha <- ifelse(df$gene_type=='In trans Up', 0.8,
  #                       ifelse(df$gene_type=='In cis Gain-Up',0.8,0.65))  
  # df$gt_alpha <- 0.8
  # df$gt_alpha <- ifelse(df$gene_type=='In trans Up' | df$gene_type=='In cis Gain-Up', 0.8, 0.9)  
  
  my_font <- "Helvetica"
  
 
  thesis_theme <- theme(  text=element_text(size = 8,family=my_font),
                         plot.title = element_text(color="black", size=11, hjust = 0, face = "bold", family=my_font),
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                         axis.text.y = element_text(size=7, hjust = 0.5, family=my_font),
                         axis.text.x = element_text(size=7, hjust = 0.5, family=my_font),
                         axis.title = element_text(size=8, hjust = 0.5, family=my_font),
                         plot.caption = element_text(size=8, hjust = 1, family=my_font),
                         legend.title = element_text(size=7, hjust = 0.5, family=my_font, angle=90),
                         legend.text = element_text(size=7, hjust = 0, family=my_font),
                         strip.text.x = element_text(color="black",size=9, family=my_font),
                         strip.text.y = element_text(color="black",size=9, family=my_font),
                         legend.spacing.x = unit(0.1, 'mm'),
                         legend.spacing.y = unit(0.1, 'mm'),
                         legend.key.height=unit(1,"line"),
                         # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.position = legend_verbose)
  df_cis <- df %>%
    dplyr::filter(grepl('In cis',gene_type))

  topGenes <- get_top_genes(df_cis, minLogFC=0.25, nbtopup=10, nbtopdown=10)
  p_cis <- ggplot(df_cis, aes(x = logFC, y = mlog10FDR)) +
    geom_point(aes(color = gene_type), size=2.5, alpha=0.6) +
    scale_color_manual(values = keyvals_colour[1:4], name = "Gene Type") +
    thesis_theme + 
    geom_text_repel(family = my_font,
                    data = df_cis[df_cis$gene_symbol %in% topGenes, ], 
                    aes(label = gene_symbol), size = 3.5, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.3, "lines"),
                    max.overlaps = Inf,
                    min.segment.length = 0,  # draw segment lines, not matter how short they are)
                    color='black', segment.alpha = 0.2)
  p_cis <- p_cis + labs(x= bquote(~Log[2] ~ ' Fold Change'), 
                        y=bquote(~-Log[10]~italic(FDR)),title = paste0(plttitle,', ','in cis genes'))#caption = capstr
  
  if(!is.null(xl)){
    p_cis <- p_cis + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p_cis <- p_cis + ylim(yl[1],yl[2])
  }
  
  
  df_trans <- df %>%
    dplyr::filter(grepl('In trans',gene_type))
  topGenes <- get_top_genes(df_trans, minLogFC=0.25, nbtopup=10, nbtopdown=10)
  p_trans <- ggplot(df_trans, aes(x = logFC, y = mlog10FDR)) +
    geom_point(aes(color = gene_type), size=2.5, alpha=0.6) +
    scale_color_manual(values = keyvals_colour[5:6], name = "Gene Type") +
    thesis_theme + 
    geom_text_repel(family = my_font,
                    data = df_trans[df_trans$gene_symbol %in% topGenes, ], 
                    aes(label = gene_symbol), size = 3.5, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.3, "lines"),
                    max.overlaps = Inf,
                    min.segment.length = 0,  # draw segment lines, not matter how short they are)
                    color='black', segment.alpha = 0.2)
  p_trans <- p_trans + labs(x= bquote(~Log[2] ~ ' Fold Change'), 
                            y=bquote(~-Log[10]~italic(FDR)),title = paste0(plttitle,', ','in trans genes'))#caption = capstr
  if(!is.null(xl)){
    p_trans <- p_trans + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p_trans <- p_trans + ylim(yl[1],yl[2])
  }
  
  # if(legend_verbose!='none'){
  #   return(p)
  # }
  
  # xdens <- cowplot::axis_canvas(p, axis = "x") +
  #   stat_density(data = df, aes(x = logFC, group = gene_type_classify,color =gene_type_classify),
  #                alpha = 1, size = 1, position="identity",geom="line") +
  #   scale_color_manual(values = keyvals_colour_classify, name = "Gene Classify") 
  
  
  # xdens <- cowplot::axis_canvas(p, axis = "x") +
  #          # geom_jitter(data = df, aes(x = logFC, fill = gene_type_classify), 
  #          #             shape=16, position=position_jitter(0.2)) + 
  #          geom_boxplot(data = df, aes(x = logFC, color = gene_type_classify),
  #                alpha = 0.6, size = 0.3) +
  #         # scale_color_grey()
  #         scale_color_manual(values = keyvals_colour_classify, name = "Gene Classify")
  #         
  # 
  # p_main <- cowplot::insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  # 
  
  # Just to get legend
  # t <- ggplot(df, aes(x = logFC, color = gene_type_classify)) +
  #      geom_boxplot(alpha = 0.6, size = 0.3) +
  #      scale_color_manual(values = keyvals_colour_classify, name = "Gene Classify") + 
  #      theme(text=element_text(size = 8,family=my_font),
  #           legend.title = element_text(size=7, hjust = 0.5, family=my_font),
  #           legend.text = element_text(size=7, hjust = 0, family=my_font),
  #           strip.text.x = element_text(color="black",size=9, family=my_font),
  #           strip.text.y = element_text(color="black",size=9, family=my_font),
  #           legend.spacing.x = unit(0.1, 'mm'),
  #           legend.spacing.y = unit(0.1, 'mm'),
  #           legend.key.height=unit(1,"line"))
  # lg <- cowplot::get_legend(t)
  # plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  
  # blank_theme <- theme_minimal()+
  #   theme(
  #     axis.title.x = element_blank(),
  #     axis.title.y = element_blank(),
  #     panel.border = element_blank(),
  #     panel.grid=element_blank(),
  #     axis.ticks = element_blank(),
  #     plot.title=element_text(size=14, face="bold")
  #   )
  # prop_df <- df %>%
  #   dplyr::group_by(gene_type)%>%
  #   dplyr::summarise(count=n())%>%
  #   as.data.frame()%>%
  #   dplyr::rename(label=gene_type)
  
  # For plotting
  shot_lbs <- c('LD','LU','GD','GU','Down','Up')
  names(shot_lbs) <- names(keyvals_colour)
  prop_cis <- df %>%
    dplyr::filter(grepl('In cis',gene_type)) %>%
    dplyr::group_by(gene_type)%>%
    dplyr::summarise(count=n())%>%
    as.data.frame()#%>%
    # dplyr::rename(label1=gene_type)
  prop_cis$label <- shot_lbs[prop_cis$gene_type]
  
  prop_trans <- df %>%
    dplyr::filter(grepl('In trans',gene_type)) %>%
    dplyr::group_by(gene_type)%>%
    dplyr::summarise(count=n())%>%
    as.data.frame()#%>%
    # dplyr::rename(label1=gene_type)
  prop_trans$label <- shot_lbs[prop_trans$gene_type]
  # dim(prop_cis)
  names(keyvals_colour) <- shot_lbs[names(keyvals_colour)]
  pcis <- viz_proportion(prop_cis, plot_label=T, keyvals_colour)
  ptrans <- viz_proportion(prop_trans, plot_label=T, keyvals_colour)
  # p1
  
  # p1 <- ggplot(df, aes(x="", y=gene_type, fill=gene_type))+
  #   geom_bar(width = 1, stat = "identity", alpha=0.8) +
  #   coord_polar("y", start=0) + 
  #   scale_fill_manual(values=keyvals_colour) + 
  #   blank_theme + 
  #   theme(axis.text.x=element_blank(), legend.position = "none")
  # 
  # caption = capstr
  # plg <- ggplot() +
  #   annotate("text", x = 8, y = 20, color="black", size=4.5, family=my_font, label = capstr) +
  #   theme_void()
  # p2
  # p1 <- cowplot::plot_grid(p2, p1, ncol=2, rel_widths = c(2,1))
  # p_total <- cowplot::plot_grid(p, p1, rel_heights = c(1,0.2), ncol=1)
  # p_total
  # tag <- ifelse(legend_verbose=='none','','_with_legend')
  # png(paste0(save_dir,"DE_",gsub(' ','_',plttitle),tag,".png"), height = 2*650, width=2*500, res = 2*72)
  # print(p_total)
  # dev.off()
  
  
  # if(save_plot){
  #   plttitle <- gsub(':','_',plttitle)
  #   plttitle <- gsub(' ','_',plttitle)
  #   saveRDS(p_total, file=paste0(save_dir,"DE_",plttitle,".rds"))
  #   ggsave(paste0(save_dir,"DE_",gsub(' ','_',plttitle),tag,".pdf"),
  #          plot = p_total,
  #          height = 5.5,
  #          width = 6.5,
  #          useDingbats=F)
  # }
  return(list(p_vol_cis=p_cis, p_vol_trans=p_trans, p_prop_cis=pcis, p_prop_trans=ptrans))
}


viz_legend <- function(df, keyvals_colour){
  p1 <- ggplot(df, aes(x = logFC, y = mlog10FDR)) +
    geom_point(aes(color = gene_type), size=2.5,alpha=0.8) +
    scale_color_manual(values = keyvals_colour, name = "Gene Type") + 
    theme(legend.position ='bottom')
  p1 <- p1 + guides(color = guide_legend(override.aes = list(size=4, nrow=3)))
  lg <- cowplot::get_legend(p1)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  return(plg)
}
# Donut plot
# Need count, label
viz_proportion <- function(data, plot_label=F, cols=NULL){
  data <- data %>%
    dplyr::filter(count>0)
  # if(is.null(cols)){
  #   cols blah blah
  # }
  
  # Compute percentages
  # data$label <- factor(data$label, levels = names(cols))
  rownames(data) <- data$label
  data <- data[names(cols)[names(cols) %in% as.character(data$label)],]
  data$fraction <- data$count / sum(data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax <- cumsum(data$fraction)
  
  # Compute the bottom of each rectangle
  data$ymin <- c(0, head(data$ymax, n=-1))
  
  # Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2
  my_font <- "Helvetica"
  # Compute a good label
  # data$label <- data$category
  # data$label <- c('1Rx','2Rx','3Rx','1RxH')
  cols_use <- cols[as.character(unique(data$label))]
  data$label_desc <- paste0(data$label, "\n (", data$count,')') #\n: N=
  if(plot_label==F){
    p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
      geom_rect() +
      # geom_label( x=3.5, aes(y=labelPosition, label=label), size=3.1, label.size=0) +
      scale_fill_manual(values = cols_use)   
  }else{
    p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
      geom_rect() +
      # geom_label(x=5.5, aes(y=labelPosition, label=label), size=3, label.size=0) +
      geom_text( x=5, aes(y=labelPosition, label=label_desc, color=label), size=3.5, family=my_font) +
      scale_fill_manual(values = cols_use) + 
      scale_color_manual(values = cols_use) 
  }
  p <- p + coord_polar(theta="y") +
    xlim(c(1.5, 5.5)) +
    theme_void() +
    theme(legend.position = "none",
          text=element_text(family=my_font))
  # plot.margin = margin(0, 0, 0, 0, "cm"),
  # axis.ticks.length = unit(0, "null")
  # panel.background = element_rect(fill='transparent'),
  # plot.background = element_rect(fill='transparent', color=NA),
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank()
  # p
  return(p)
  
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
  View(pair_df)
  dim(meta)
  de <- de[!duplicated(de$desc),]
  de <- de %>%
    dplyr::select(datatag, desc)
  de$de_desc <- gsub(paste0(datatag,'_'),'',de$desc)
  # View( head(de))
}


get_gene_type_edgeR_v2 <- function(pair_groups, cnv_mat, 
                                  datatag, input_dir, save_dir,
                                  cancer_ref_genes_fn=NULL, 
                                  FDR_cutoff=0.01, minLogFC=0.25, pValueThrs=0.05){  #FDR_cutoff=0.01, FDR_cutoff=0.05
  
  
  # if(is.null(cancer_ref_genes_fn)){
  #   ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  #   cancer_ref_genes_fn <- paste0(ref_dif,'Behan_CFgenes.csv')
  # }
  # cancer_ref_genes_df <- read.csv(cancer_ref_genes_fn, stringsAsFactors=F, check.names = F)
  # print(dim(cancer_ref_genes_df))
  # colnames(cancer_ref_genes_df)[which(colnames(cancer_ref_genes_df) == "ADAM PanCancer Core-Fitness genes")] <- "PanCancer_Fitness_genes"
  
  if (!file.exists(save_dir)){
    dir.create(save_dir, recursive=T)
  }
  
  # pair_groups <- read.csv(pair_groups_fn, header=T, check.names=F, stringsAsFactors=F)
  # pair_groups <- pair_groups[pair_groups$datatag==datatag,]
  
  # View(pair_groups1)
  # if(!is.null(subtag)){
  #   pair_groups <- pair_groups[pair_groups$title==subtag,]
  #   # if(subtag=='SA535_cisplatin'){
  #   #   pair_groups$desc[6] <- paste0(pair_groups$desc[6],'_2')
  #   # }
  # }
  # print(dim(pair_groups))
  # print(pair_groups)
  # minLogFC <- 0.5
  
  rownames(pair_groups) <- pair_groups$desc
  frac_dlp <- list()
  genes_summary_stat <- list()
  genes_summary <- list()
  c <- 0
  for(de in pair_groups$desc){
    output_dir <- paste0(save_dir,de,"/")
    if (!file.exists(output_dir)){
      dir.create(output_dir, recursive=T)
    }
    # signif_genes <- read.csv(paste0(output_dir,'edgeR_significant_genes.csv'), check.names=F, stringsAsFactors=F)
    signif_genes <- read.csv(paste0(input_dir,pair_groups[de,'result_fn']), check.names=F, stringsAsFactors=F)
    
    print(dim(signif_genes))
    # colnames(signif_genes)
    signif_genes <- signif_genes %>%
      dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
    print('After genes filtering: ')
    print(dim(signif_genes))
    # colnames(signif_genes)[which(names(signif_genes)=='gene_id')] <- 'ensembl_gene_id'
    # print(colnames(signif_genes))
    # signif_genes <- signif_genes[abs(signif_genes$logFC)>0.25,]
    # signif_genes <- signif_genes %>% left_join(cancer_ref_genes_df, by = c("ensembl_gene_id"="ensemble_id"))
    s1 <- as.character(pair_groups[de,'clone1'])
    s2 <- as.character(pair_groups[de,'clone2'])
    if(sum(colnames(cnv_mat) %in% c(s1,s2))!=2){
      print(de)
      stop("DEBUG: do not exist 2 clones in cnv median CN profiles")
    }
    
    cnv_mat <- cnv_mat[rownames(cnv_mat) %in% signif_genes$ensembl_gene_id,]
    cnv_mat_tmp <- cnv_mat[,c(s1, s2)]
    print(paste0('Raw cnv: ',dim(cnv_mat_tmp)[1]))
    # View(cnv_mat_tmp[1:100,])
    # View(head(cnv_mat_tmp))
    cnv_mat_tmp <- cnv_mat_tmp[rowSums(is.na(cnv_mat_tmp))==0,]
    print(paste0('Removed NA values from cnv: ',dim(cnv_mat_tmp)[1]))
    
    var_genes <- apply(cnv_mat_tmp, 1, var)
    cnv_mat_tmp <- cnv_mat_tmp[var_genes > 0,]
    print(paste0('Nb bins with CN change: ',dim(cnv_mat_tmp)[1]))
    head(cnv_mat_tmp)
    t <- cnv_mat_tmp[,s1]-cnv_mat_tmp[,s2]
    cnv_mat_tmp$classified_gene_dlp <- ifelse(t>0, 'Increase', 
                                              ifelse(t<0, 'Decrease','No_variance'))
    # for(g in rownames(cnv_mat_tmp)){
    #   t <- cnv_mat_tmp[g,s1]-cnv_mat_tmp[g,s2]
    #   if(t > 0){
    #     cnv_mat_tmp[g,'classified_gene_dlp'] <- 'Increase'
    #   }else if(t < 0){
    #     cnv_mat_tmp[g,'classified_gene_dlp'] <- 'Decrease'
    #   }else{
    #     cnv_mat_tmp[g,'classified_gene_dlp'] <- 'No_variance'
    #   }
    # }  
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
    
    # signif_genes_outliers <- signif_genes[abs(signif_genes$logFC)>outlier_FC_thrs,]
    # if(nrow(signif_genes_outliers)>0){
    #   # write.csv(signif_genes_outliers, paste0(output_dir,'signif_genes_outliers.csv'), quote = F, row.names = F)
    # }
    # signif_genes <- signif_genes[abs(signif_genes$log2FoldChange) <= outlier_FC_thrs,]
    
    
    signif_genes$classified_gene_rna <- ifelse(signif_genes$logFC>0,'UpRegulated',
                                               ifelse(signif_genes$logFC<0,'DownRegulated','Invalid'))
    signif_genes$Gene_Type <- ifelse(!is.na(signif_genes$classified_gene_dlp),
                                     paste0('In_cis_',signif_genes$classified_gene_dlp, '_',signif_genes$classified_gene_rna),
                                     paste0('In_trans_',signif_genes$classified_gene_rna))
    
    signif_genes$Gene_Type <- as.factor(signif_genes$Gene_Type)
    
    # signif_genes$is_fitness_gene <- as.factor(!is.na(signif_genes$PanCancer_Fitness_genes))
    
    # plot_DE(de, signif_genes, pair_groups, output_dir,
    #         minLogFC=0.25, pValueThrs=0.05, nbtopup=30, nbtopdown=30)
    
    # observed_signif_genes <- signif_genes[abs(signif_genes$logFC)>0.25 
    #                                       & signif_genes$FDR<0.01
    #                                       & signif_genes$PValue<0.05,]
    observed_signif_genes <- signif_genes
    print(paste0('Remove nb genes with logFC < 0.25: ',sum(signif_genes$logFC<=0.25)))
    
    # print(paste0('Comparison: ',de,"  Number of fitness genes: "))
    # print(summary(as.factor(observed_signif_genes$is_fitness_gene)))
    
    # Summary incis genes
    # nb_incr_incis <- sum(observed_signif_genes$classified_gene_dlp=='Increase', na.rm = T)
    # nb_decr_incis <- sum(observed_signif_genes$classified_gene_dlp=='Decrease', na.rm = T)
    # frac_dlp[[as.character(de)]] <- c(pair_groups[de,'result_fn'],
    #                                   paste0(datatag,'-',s1,'_vs_',s2),
    #                                   nrow(cnv_mat_tmp),length(in_cis),
    #                                   nb_incr_incis,
    #                                   nb_decr_incis,
    #                                   round(length(in_cis)/nrow(cnv_mat_tmp) * 100,2))
    # 
    # st <- summary(as.factor(observed_signif_genes$Gene_Type))
    # print(st)
    # s <- as.data.frame(st)
    # # summary(-log10(signif_genes$pvalue))
    # # summary(signif_genes$padj)
    # colnames(s)[1] <- 'nb_genes'
    # # View(s)
    # s$pct_genes <- round(s$nb_genes/colSums(s), 3) * 100
    # s$datatag <- pair_groups[de,'datatag']
    # s$desc <- as.character(de)
    # s$result_fn <- pair_groups[de,'result_fn']
    # s$de_analysis <- paste0(datatag,'-',s1,'_vs_',s2)
    # # colnames(s)[which(colnames(s) == "pct_genes")] 
    # # colnames(s)[which(colnames(s) == "nb_genes")] <- paste0(datatag,'-',s1,'_vs_',s2,'_nbgenes')
    # c <- c + 1
    # s$gene_type <- rownames(s)
    # genes_summary_stat[[c]] <- s
    
    # write.csv(s, paste0(output_dir,de,'_FDR_',FDR_cutoff,'.csv'), quote = F, row.names = F)
    data.table::fwrite(signif_genes, paste0(output_dir, 'signif_genes_FDR_',FDR_cutoff,'_minLogFC_',minLogFC,'.csv'), quote = F, row.names = F)
    # genes_summary[[paste0(de,'_signif_genes')]] <- signif_genes
    # genes_summary[[paste0(de,'_stat')]] <- s
    
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
  return(stat)
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

