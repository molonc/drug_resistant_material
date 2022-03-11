library(parallel)
# library(snpEnrichment) # mclapply2
# TO DO 
# + Add 530 and 501 at the beginning
# Get CNV and segment CNV for each clone, 530 and 501 first
# + Load data, run cis trans
# + Save output in same folder 
# + Plot output, keep csv 
# + Plot statistical values 

source("/home/htran/Projects/farhia_project/rnaseq/cis_trans/in_cis_trans_utils.R")

get_info <- function(cell_ids, idx=1) {
  labels <- lapply(strsplit(cell_ids, "_"), function(x) {
    return(x[idx])
  })
  return(as.character(labels))
}

# data_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'

create_metasample <- function(data_dir){
  fns <- list.files(data_dir)
  fns <- fns[grepl('logfc_results.csv', fns)]
  pair_groups <- data.frame(result_fn=fns, stringsAsFactors=F)
  pair_groups$datatag <- get_info(pair_groups$result_fn, idx=2)
  pair_groups$title <- pair_groups$series
  pair_groups$clone1 <- get_info(pair_groups$result_fn, idx=6)
  pair_groups$clone2 <- get_info(pair_groups$result_fn, idx=8)
  # View(pair_groups)
  
  pair_groups$desc <- gsub('_logfc_results.csv$','', pair_groups$result_fn)
  pair_groups$desc <- gsub('^scrande_','', pair_groups$desc)  
  pair_groups$comp_type <- ifelse(!grepl('T',pair_groups$desc),'untreated','treated_vs_untreated')
  data.table::fwrite(pair_groups, paste0(data_dir,'comparisons_drug_res.csv'), quote=F)
  return(pair_groups)
}
# pair_groups_backup <- pair_groups
# dim(pair_groups)
# pair_groups <- pair_groups_backup
calculate_cis_trans_proportions <- function(input_dir, FDR_cutoff=0.01, 
                                            minLogFC=0.5, pValueThrs=0.05){
  # datatag <- 'SA609'
  # get_gene_type(datatag, minLogFC, pair_groups_fn=NULL,
  #               input_dir=NULL, dlp_dir=NULL, save_dir=NULL, FDR_cutoff, pValueThrs)
  # datatag <- 'SA1035'
  # get_gene_type(datatag, minLogFC, pair_groups_fn=NULL,
  #               input_dir=NULL, dlp_dir=NULL, save_dir=NULL, FDR_cutoff, pValueThrs)
  # 
  # datatag <- 'SA535'
  # get_gene_type(datatag, minLogFC, pair_groups_fn=NULL,
  #               input_dir=NULL, dlp_dir=NULL, save_dir=NULL, FDR_cutoff, pValueThrs)
  # 
  # datatag <- 'SA530'
  # get_gene_type(datatag, minLogFC, pair_groups_fn=NULL,
  #               input_dir=NULL, dlp_dir=NULL, save_dir=NULL, FDR_cutoff, pValueThrs)
  # datatag <- 'SA501'
  # get_gene_type(datatag, minLogFC, pair_groups_fn=NULL,
  #               input_dir=NULL, dlp_dir=NULL, save_dir=NULL, FDR_cutoff, pValueThrs)
  # 
  # datatag <- 'SA604'
  # get_gene_type(datatag, minLogFC, pair_groups_fn=NULL,
  #               input_dir=NULL, dlp_dir=NULL, save_dir=NULL, FDR_cutoff, pValueThrs)
        
  datatags <- c('SA501','SA530','SA604','SA609','SA535','SA1035')
  stat <- tibble::tibble()
  for(datatag in datatags){
    stat_out <- get_gene_type(datatag, minLogFC, pair_groups_fn=NULL,
                  input_dir=NULL, dlp_dir=NULL, save_dir=NULL, FDR_cutoff, pValueThrs)
    # tmp <- data.table::fread(paste0(input_dir,datatag,'_cistrans_summary.csv'))
    stat <- dplyr::bind_rows(stat, stat_out)
  }
  stat <- as.data.frame(stat)
  dim(stat)
  # View(stat)
  data.table::fwrite(stat, paste0(input_dir,'wholedata_summary.csv'), quote=F)
  stat$PDX <- stat$datatag
  stat_backup <-stat
  stat1 <- stat
  stat1$plt_desc <- stat1$desc
  
  plot_intrans_incis_prevalence_v2()
    
  
}

get_gene_type <- function(datatag, minLogFC=0.5, pair_groups_fn=NULL,
                          input_dir=NULL, dlp_dir=NULL, save_dir=NULL, FDR_cutoff=0.01, pValueThrs=0.05){
  if(is.null(dlp_dir)){
    dlp_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/clonealign/')  
  }
  if(is.null(save_dir)){
    save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/cis_trans_scran/')
  }  
  if(is.null(input_dir)){
    input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
  }
  if(is.null(pair_groups_fn)){
    pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
  }  
  
  # rm(cnv_mat)
  cnv_mat <- data.table::fread(paste0(dlp_dir,datatag,'_cnv_mat.csv'))%>% as.data.frame()
  print(head(cnv_mat))
  rownames(cnv_mat) <- cnv_mat$V1
  cnv_mat$V1 <- NULL
  print(dim(cnv_mat))
  # cnv_mat <- cnv_mat %>% 
  #   tibble::column_to_rownames(var="ensembl_gene_id")
  pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  pair_groups <- pair_groups[pair_groups$datatag==datatag,]
  stat_out <- get_gene_type_edgeR_v2(pair_groups, cnv_mat, 
                         datatag, input_dir, save_dir,
                         cancer_ref_genes_fn=NULL, 
                         FDR_cutoff, minLogFC, pValueThrs)
  return(stat_out)
  
}

# input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
viz_dotplot_proportions <- function(input_dir){
  datatags <- c('SA501','SA530','SA604','SA609','SA535','SA1035')
  stat <- tibble::tibble()
  for(datatag in datatags){
    tmp <- data.table::fread(paste0(input_dir,datatag,'_cistrans_summary.csv'))
    stat <- dplyr::bind_rows(stat, tmp)
  }
  stat <- as.data.frame(stat)
  dim(stat)
  # View(stat)
  data.table::fwrite(stat, paste0(input_dir,'wholedata_summary.csv'), quote=F)
}


# abs(logFC) > 0.5
# 2555 SA609_1_SA609_UT_R_UU_H
# 1424 SA609_10_SA609_UUUU_H_UUUU_C
# 3938 SA609_2_SA609_UTT_R_UUU_H
# 2697 SA609_3_SA609_UTTT_R_UUUU_H
# 2404  SA609_4_SA609_UTTTT_R_UUUUU_H
# 
# abs(logFC) > 0.25
# 4413 SA609_1_SA609_UT_R_UU_H
# 4311 SA609_10_SA609_UUUU_H_UUUU_C
# 7024 SA609_2_SA609_UTT_R_UUU_H
# 5557 SA609_3_SA609_UTTT_R_UUUU_H
# 5148  SA609_4_SA609_UTTTT_R_UUUUU_H


