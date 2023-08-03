

source("/home/htran/Projects/farhia_project/rnaseq/cis_trans/in_cis_trans_utils.R")
source("/home/htran/Projects/farhia_project/rnaseq/cis_trans/viz_utils.R")


# 
# pts <- c('Pt4','Pt5','Pt6')
# names(pts) <- series


series=c("SA501","SA530","SA604","SA609","SA535","SA1035")
pts=paste0("Pt",rep(1:6,1))
names(pts) <- series

FDR_cutoff=0.01
minLogFC=0.5
pValueThrs=0.05
# data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
data_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/cis_trans/'
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v5.csv')
save_dir <- paste0(data_dir,'signif_genes/')
# datatag <- 'SA609'
save_csv_dir <- paste0(data_dir,'signif_genes/csv_submission/')
# dir.create(save_csv_dir)
comp_treated <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_treated.csv')) %>% as.data.frame()
comp_treated <- comp_treated %>%
  dplyr::filter(file_header!='scrande_SA1035_41') %>% # early time and do not contain cis genes. 
  dplyr::select(file_header, order)

save_to_csv_treated <- function(){
  
  # View(head(comp_treated))
  # pair_groups$desc <- paste0(pair_groups$patient,gsub(' ','_',pair_groups$labels_detailed))
  treated_series <- c('SA609','SA535','SA1035')
  total <- tibble::tibble()
  for(datatag in treated_series){
    pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
    dim(pair_groups)
    # View(head(pair_groups))
    unique(pair_groups$datatag)
    pair_groups <- pair_groups[pair_groups$datatag==datatag,] 
    rownames(pair_groups) <- pair_groups$file_header
    dim(pair_groups)
    pair_groups$patient <- pts[pair_groups$datatag]
    pair_groups$desc <- paste0(pair_groups$patient,' ',pair_groups$labels_detailed)
    # View(pair_groups)
    pair_groups <- pair_groups %>%
      dplyr::filter(file_header %in% comp_treated$file_header)
    
    
    for(de in pair_groups$file_header){
      fn <- paste0(save_dir,pair_groups[de,'file_header'],"/signif_genes.csv")
      if(file.exists(fn)){
        df <- data.table::fread(fn)
        df <- df %>%
          dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
        # df$patient <- pts[datatag]
        print(dim(df))
        df <- df %>%
          dplyr::select(ensembl_gene_id, gene_symbol, logFC, PValue, FDR, Gene_Type)
        df$DE_comparison <- pair_groups[de,'desc']
        df$logFC <- round(df$logFC, 2)
        total <- dplyr::bind_rows(total, df)
      }
    }
    
    
    
    
  }
  lb <- paste(pts[treated_series], collapse = '_')
  data.table::fwrite(total, paste0(save_csv_dir, 'DE_cis_trans_genes_',lb,'.csv.gz'))
  dim(total)
  
}


save_to_csv_untreated <- function(){
  
  # View(head(comp_treated))
  # pair_groups$desc <- paste0(pair_groups$patient,gsub(' ','_',pair_groups$labels_detailed))
  untreated_series <- c("SA501","SA530","SA604")
  pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  dim(pair_groups)
  # View(head(pair_groups))
  unique(pair_groups$datatag)
  # pair_groups <- pair_groups[pair_groups$datatag==datatag,] 
  rownames(pair_groups) <- pair_groups$file_header
  dim(pair_groups)
  pair_groups$patient <- pts[pair_groups$datatag]
  pair_groups$desc <- paste0(pair_groups$patient,' ',pair_groups$labels_detailed)
  # View(pair_groups)
  comp_untreated_DH <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_untreated_DH.csv')) %>% as.data.frame()
  pair_groups <- pair_groups %>%
    dplyr::filter(file_header %in% comp_untreated_DH$file_header)
  dim(pair_groups)
  
  total <- tibble::tibble()
  for(datatag in untreated_series){
    for(de in pair_groups$file_header){
      fn <- paste0(save_dir,pair_groups[de,'file_header'],"/signif_genes.csv")
      if(file.exists(fn)){
        df <- data.table::fread(fn)
        df <- df %>%
          dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
        # df$patient <- pts[datatag]
        print(dim(df))
        df <- df %>%
          dplyr::select(ensembl_gene_id, gene_symbol, logFC, PValue, FDR, Gene_Type)
        df$DE_comparison <- pair_groups[de,'desc']
        df$logFC <- round(df$logFC, 2)
        total <- dplyr::bind_rows(total, df)
      }
    }
    
    
    
    
  }
  lb <- paste(pts[untreated_series], collapse = '_')
  data.table::fwrite(total, paste0(save_csv_dir, 'DE_cis_trans_genes_untreated_drugHoliday_',lb,'.csv.gz'))
  dim(total)
  
}

# View(head(total))

save_to_csv_untreated <- function(){
  ## For untreated clones and drug holiday 
  pair_groups <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  rownames(pair_groups) <- pair_groups$file_header
  dim(pair_groups)
  # View(pair_groups)
  
  comp_untreated_DH <- data.table::fread(paste0(data_dir, 'comparisons_drug_res_v5_untreated_DH.csv')) %>% as.data.frame()
  # View(comp_untreated_DH)
  pair_groups <- pair_groups %>%
    dplyr::filter(file_header %in% comp_untreated$file_header)
  dim(pair_groups)
  pair_groups <- pair_groups[pair_groups$datatag==datatag,] 
  rownames(pair_groups) <- pair_groups$file_header
  dim(pair_groups)
  pair_groups$desc <- paste0(pair_groups$patient,' ',pair_groups$labels_detailed)
  
  treated_series <- c('SA609','SA535','SA1035')
  total <- tibble::tibble()
  for(datatag in treated_series){
  for(de in pair_groups$file_header){
    fn <- paste0(save_dir,pair_groups[de,'file_header'],"/signif_genes.csv")
    if(file.exists(fn)){
      df <- data.table::fread(fn)
      df <- df %>%
        dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
      df$patient <- pair_groups[de,'datatag']
      print(dim(df))
      df <- df %>%
        dplyr::select(ensembl_gene_id, gene_symbol, logFC, PValue, FDR, Gene_Type, patient)
      df$file_header <- de
      total <- dplyr::bind_rows(total, df)  
    }
  }  
  }  
  dim(total)  
  pair_groups <- pair_groups %>%
    dplyr::select(file_header, clone1, clone2, series)
  total <- total %>%
    dplyr::left_join(pair_groups, by='file_header') %>%
    dplyr::select(-file_header)
  print(datatag)
  print(dim(total))
  data.table::fwrite(total, paste0(save_csv_dir, 'DE_cis_trans_genes_Pt1-6_untreated_DH_comparisons.csv.gz'))
  # View(head(total_609))
  data.table::fwrite(total, paste0(save_csv_dir, 'DE_cis_trans_genes_', pts[datatag],'.csv.gz'))
  
}
