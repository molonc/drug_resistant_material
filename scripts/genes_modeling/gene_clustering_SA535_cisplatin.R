  script_dir <- '/home/htran/Projects/farhia_project/rscript/genes_modeling/'
  source(paste0(script_dir, "gene_utils.R"))
  
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  # gene_type <- 'increase'
  gene_type <- 'decrease'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/gene_regression_',gene_type,'_cis/')
  save_dir
  
  datatag <- 'SA535'
  obs_clones <- c('R','S','T')  #
  # clone_aware <- FALSE
  clone_aware <- TRUE
  obs_treatment_st <- c('UUT','UUTT','UUTTT','UUTTTT') #,'UUTTTTT'
  obs_clones_untreated <- c('J')
  obs_untreated_st <- c('UU','UUU','UUUU','UUUUU')
  
  dir.create(save_dir, recursive = T, showWarnings = F)
  # Define clone R: R, clone S: H
  # Read csv files of DE genes
  de_x4 = read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/pathway4Tvs4U_SA535_UUTTTT_S_T_UUUUUU_J_Q_logfc_results.csv')) #Q here no make sense
  # de_x4 = read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UUTTTT_S_T_UUUUUU_J_Q/signif_genes.csv'))
  de_x3 = read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UUTTT_S_T_UUUUU_J/signif_genes.csv'))
  # de_x4 = read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UUTTT_T_UUUUU_J/signif_genes.csv'))
  dim(de_x3)
  dim(de_x4)
  
  FDR_cutoff <- 0.01
  minLogFC <- 0.25
  pValueThrs <- 0.05
  de_x4 <- de_x4 %>%
    dplyr::filter(abs(logFC)>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
  
  
  # de_x3 <- get_top_up_genes(de_x3, ntop = 500)
  # de_x4 <- get_top_up_genes(de_x4, ntop = 500)
  # print(dim(de_x3))
  # sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA535-v6/sce_cis_clones.rds'))
  sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA535-v6/SA535_cisplatin_sctransform_normalized.rds')) #total_sce_clones.rds
  # rownames(sce) <- rowData(sce)$ensgene
  print(dim(sce))
  rownames(sce)[1]
  print(summary(as.factor(sce$clone)))
  if(is.null(rowData(sce)$Symbol)){
    genes_symb_df <- read.csv(paste0(base_dir,'biodatabase/meta_genes.csv'), check.names = F, stringsAsFactors = F)
    dim(genes_symb_df)
    rownames(genes_symb_df) <- genes_symb_df$gene_ens
    rowData(sce)$Symbol <- genes_symb_df[rownames(sce),'gene_symb']
  }
  
  # Select only resistant genes 
  # For testing, get top genes first
  if(gene_type=='increase'){
    print('option increase')
    observed_genes <- intersect(de_x3[de_x3$logFC>0,'ensembl_gene_id'], de_x4[de_x4$logFC>0,'ensembl_gene_id'])
  }else{
    print('option decrease')
    observed_genes <- intersect(de_x3[de_x3$logFC<0,'ensembl_gene_id'], de_x4[de_x4$logFC<0,'ensembl_gene_id'])
  }
  
  # sce$clone <- ifelse(sce$clone=='TRUE','T',sce$clone)
  # sum(sce$clone=='TRUE')
  
  # Select only resistant genes 
  # For testing, get top genes first
  # observed_genes <- de_x3$ensembl_gene_id
  observed_genes <- intersect(observed_genes, rownames(sce))
  length(observed_genes)
  print(paste0("Number of observed genes: ", length(observed_genes)))
  # de_x3 <- de_x3[de_x3$ensembl_gene_id %in% observed_genes,]
  # de_x3 <- de_x3 %>%
  #   dplyr::filter(ensembl_gene_id %in% observed_genes)
  # print(dim(de_x3))
  # # de_x3 <- get_top_up_genes(de_x3, ntop = 500)
  # observed_genes <- as.character(de_x3$ensembl_gene_id)
  
  

  meta_cells_df <- as.data.frame(colData(sce))
  # meta_cells_df$cell_id <- rownames(meta_cells_df)
  # meta_cells_df$treatmentSt <- get_treatment_status(meta_cells_df$series)
  print(unique(meta_cells_df$treatmentSt))
  print(summary(as.factor(meta_cells_df$clone)))
  meta_cells_df <- meta_cells_df %>%
    dplyr::select(cell_id, treatmentSt, timepoint, clone) %>%
    dplyr::mutate(treatment_status=paste0(timepoint,'_',treatmentSt))
  
  if(clone_aware){
    print("With clone aware")
    meta_cells_untreated <- meta_cells_df %>%
      dplyr::filter(clone %in% obs_clones_untreated & treatmentSt %in% obs_untreated_st)#  
    
    
    meta_cells_df <- meta_cells_df %>%
      dplyr::filter(clone %in% obs_clones & treatmentSt %in% obs_treatment_st)
    
  }else{
    print("Without clone aware")
    meta_cells_untreated <- meta_cells_df %>%
      dplyr::filter(treatmentSt %in% obs_untreated_st)
   
    meta_cells_df <- meta_cells_df %>%
      dplyr::filter(treatmentSt %in% obs_treatment_st)
  }
  print(summary(meta_cells_untreated$treatmentSt))
  print(summary(meta_cells_df$treatmentSt))
  print(dim(meta_cells_df))
  print(dim(meta_cells_untreated))
  # meta_cells_df <- downsample_cells(meta_cells_df, downsample_ratio=0.3,
  #                              thres_small_clone=500)
  dim(meta_cells_df)
  
  
  table(meta_cells_df$clone, meta_cells_df$treatmentSt)
  table(meta_cells_untreated$clone, meta_cells_untreated$treatmentSt)
  
  observed_sce_untreated <- sce[observed_genes, meta_cells_untreated$cell_id]
  observed_sce <- sce[observed_genes, meta_cells_df$cell_id]
  dim(observed_sce)
  dim(observed_sce_untreated)
  # observed_sce <- observed_sce[,meta_cells_df$cell_id]
  # View(rowData(observed_sce))
  # norm_data <- logcounts(observed_sce)
  # norm_data <- as.data.frame(as.matrix(norm_data))
  # norm_data$gene <- rowData(observed_sce)$Symbol
  # dim(norm_data)
  
  
  meta_genes <- data.frame(gene=rowData(observed_sce)$Symbol, ens_gene=rownames(observed_sce), stringsAsFactors=F)
  de_x3 <- de_x3 %>%
    dplyr::filter(ensembl_gene_id %in% meta_genes$ens_gene) %>%
    dplyr::select(logFC, ensembl_gene_id, Gene_Type) %>%
    dplyr::rename(logFC_UUTTT_S_T_UUUUU_J=logFC)
  
  de_x4 <- de_x4 %>%
    dplyr::filter(ensembl_gene_id %in% meta_genes$ens_gene) %>%
    dplyr::select(logFC, ensembl_gene_id) %>%
    dplyr::rename(logFC_UUTTTT_S_T_UUUUUU_J=logFC)
  meta_genes <- meta_genes %>% left_join(de_x3, by=c('ens_gene'='ensembl_gene_id'))
  meta_genes <- meta_genes %>% left_join(de_x4, by=c('ens_gene'='ensembl_gene_id'))
  dim(meta_genes)
  summary(as.factor(meta_genes$Gene_Type))
  # gene_clustering(meta_genes, norm_data, meta_cells_df, obs_clones, 
  #                 save_dir, datatag, TRUE)
  # 
  
  # Input data to WGCNA
  dim(observed_sce)
  norm_data <- logcounts(observed_sce)
  datExpr <- as.data.frame(as.matrix(t(norm_data))) # # now cells are rows and genes are columns
  tag <- 'Cisplatin'
  res <- remove_outliers(datExpr, save_dir, datatag, tag)
  datExpr <- res$datExpr
  print(dim(datExpr))
  norm_data <- as.data.frame(as.matrix(norm_data))
  norm_data[norm_data==0] <- NA
  norm_data$gene <- rowData(observed_sce)$Symbol
  dim(norm_data)
  
  norm_data_untreated <- logcounts(observed_sce_untreated)
  norm_data_untreated <- as.data.frame(as.matrix(norm_data_untreated))
  norm_data_untreated[norm_data_untreated==0] <- NA
  norm_data_untreated$gene <- rowData(observed_sce_untreated)$Symbol
  dim(norm_data_untreated)
  
  
  
  metacells_df <- encode_metadata(meta_cells_df, save_dir, datatag)
  table(rownames(metacells_df)==rownames(datExpr))
  
  select_soft_threshold_power(datExpr, save_dir, option=1)
  softPowerThres <- 12
  minClusterSz <- 40
  results_mtx <- compute_adjacency_correlation_mtx(datExpr, save_dir, softPowerThres, 1, T)
  
  results_mtx$diss_mtx[1:3,1:3]
  results <- generate_cluster_gene_tree_v1(meta_genes, datExpr, results_mtx$diss_mtx, save_dir, softPowerThres, minClusterSz, T)  
  # results <- generate_cluster_gene_tree_v2(datExpr, results_mtx$diss_mtx, save_dir, minClusterSz)  
  # correlate_genes_with_treatment_effect(results, datExpr, metacells_df, save_dir)
  # gene_cluster <- results$gene_cluster
  # gene_cluster[gene_cluster$cluster %in% c('grey'),'cluster'] <- 'blue'
  # summary(as.factor(gene_cluster$cluster))
  # write.csv(gene_cluster, paste0(save_dir,'genes_clusters.csv'), quote=F, row.names = F)
  # results$gene_cluster <- gene_cluster
  # save(results, file= paste0(save_dir, "Network_v1.RData"))
  
  
  if(gene_type=='increase'){
    plttitle <- paste0(datatag,': up-regulated DEG in Rx vs UnRx')
  }else{
    plttitle <- paste0(datatag,': down-regulated DEG in Rx vs UnRx')
  }
  plttitle
  plot_genes_exp(results$gene_cluster, meta_genes, norm_data, norm_data_untreated,meta_cells_untreated,
                 meta_cells_df, obs_clones, save_dir, datatag, clone_aware=TRUE, plttitle)
  


  
  
  plot_actual_genes_exp(results$gene_cluster, meta_genes, norm_data, norm_data_untreated,meta_cells_untreated,
                        meta_cells_df, obs_clones, save_dir, datatag, clone_aware)
  
  
  clusters_ls <- unique(results$gene_cluster$cluster)
  for(cluster_use in clusters_ls){
    net <- extract_genes_correlation_network(results_mtx, results, datExpr, save_dir, datatag, softPowerThres, cluster_use)
    construct_graph(net$edges_df, net$node_df, save_dir, datatag, cluster_use)
  }
  
  
  load(paste0(save_dir, "Network_v1.RData"))  # load results objects
  tag <- 'Cisplatin'
  datExpr <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_datExpr.csv.gz'))
  datExpr <- as.data.frame(datExpr)
  rownames(datExpr) <- datExpr$V1
  datExpr$V1 <- NULL
  metacells_df <- read.csv(paste0(save_dir, datatag,'_meta_cells.csv'), check.names = F, stringsAsFactors = F)
  
  rownames(metacells_df) <- metacells_df$cell_id
  metacells_df$cell_id <- NULL
  correlate_genes_with_treatment_effect(results, datExpr, metacells_df, save_dir)
  
  
  
  
  
  # base_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA501_v2/'
  # mat <- data.table::fread(paste0(base_dir,'total_merged_filtered_states.csv'))  
  # View(mat[1:5,1:5])
  