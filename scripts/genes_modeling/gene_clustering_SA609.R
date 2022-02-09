    script_dir <- '/home/htran/Projects/farhia_project/rnaseq/genes_modeling/'
    source(paste0(script_dir, "gene_utils.R"))
    
    base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
    # save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/gene_regression_decrease/'
    gene_type <- 'increase'
    # gene_type <- 'decrease'
    save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/gene_regression_',gene_type,'_v2/')
    save_dir
    datatag <- 'SA609'
    obs_clones <- c('R')
    # clone_aware <- FALSE
    clone_aware <- TRUE
    obs_treatment_st <- c('UT','UTT','UTTT','UTTTT','UTTTTT')
    obs_clones_untreated <- c('H')
    obs_untreated_st <- c('UU','UUU','UUUU','UUUUU')
    # input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
    dir.create(save_dir, recursive = T, showWarning=F)
    # Define clone R: R, clone S: H
    # Read csv files of DE genes
    de_x3 = read.csv(paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTT_R_UUUU_H/signif_genes.csv'))
    de_x4 = read.csv(paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTTT_R_UUUUU_H/signif_genes.csv'))
    
    # summary(as.factor(de_x4$Gene_Type))
    
    # de_x3 <- get_top_up_genes(de_x3, ntop = 500)
    # de_x4 <- get_top_up_genes(de_x4, ntop = 500)
    # print(dim(de_x3))
    sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_normalized.rds'))
    # sce <- corrected_sce
    # sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/total_sce_clones.rds'))
    print(dim(sce))
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
    
    # 
    length(observed_genes)
    print(summary(as.factor(de_x4[de_x4$ensembl_gene_id %in% observed_genes,'Gene_Type'])))
    observed_genes <- intersect(observed_genes, rownames(sce))
    length(observed_genes)
    # de_x3 <- de_x3 %>%
    #          dplyr::filter(ensembl_gene_id %in% observed_genes)
    # print(dim(de_x3))
    # # de_x3 <- get_top_up_genes(de_x3, ntop = 500)
    # observed_genes <- as.character(de_x3$ensembl_gene_id)
    print(paste0("Number of observed genes: ", length(observed_genes)))
    
    meta_cells_df <- as.data.frame(colData(sce))
    # meta_cells_df$cell_id <- rownames(meta_cells_df)
    meta_cells_df <- meta_cells_df %>%
      dplyr::select(cell_id, treatmentSt, timepoint, clone) %>%
      dplyr::mutate(treatment_status=paste0(timepoint,'_',treatmentSt))
    if(clone_aware){
      print("With clone aware")
      meta_cells_untreated <- meta_cells_df %>%
        dplyr::filter(clone %in% obs_clones_untreated & treatmentSt %in% obs_untreated_st)
      print(dim(meta_cells_untreated))
      
      meta_cells_df <- meta_cells_df %>%
        dplyr::filter(clone %in% obs_clones & treatmentSt %in% obs_treatment_st)
      print(dim(meta_cells_df))
    }else{
      print("Without clone aware")
      meta_cells_untreated <- meta_cells_df %>%
        dplyr::filter(treatmentSt %in% obs_untreated_st)
      print(dim(meta_cells_untreated))
      
      meta_cells_df <- meta_cells_df %>%
        dplyr::filter(treatmentSt %in% obs_treatment_st)
      print(dim(meta_cells_df))
    }
    print(unique(meta_cells_untreated$treatmentSt))
    print(unique(meta_cells_df$treatmentSt))
    dim(meta_cells_df)
    # meta_cells_df <- downsample_cells(meta_cells_df, downsample_ratio=0.3,
    #                              thres_small_clone=500)
    # colnames(meta_cells_df)
    # print(dim(meta_cells_df))
    # meta_cells_df$timepoint[1:3]
    
    observed_sce_untreated <- sce[observed_genes, meta_cells_untreated$cell_id]
    observed_sce <- sce[observed_genes, meta_cells_df$cell_id]
    dim(observed_sce_untreated)
    dim(observed_sce)
    
    # View(rowData(observed_sce))
    # norm_data <- logcounts(observed_sce)
    # norm_data <- as.data.frame(as.matrix(norm_data))
    # norm_data$gene <- rowData(observed_sce)$Symbol
    # dim(norm_data)
    
    meta_genes <- data.frame(gene=rowData(observed_sce)$Symbol,ens_gene=rownames(observed_sce), stringsAsFactors=F)
    # meta_genes <- meta_genes %>% left_join(genes_symb_df, by=c('ens_gene'='gene_ens'))
    
    # de_x3 <- de_x3 %>%
    #   dplyr::filter(ensembl_gene_id %in% meta_genes$ens_gene) %>%
    #   dplyr::select(logFC, ensembl_gene_id) %>%
    #   dplyr::rename(logFC_UTTT_R_UUUU_H=logFC)
    # 
    # de_x4 <- de_x4 %>%
    #   dplyr::filter(ensembl_gene_id %in% meta_genes$ens_gene) %>%
    #   dplyr::select(logFC, ensembl_gene_id, Gene_Type) %>%
    #   dplyr::rename(logFC_UTTTT_R_UUUUU_H=logFC)
    # meta_genes <- meta_genes %>% left_join(de_x3, by=c('ens_gene'='ensembl_gene_id'))
    # meta_genes <- meta_genes %>% left_join(de_x4, by=c('ens_gene'='ensembl_gene_id'))
    # head(meta_genes)
    
    # gene_clustering(meta_genes, norm_data, meta_cells_df, obs_clones, 
    #                 save_dir, datatag, TRUE)
    
    # Input data to WGCNA
    norm_data <- logcounts(observed_sce)
    dim(observed_sce)
    # norm_data <- logcounts(sce[observed_genes,])  #include all untreated and treated cells to comparison
    dim(norm_data)
    datExpr <- as.data.frame(as.matrix(t(norm_data))) # # now cells are rows and genes are columns
    tag <- ''
    res <- remove_outliers(datExpr, save_dir, datatag, tag)
    datExpr <- res$datExpr
    print(dim(datExpr))
    norm_data <- as.data.frame(as.matrix(norm_data))
    norm_data[norm_data==0] <- NA
    norm_data$gene <- rowData(observed_sce)$Symbol
    # norm_data$gene <- rownames(observed_sce)
    dim(norm_data)
    
    norm_data_untreated <- logcounts(observed_sce_untreated)
    norm_data_untreated <- as.data.frame(as.matrix(norm_data_untreated))
    norm_data_untreated[norm_data_untreated==0] <- NA
    norm_data_untreated$gene <- rowData(observed_sce_untreated)$Symbol
    dim(norm_data_untreated)
    # sum(is.na(norm_data_untreated))
    
    
    metacells_df <- encode_metadata(meta_cells_df, save_dir, datatag)
    table(rownames(metacells_df)==rownames(datExpr))
    
    select_soft_threshold_power(datExpr, save_dir, option=1)
    softPowerThres <- 10
    minClusterSz <- 10 #10
    results_mtx <- compute_adjacency_correlation_mtx(datExpr, save_dir, softPowerThres, 1, T)
    # rm(results)
    results <- generate_cluster_gene_tree_v1(meta_genes, datExpr, results_mtx$diss_mtx, save_dir, softPowerThres, minClusterSz, T)
    # results <- generate_cluster_gene_tree_v2(datExpr, results_mtx$diss_mtx, save_dir, minClusterSz)
    # View(head(results$dynamicColors))
    # summary(as.factor(results$dynamicColors))
    # gene_cluster$cluster <- results$dynamicColors
    
    ##Reload results
    load(paste0(save_dir, "Network_v1.RData"))  # load results objects
    tag <- ''
    # datExpr <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_datExpr.csv.gz'))
    # datExpr <- as.data.frame(datExpr)
    # rownames(datExpr) <- datExpr$V1
    # datExpr$V1 <- NULL
    # metacells_df <- read.csv(paste0(save_dir, datatag,'_meta_cells.csv'), check.names = F, stringsAsFactors = F)
    # rownames(metacells_df) <- metacells_df$cell_id
    # metacells_df$cell_id <- NULL
    
    
    # correlate_genes_with_treatment_effect(results, datExpr, metacells_df, save_dir)
    gene_cluster <- results$gene_cluster
    dim(gene_cluster)
    unique(gene_cluster$cluster)
    View(head(gene_cluster))
    gene_cluster[gene_cluster$cluster %in% c('brown'),'cluster'] <- 'grey'
    # gene_cluster[gene_cluster$cluster %in% c('brown'),'cluster'] <- 'yellow'
    # gene_cluster[gene_cluster$cluster %in% c('turquoise'),'cluster'] <- 'green'
    # # gene_cluster[gene_cluster$cluster_id==0,'cluster_id'] <- 2
    # summary(as.factor(gene_cluster$cluster))
    # write.csv(gene_cluster, paste0(save_dir,'genes_clusters.csv'), quote=F, row.names = F)
    # results$gene_cluster <- gene_cluster
    # save(results, file= paste0(save_dir, "Network_v1.RData"))
    if(gene_type=='increase'){
      plttitle <- paste0(datatag,': up-regulated DEG in Rx vs UnRx')
    }else{
      plttitle <- paste0(datatag,': down-regulated DEG in Rx vs UnRx')
    }
    # plttitle
    plot_genes_exp(gene_cluster, meta_genes, norm_data, norm_data_untreated,meta_cells_untreated,
                   meta_cells_df, obs_clones, save_dir, datatag, clone_aware=TRUE, plttitle)
    plot_actual_genes_exp(gene_cluster, meta_genes, norm_data, norm_data_untreated,meta_cells_untreated,
                   meta_cells_df, obs_clones, save_dir, datatag, clone_aware=TRUE)

    
    
    clusters_ls <- unique(results$gene_cluster$cluster)
    for(cluster_use in clusters_ls){
      net <- extract_genes_correlation_network(results_mtx, results, datExpr, save_dir, datatag, softPowerThres, cluster_use)
      construct_graph(net$edges_df, net$node_df, save_dir, datatag, cluster_use)
    }
    
    
    
    # linearMod <- lm(mean_exp_scaled ~ treatment_status, data=gexp_cls1)
    # gexp_cls1$slope <- linearMod$coefficients[2]  
    # mon <- (gexp_cls1[gexp_cls1$treatment_status=="X4_UT","mean_exp_scaled"]-gexp_cls1[gexp_cls1$treatment_status=="X5_UTT","mean_exp_scaled"]) * 
    #   (gexp_cls1[gexp_cls1$treatment_status=="X6_UTTT","mean_exp_scaled"]-gexp_cls1[gexp_cls1$treatment_status=="X5_UTT","mean_exp_scaled"])
    # if (mon > 0) {
    #   m <- TRUE
    # } else {
    #   m <- FALSE
    # }
    # genedf[genedf$name==gene,"monotonic"] <- m
    # 
    # gexp_cls1$treatment_status
    # p <- ggplot(gexp_cls1, aes(as.factor(treatment_status), mean_exp_scaled)) +
    #   geom_point()
    # p + geom_smooth(method = lm)
    
    
    
    