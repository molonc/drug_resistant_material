# pkgs <- c("factoextra",  "NbClust","Rfast") 
# install.packages(pkgs)
# BiocManager::install("Rfast")
# library(factoextra)
# library(NbClust)
# library(Rfast)
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(SingleCellExperiment)
})

gene_clustering <- function(meta_genes, norm_data, meta_cells_df, obs_clones, 
                            save_dir, datatag, clone_aware=TRUE){
  tag <- paste(obs_clones, collapse = '_') # tag <- 'H_A_K_I_J' SA1035
  
  exp_mean_df <- norm_data %>% 
    # convert to long format
    pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
    # join with sample info table
    full_join(meta_cells_df, by = ("cell_id")) %>% 
    # filter to retain only genes of interest
    # filter(gene %in% observed_genes) %>% 
    # for each gene
    group_by(gene) %>% 
    # scale the cts column
    mutate(exp_scaled = (exp - mean(exp, na.rm=T))/sd(exp)) 
  
  print(summary(exp_mean_df$exp_scaled))
 
  print(summary(exp_mean_df$exp_scaled))
  data.table::fwrite(x = exp_mean_df, paste0(save_dir, datatag,'_',tag,'_exp_mean_df.csv.gz'))
  data.table::fwrite(x = norm_data, paste0(save_dir, datatag,'_',tag,'_norm_data.csv.gz'))
  data.table::fwrite(x = meta_cells_df, paste0(save_dir, datatag,'_',tag,'_meta_cells.csv.gz'))
  # 
  gexp_cls <- exp_mean_df %>% 
    # for each gene, summary by treatment_status
    group_by(gene, treatment_status) %>%
    # calculate the mean (scaled) cts
    summarise(mean_exp_scaled = mean(exp_scaled, na.rm=T),
              nrep = n()) %>% 
    ungroup()
  
  
  dim(gexp_cls)
  length(unique(gexp_cls$gene))
  print(summary(gexp_cls$mean_exp_scaled))
  
  # Create a matrix
  print("Compute distance matrix:")
  hclust_mtx <- norm_data %>%
    select(-gene) %>%
    as.matrix()
  
  # assign rownames
  rownames(hclust_mtx) <- norm_data$gene
  dim(hclust_mtx)
  # hclust_mtx <- hclust_mtx[observed_genes, ]
  hclust_mtx <- hclust_mtx %>% 
    # transpose the matrix so genes are as columns
    t() %>% 
    # apply scalling to each column of the matrix (genes)
    scale() %>% 
    # transpose back so genes are as rows again
    t()
  
  
  print("Gene Regression:")
  gene_dist <- Rfast::Dist(hclust_mtx, method = "euclidean")
  data.table::fwrite(x = gene_dist, paste0(save_dir, datatag,'_',tag,'_gene_dist.csv.gz'))
  
  fitclusters <- NbClust::NbClust(hclust_mtx, diss=gene_dist, distance = NULL, 
                                  min.nc = 5, max.nc = 10, method = "kmeans", index = "silhouette") #index = "all"
  
  
  saveRDS(fitclusters, file=paste0(save_dir, datatag,'_',tag,'_fit_clusters.rds'))
  
  
  # gene_hclust <- hclust(gene_dist, method = "complete")
  # The default `plot()` function can be used to produce a simple dendrogram
  # plot(gene_hclust, labels = FALSE)
  # abline(h = 110, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
  # nbclusters <- 10
  # gene_cluster <- cutree(gene_hclust, k = nbclusters) %>% 
  #   # turn the named vector into a tibble
  #   enframe() 
  
  
  gene_cluster <- fitclusters$Best.partition %>% tibble::enframe()
  gene_cluster <- as.data.frame(gene_cluster)
  colnames(gene_cluster) <- c('gene', 'cluster')
  
  print(summary(as.factor(gene_cluster$cluster)))
  gene_cluster <- gene_cluster %>% left_join(meta_genes, by=c("gene"))
  data.table::fwrite(x = gene_cluster, paste0(save_dir, datatag,'_',tag,'_gene_cluster.csv'))
  # data.table::fwrite(x = gene_cluster, paste0(save_dir, datatag,'_',tag,'_gene_cluster_classified.csv'))
  
  gexp_cls <- gexp_cls %>% inner_join(gene_cluster, by = "gene")
  data.table::fwrite(x = gexp_cls, paste0(save_dir, datatag,'_',tag,'_gexp_cls.csv'))
  # saveRDS(gexp_cls, file=paste0(save_dir, datatag,'_',tag,'_gexp_cls.rds'))
  print("Plot gene clusters:")
  plot_gene_clusters(gexp_cls, tag, datatag, obs_clones, save_dir, clone_aware, 
                     xl="Time point - treatment status", yl="Mean z-score value", trim_data=F)
  
  print(dim(gexp_cls))
  print(summary(gexp_cls$mean_exp_scaled))
  print("Gene Regression:")
  # gene_regression_SA1035(exp_mean_df, gene_cluster, save_dir, paste0(datatag,'_',tag))
  # exp_mean_df <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_exp_mean_df.csv.gz'))
  # # View(head(exp_mean_df))
  # gene_cluster <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_gene_cluster.csv')) #
  # meta_cells_df <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_meta_cells.csv.gz'))
  # meta_cells_df <- as.data.frame(meta_cells_df)
  # gene_cluster <- as.data.frame(gene_cluster)
  # View(head(meta_cells_df))
  # length(intersect(unique(exp_mean_df$gene), gene_cluster$gene))
  
  # gene_cluster <- gene_regression_SA1035(exp_mean_df, gene_cluster, save_dir, paste0(datatag,'_',tag))
  gene_cluster <- gene_regression_SA609(exp_mean_df, gene_cluster, save_dir, paste0(datatag,'_',tag))
  table(gene_cluster$classify_gene, gene_cluster$cluster)
  # obs_clusters <- gene_cluster[gene_cluster$classify_gene %in% c('mono_incrc','mono_decrc'),'cluster']
  # print("Plot gene important clusters:")
  # if(!is.null(obs_clusters)){
  #   plot_gene_each_cluster(gexp_cls, datatag, obs_clones, as.character(obs_clusters), save_dir,
  #                          xl="Time point - treatment status", yl="Mean z-score value")
  # }
  
  
  
  # head(trans_cts_cluster)
  
}

plot_heatmap <- function(gene_cluster, norm_data, meta_cells_df){
  unique(gene_cluster$cluster)
  gene_cluster$cluster <- paste0('Cluster ',gene_cluster$cluster)
  rownames(gene_cluster) <- gene_cluster$ens_gene
  rownames(meta_cells_df) <- meta_cells_df$cell_id
  gene_cluster$ens_gene[1:3]
  norm_data <- norm_data %>%
    select(-gene) %>%
    as.matrix()
  row_cls <- gene_cluster[rownames(norm_data),'cluster']
  col_cls <- meta_cells_df[colnames(norm_data),'treatment_status']
  p <- ComplexHeatmap::Heatmap(norm_data, na_col = "white",
                               show_column_names=F,
                               show_row_names = F,
                               cluster_rows=F,cluster_columns=F,
                               name = "Genes Clusters", 
                               # row_order = sort(rownames(test)),
                               row_split= row_cls,
                               row_title_rot = 0,
                               # row_gap = unit(2, "mm"),
                               column_split = col_cls, 
                               column_title = "Resistant Genes Regression across time",
                               # column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 10),
                               row_names_gp = grid::gpar(fontsize = 10),
                               show_heatmap_legend = T,
                               # top_annotation=top_anno,
                               # left_annotation = left_anno,
                               # cell_fun = cell_func,
                               row_dend_reorder=F)
  # p
  
  
}

gene_regression_SA1035 <- function(exp_mean_df, gene_cluster, save_dir, datatag){
  # gene_cluster <- as.data.frame(gene_cluster)
  # exp_mean_df <- as.data.frame(exp_mean_df)
  epsilon <- 0.05 # error rate between different batches
  rownames(gene_cluster) <- gene_cluster$gene
  gene_cluster$classify_gene <- ''
  gene_cluster$slope <- NA
  for(cl in unique(gene_cluster$cluster)){
    print(paste0("Observe clone: ", cl))
    genes_use <- gene_cluster[gene_cluster$cluster==cl,]$gene
    tmp <- exp_mean_df[exp_mean_df$gene %in% genes_use,]
    # tmp <- exp_mean_df %>%
    #               dplyr::filter(gene %in% as.character(genes_use))
    tmp$treatment_status <- as.factor(tmp$treatment_status)
    dim(tmp)
    rg <- lm(exp_scaled ~ treatment_status, data = tmp)  
    saveRDS(rg, paste0(save_dir, datatag,'_cluster_',cl,'_lm.rds'))
    gene_cluster[genes_use,'slope'] <- rg$coefficients[2] 
    print(summary(rg))
    coeffs <- coef(rg)
    # Get means in each treatment
    ts_med <- tapply(tmp$exp_scaled, tmp$treatment_status, median)
    # print(ts_means)
    c2 <- coeffs[2] > 0
    c3 <- coeffs[3] > 0
    c21 <- ts_med[2] - ts_med[1] >= -epsilon
    c32 <- ts_med[3] - ts_med[2] >= -epsilon
    c31 <- ts_med[3] - ts_med[1] >= -epsilon
    # c43 <- ts_med[4] - ts_med[3] >= -epsilon
    
    c12 <- ts_med[1] - ts_med[2] >= -epsilon
    c23 <- ts_med[2] - ts_med[3] >= -epsilon
    c13 <- ts_med[1] - ts_med[3] >= -epsilon
    # c34 <- ts_med[3] - ts_med[4] >= -epsilon
    if(c21 & c31){
      if(c32){
        gene_type <- 'mono_incrc'
      } else{
        gene_type <- 'incrc'
      }
      print(gene_type)
    } else if(!c2 & !c3 & c12 & c13){
      if(c23){
        gene_type <- 'mono_decrc'
      }else{
        gene_type <- 'decrc'
      }
      print(gene_type)
      
    } else{
      gene_type <- 'undefined'
    }
    
    gene_cluster[genes_use,'classify_gene'] <- gene_type
    print(paste0("Cluster ",cl," and gene type is: ",gene_type))
    print(" ")
  }
  
  write.csv(gene_cluster, paste0(save_dir, datatag,'_gene_clusters_classified.csv'), quote = F, row.names = F)
  return(gene_cluster)
}

gene_regression_SA535 <- function(exp_mean_df, gene_cluster, save_dir, datatag){
  # gene_cluster <- as.data.frame(gene_cluster)
  # exp_mean_df <- as.data.frame(exp_mean_df)
  epsilon <- 0.05 # error rate between different batches
  rownames(gene_cluster) <- gene_cluster$gene
  gene_cluster$classify_gene <- ''
  gene_cluster$slope <- NA
  for(cl in unique(gene_cluster$cluster)){
    print(paste0("Observe clone: ", cl))
    genes_use <- gene_cluster[gene_cluster$cluster==cl,]$gene
    tmp <- exp_mean_df[exp_mean_df$gene %in% genes_use,]
    # tmp <- exp_mean_df %>%
    #               dplyr::filter(gene %in% as.character(genes_use))
    tmp$treatment_status <- as.factor(tmp$treatment_status)
    dim(tmp)
    rg <- lm(exp_scaled ~ treatment_status, data = tmp)  
    saveRDS(rg, paste0(save_dir, datatag,'_cluster_',cl,'_lm.rds'))
    gene_cluster[genes_use,'slope'] <- rg$coefficients[2] 
    print(summary(rg))
    coeffs <- coef(rg)
    # Get means in each treatment
    ts_med <- tapply(tmp$exp_scaled, tmp$treatment_status, median)
    # print(ts_means)
    c2 <- coeffs[2] > 0
    c3 <- coeffs[3] > 0
    c21 <- ts_med[2] - ts_med[1] >= -epsilon
    c32 <- ts_med[3] - ts_med[2] >= -epsilon
    c31 <- ts_med[3] - ts_med[1] >= -epsilon
    # c43 <- ts_med[4] - ts_med[3] >= -epsilon
    
    c12 <- ts_med[1] - ts_med[2] >= -epsilon
    c23 <- ts_med[2] - ts_med[3] >= -epsilon
    c13 <- ts_med[1] - ts_med[3] >= -epsilon
    # c34 <- ts_med[3] - ts_med[4] >= -epsilon
    if(c21 & c31){
      if(c2 & c3 & c32){
        gene_type <- 'mono_incrc'
      } else{
        gene_type <- 'incrc'
      }
      print(gene_type)
    } else if(!c2 & !c3 & c12 & c13){
      if(c23){
        gene_type <- 'mono_decrc'
      }else{
        gene_type <- 'decrc'
      }
      print(gene_type)
      
    } else{
      gene_type <- 'undefined'
    }
    
    
    gene_cluster[genes_use,'classify_gene'] <- gene_type
    print(paste0("Cluster ",cl," and gene type is: ",gene_type))
    print(" ")
  }
  
  
  
  write.csv(gene_cluster, paste0(save_dir, datatag,'_gene_clusters_classified.csv'), quote = F, row.names = F)
  return(gene_cluster)
}


gene_regression_SA609 <- function(exp_mean_df, gene_cluster, save_dir, datatag){
  # gene_cluster <- as.data.frame(gene_cluster)
  # exp_mean_df <- as.data.frame(exp_mean_df)
  epsilon <- 0.05 # error rate between different batches
  rownames(gene_cluster) <- gene_cluster$gene
  gene_cluster$classify_gene <- ''
  gene_cluster$slope <- NA
  for(cl in unique(gene_cluster$cluster)){
    print(paste0("Observe clone: ", cl))
    genes_use <- gene_cluster[gene_cluster$cluster==cl,]$gene
    tmp <- exp_mean_df[exp_mean_df$gene %in% genes_use,]
    # tmp <- exp_mean_df %>%
    #               dplyr::filter(gene %in% as.character(genes_use))
    tmp$treatment_status <- as.factor(tmp$treatment_status)
    dim(tmp)
    rg <- lm(exp_scaled ~ treatment_status, data = tmp)  
    saveRDS(rg, paste0(save_dir, datatag,'_cluster_',cl,'_lm.rds'))
    gene_cluster[genes_use,'slope'] <- rg$coefficients[2] 
    print(summary(rg))
    coeffs <- coef(rg)
    # Get means in each treatment
    ts_med <- tapply(tmp$exp_scaled, tmp$treatment_status, median)
    # print(ts_means)
    c2 <- coeffs[2] > 0
    c3 <- coeffs[3] > 0
    c4 <- coeffs[4] > 0
    c21 <- ts_med[2] - ts_med[1] >= -epsilon
    c32 <- ts_med[3] - ts_med[2] >= -epsilon
    c31 <- ts_med[3] - ts_med[1] >= -epsilon
    c41 <- ts_med[4] - ts_med[1] >= -epsilon
    c43 <- ts_med[4] - ts_med[3] >= -epsilon
    
    c12 <- ts_med[1] - ts_med[2] >= -epsilon
    c23 <- ts_med[2] - ts_med[3] >= -epsilon
    c13 <- ts_med[1] - ts_med[3] >= -epsilon
    c34 <- ts_med[3] - ts_med[4] >= -epsilon
    if(c4 & c21 & c31){
      if(c32 & c43){
        gene_type <- 'mono_incrc'
      } else{
        gene_type <- 'incrc'
      }
      print(gene_type)
    } else if(!c4 & c12 & c13){
      if(c23 & c34){
        gene_type <- 'mono_decrc'
      }else{
        gene_type <- 'decrc'
      }
      print(gene_type)
      
    } else{
      gene_type <- 'undefined'
    }
    
    
    gene_cluster[genes_use,'classify_gene'] <- gene_type
    print(paste0("Cluster ",cl," and gene type is: ",gene_type))
    print(" ")
  }
  
  
  
  write.csv(gene_cluster, paste0(save_dir, datatag,'_gene_clusters_classified.csv'), quote = F, row.names = F)
  return(gene_cluster)
  
}

plot_gene_clusters <- function(gexp_cls, tag, datatag, obs_clones, save_dir, clone_aware,
                               xl="Time point - treatment status", yl="Mean z-score value", trim_data=TRUE){
  print(summary(gexp_cls$mean_exp_scaled))
  if(trim_data){
    gexp_cls$mean_exp_scaled <- ifelse(as.numeric(gexp_cls$mean_exp_scaled)>-1.5,-1.5, gexp_cls$mean_exp_scaled)
  }
  nbclusters <- length(unique(gexp_cls$cluster))
  if(clone_aware){
    if(is.null(tag)){
      tag <- paste0(' in clone ',paste(obs_clones, collapse = '_'))
    }
    tag1 <- paste0('_clone_aware_',paste(obs_clones, collapse = '_'))
  }else{
    tag <- ''
    tag1 <- '_without_clone_aware'
  }
  plttitle <- paste0(datatag, ": top 500 resistant genes across time ",tag)
  # lowpass.spline <- smooth.spline(gexp_cls$treatment_status, gexp_cls$mean_exp_scaled, spar = 0.6) ## Control spar for amount of smoothing
  gexp_cls$cluster <- paste0('Cluster ',gexp_cls$cluster)
  pcm <- gexp_cls %>% 
    ggplot(aes(treatment_status, mean_exp_scaled)) +
    geom_boxplot() 
  
  pcm <- gexp_cls %>% 
    ggplot(aes(treatment_status, mean_exp_scaled)) +
    geom_line(aes(group = gene), alpha = 0.3) +  #geom_line  #
    # geom_boxplot() +
    # geom_line(predict(lowpass.spline, gexp_cls$treatment_status), col = "red", lwd = 2) + 
    geom_line(stat = "summary", fun = "median", colour = "brown", size = 1,
              aes(group = 1)) +
    facet_grid(cols = vars(cluster))
  
  pcm <- pcm + labs(x=xl, y=yl, title=plttitle)
  pcm <- pcm + theme(plot.title = element_text(color="black", size=13, hjust = 0.5),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.y = element_text(color="black", size=10),
                     axis.text.x = element_text(color="black", size=8, angle=90),
                     axis.title = element_text(color="black", size=12)
  )#egend.text = element_text(color="black", size=6)
  png(paste0(save_dir,datatag,tag1 ,"_gene_regression.png"), height =2*300, width= 2*100*nbclusters+250,res = 2*72)
  print(pcm)
  dev.off()
  
}


plot_gene_each_cluster <- function(gexp_cls, datatag, obs_clones, obs_clusters, save_dir,
                                   xl="Time point - treatment status", yl="Mean z-score value"){
  tag <- paste0(' in clone ',paste(obs_clones, collapse = ' '))
  # if(clone_aware){
  #   tag <- paste0(' in clone ',paste(obs_clones, collapse = ' '))
  #   tag1 <- paste0('_clone_aware_',paste(obs_clones, collapse = '_'))
  # }else{
  #   tag <- ''
  #   tag1 <- '_without_clone_aware'
  # }
  plttitle <- paste0("Resistant genes across treatment time points ",tag)
  # lowpass.spline <- smooth.spline(gexp_cls$treatment_status, gexp_cls$mean_exp_scaled, spar = 0.6) ## Control spar for amount of smoothing
  plot_ls <- list()
  for(cl in obs_clusters){
    pltsubtitle <- paste0("Genes Cluster ",cl)
    gexp_cls_tmp <- gexp_cls %>%
      dplyr::filter(cluster==cl)
    pcm <- gexp_cls_tmp %>% 
      ggplot(aes(treatment_status, mean_exp_scaled)) +
      geom_line(aes(group = gene, color=gene), alpha = 0.3) +  #geom_line  #
      # geom_smooth(method = "loess") +
      # geom_line(predict(lowpass.spline, gexp_cls$treatment_status), col = "red", lwd = 2) + 
      geom_line(stat = "summary", fun = "median", colour = "brown", size = 1,
                aes(group = 1)) #+
    # facet_grid(cols = vars(cluster))
    
    pcm <- pcm + labs(x=xl, y=yl, title=pltsubtitle)
    pcm <- pcm + theme(plot.title = element_text(color="black", size=13, hjust = 0.5),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                       axis.text.y = element_text(color="black", size=10),
                       axis.text.x = element_text(color="black", size=8, angle=90),
                       axis.title = element_text(color="black", size=12),
                       legend.text = element_text(color="black", size=5)
    )#
    plot_ls[[cl]] <- pcm
  }
  main_plot <- cowplot::plot_grid(plotlist = plot_ls,
                                  ncol = 1,
                                  align = 'v'
  )
  
  
  png(paste0(save_dir,datatag,tag ,"_gene_regression.png"), height = 2*600, width=2*1200,res = 2*72)
  print(main_plot)
  dev.off()
  
}


get_top_up_genes <- function(de_genes, ntop=NULL, minLogFC=0.25, FDR_cutoff=0.01, pValueThrs=0.5){
  de_genes <- de_genes %>%
    dplyr::filter(logFC>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
  de_genes <- de_genes[order(de_genes$logFC,decreasing = T),] 
  if(ntop > nrow(de_genes)){
    ntop = nrow(de_genes)
  }
  if(!is.null(ntop)){
    return(de_genes[1:ntop,])  
  }else{  # select all resistant genes with logFC > 0
    return(de_genes)
  }
  
}

downsample_cells <- function(meta_cells_df,downsample_ratio=0.5,
                             thres_small_clone=500){
  # meta_cells_df <- as.data.frame(meta_cells_df)
  ts <- unique(meta_cells_df$treatmentSt)
  thres_minority <- 20
  ext_cells <- c()
  for(t in ts){
    tmp <- meta_cells_df %>%
      dplyr::filter(treatmentSt==t)
    
    if(!is.null(tmp) & nrow(tmp)<thres_small_clone & nrow(tmp)>thres_minority){   #if clone contain less than 15 cells, remove this clone from consideration
      cells_to_sample <- tmp$cell_id
    } else if(!is.null(tmp) & nrow(tmp) >= thres_small_clone){
      cells_to_sample <- sample(tmp$cell_id, floor(nrow(tmp)*downsample_ratio),replace = FALSE)
    } else{
      cells_to_sample <- NULL
      
      if(!is.null(tmp)){
        print(paste0('DEBUG: double check this step, treatment st: ',t, ' nb cells: ',nrow(tmp)))
      } else{
        print(paste0('DEBUG: double check this step, treatment st: ',t))
      }
    }
    print(paste0('Extract ',length(cells_to_sample),' from treatment st: ', t))
    
    if(length(cells_to_sample)>0){
      ext_cells <- c(ext_cells, cells_to_sample)
    }
    
  }
  print(length(ext_cells))
  return(meta_cells_df[ext_cells,])
}


get_treatment_status <- function(status) {
  labels <- sapply(strsplit(status, "-"), function(x) {
    return(x[3])
  })
  return(as.character(labels))
}