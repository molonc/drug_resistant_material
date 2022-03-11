# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')

suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

get_earliest_principal_node <- function(cds, time_bin="UT", feature_use="treatmentSt"){
  cell_ids <- which(colData(cds)[, feature_use] == time_bin)
  # print(length(cell_ids))
  # cell_ids_p1 <- which(partitions(cds)==1)
  # cell_ids <- intersect(cell_ids, cell_ids_p1)
  print('Number of used cells:')
  print(length(cell_ids))
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
create_monocle_obj_from_sce <- function(tsce, save_dir=NULL, exp='logcounts', preprocess=F, save_data=F){
  expression_matrix <- as.matrix(assay(tsce, exp))
  cell_metadata <- as.data.frame(colData(tsce))
  gene_metadata <- as.data.frame(rowData(tsce))
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata=cell_metadata,
                           gene_metadata = gene_metadata)
  print(dim(cds))
  if(preprocess){
    cds <- preprocess_cds(cds, num_dim = 50, norm_method='none')  # data have already processed
    cds <- reduce_dimension(cds, preprocess_method = 'PCA')
  }
  if(save_data){
    print('Save data into folder...')
    saveRDS(cds, paste0(save_dir,'cds_Rx.rds'))
  }
  return(cds)
}  

plot_raw_data <- function(cds, save_dir, cls_resolution=0.3){
  cds <- cluster_cells(cds, resolution=cls_resolution, cluster_method="leiden", reduction_method = 'UMAP')
  p1 <- plot_cells(cds, color_cells_by = "partition",
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE)
  
  # cds <- learn_graph(cds)
  p2 <- plot_cells(cds,
                   color_cells_by = "treatmentSt",
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE)
  tag <- ''
  tag <- paste(unique(colData(cds)$treatmentSt),collapse = '_')
  p <- cowplot::plot_grid(p1, p2, ncol=2, align = 'h') #label=c('Drug Cycle Treatment','Untreated')
  png(paste0(save_dir, "umap_inputdata_",tag,".png"), height =2*350, width= 2*800,res = 2*72)
  print(p)
  dev.off()
  
  return(list(p=p, cds=cds))
}

downsample_data <- function(meta_data, downsample_ratio=0.5){
  thres_small_group <- 500
  thres_minority <- 20
  ext_cells <- c()
  for(c in unique(meta_data$treatmentSt)){
    
    tmp <- meta_data[meta_data$treatmentSt==c,]
    if(!is.null(tmp) & nrow(tmp)<thres_small_group & nrow(tmp)>thres_minority){   #if clone contain less than 15 cells, remove this clone from consideration
      cells_to_sample <- tmp$cell_id
    } else if(!is.null(tmp) & nrow(tmp) >= thres_small_group){
      cells_to_sample <- sample(tmp$cell_id, floor(nrow(tmp)*downsample_ratio),replace = FALSE)
    } else{
      cells_to_sample <- NULL
      
      if(!is.null(tmp)){
        print(paste0('DEBUG: double check this step, ts: ',c, ' nb cells: ',nrow(tmp)))
      } else{
        print(paste0('DEBUG: double check this step, ts: ',c))
      }
    }
    print(paste0('Extract ',length(cells_to_sample),' from treatmentSt: ', c))
    
    if(length(cells_to_sample)>0){
      ext_cells <- c(ext_cells, cells_to_sample)
    }
    
  }
  print(length(ext_cells))
  return(ext_cells)
}  

# FM is normalized data
preprocess_cds_v2 <- function(cds, FM, method ='PCA',
                           num_dim=50,
                           scaling = TRUE,
                           verbose=TRUE
                           ) {
  cds <- cds[rownames(FM),colnames(FM)]
  print(dim(cds))
  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
  
  if(method == 'PCA') {
    if (verbose) message("Remove noise by PCA ...")
    irlba_res <- monocle3:::sparse_prcomp_irlba(Matrix::t(FM),
                                     n = min(num_dim,min(dim(FM)) - 1),
                                     center = scaling, scale. = scaling)
    preproc_res <- irlba_res$x
    row.names(preproc_res) <- colnames(cds)
    
    irlba_rotation <- irlba_res$rotation
    row.names(irlba_rotation) <- rownames(FM)
    cds@preprocess_aux$gene_loadings <- irlba_rotation %*% diag(irlba_res$sdev)
    cds@preprocess_aux$prop_var_expl <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)
    
  } 
  print('Fill in data')
  row.names(preproc_res) <- colnames(cds)
  
  reducedDims(cds)[[method]] <- as.matrix(preproc_res)
  cds@preprocess_aux$beta = NULL
  
  cds
}  
# library(monocle3)
# monocle3::find_gene_modules
find_coregulated_genes_v2 <- function(cds, genes_use, output_dir, datatag, cls_res=0.05){
  gene_module_df <- monocle3::find_gene_modules(cds[genes_use,], resolution=cls_res, max_components = 15)
  gene_module_df <- monocle3::find_gene_modules(cds[rowData(cds)$gene_short_name %in% genes_use,], resolution=0.05, max_components = 15)
  data.table::fwrite(gene_module_df, paste0(output_dir,datatag,'gene_module_treated.csv.gz'))
  print(dim(gene_module_df))
  print(length(unique(gene_module_df$module)))
  unique(colData(cds_subset)$cluster_label)
  cds_subset <- cds[,!colData(cds)$cluster_label %in% c('Cls_6','Cls_4')]
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)), 
                                  cell_group=paste0(colData(cds_subset)$treatmentSt,'_',colData(cds_subset)$cluster_label))
  
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)), 
                                  cell_group=colData(cds_subset)$treatmentSt)
  agg_mat <- monocle3::aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df, max_agg_value = 3, min_agg_value = -3)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  
  # p4
  
  p <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
  png(paste0(output_dir, datatag,"_genes_modules_without_supermodules.png"), height =2*700, width= 2*700,res = 2*72)
  print(p)
  dev.off()
  saveRDS(agg_mat,paste0(output_dir,datatag,'_coregulated_genes_module_treated.rds'))
  
  # library(pheatmap)
  p1 <- pheatmap::pheatmap(agg_mat, cluster_rows=T, cluster_cols=TRUE,
                           scale="column", clustering_method="ward.D2",
                           fontsize=12)
  
  module_cls <- cutree(p1$tree_row, k = 2)
  cls <- data.frame(genes_super_module=module_cls, row.names = names(module_cls))
  cls$genes_super_module <- as.factor(cls$genes_super_module)
  
  p1 <- pheatmap::pheatmap(agg_mat, cluster_rows=T, cluster_cols=TRUE,
                           scale="column", clustering_method="ward.D2",
                           fontsize=14,
                           annotation_row = cls,
                           cutree_rows = 2, cutree_cols = 2)
  png(paste0(output_dir, datatag, "_genes_modules_treated.png"), height =2*850, width= 2*900,res = 2*72)
  print(p1)
  dev.off()
}  
find_coregulated_genes <- function(cds, cls_res=0.2, save_dir){
  dim(cds)
  # preprocess_mat <- cds@preprocess_aux$gene_loadings
  # dim(preprocess_mat)
  gene_module_df <- find_gene_modules(cds, resolution=cls_res, max_components = 15)
  # gene_module_df <- data.table::fread(paste0(save_dir,'gene_module.csv.gz')) %>% as.data.frame()
  print(dim(gene_module_df))
  unique(gene_module_df$module)
  metacells <- colData(cds) %>% as.data.frame()
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=metacells$treatmentSt) #[colnames(cds)] #partitions(cds)
  agg_mat <- aggregate_gene_expression_v2(cds, gene_module_df, cell_group_df, max_agg_value = 3, min_agg_value = -3)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  # colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
  print(dim(agg_mat))
  library(pheatmap)
  p1 <- pheatmap::pheatmap(agg_mat, cluster_rows=T, cluster_cols=TRUE,
                           scale="column", clustering_method="ward.D2",
                           fontsize=12)
  
  module_cls <- cutree(p1$tree_row, k = 2)
  cls <- data.frame(genes_super_module=module_cls, row.names = names(module_cls))
  cls$genes_super_module <- as.factor(cls$genes_super_module)
  
  cells_cls <- data.frame(ts=agg_mat@Dimnames[[2]])
  cells_cls$treatment <- ifelse(grepl('*T',cells_cls$ts),'Rx','UnRx')
  cells_cls <- cells_cls %>% tibble::column_to_rownames(var="ts")
  p1 <- pheatmap::pheatmap(agg_mat, cluster_rows=T, cluster_cols=TRUE,
                           scale="column", clustering_method="ward.D2",
                           fontsize=14,
                           annotation_row = cls,
                           annotation_col = cells_cls,
                           cutree_rows = 2, cutree_cols = 2)
  png(paste0(save_dir, "genes_modules.png"), height =2*850, width= 2*900,res = 2*72)
  print(p1)
  dev.off()
  
  
  cls$module <- rownames(cls)
  gene_module_df$module <- paste0('Module ',gene_module_df$module)
  gene_module_df <- gene_module_df %>% inner_join(cls, by=c('module'))
  dim(gene_module_df)
  data.table::fwrite(gene_module_df, paste0(save_dir,'gene_module.csv.gz'))
  # View(head(cells_cls))
  # p2 <- plot_cells(cds, 
  #            genes=gene_module_df %>% filter(module %in% c(2)),
  #            group_cells_by="cluster",
  #            color_cells_by="treatmentSt",
  #            show_trajectory_graph=FALSE)
  # obs_modules <- c("Module 6","Module 1","Module 16")
  # gene_module_ext <- gene_module_df %>%
  #   dplyr::filter(module %in% obs_modules)
  # # dim(gene_module_ext)
  # dim(gene_module_df)
  # data.table::fwrite(gene_module_ext, paste0(save_dir,'gene_module_increase_resistant.csv.gz'))
  # gene_module_ext$id
  # # View(gene_module_ext[1:5,])
  # # summary(as.factor(gene_module_df$module))
  # minexp <- min(SingleCellExperiment::counts(cds[rownames(cds) %in% gene_module_ext$id[1:6],]))
  # p <- plot_genes_violin(cds[rownames(cds) %in% gene_module_ext$id[1:6],], group_cells_by='treatmentSt', ncol=3, normalize=F, log_scale=F, min_expr = minexp) +
  #   theme(axis.text.x=element_text(angle=45, hjust=1))
  # png(paste0(output_dir, "genes_module6.png"), height =2*500, width= 2*800,res = 2*72)
  # print(p)
  # dev.off()
  # 
  
  saveRDS(agg_mat,paste0(save_dir,'coregulated_genes_module_treated.rds'))
  # data.table::fwrite(agg_mat, paste0(save_dir,'coregulated_genes_treated.csv'), quote=F)
  
}
extract_DE_genes_clusters <- function(cds, genes_df, output_dir, cls_ls, feature_use='cluster_label'){
  fm <- paste0('~',feature_use)
  print(fm)
  cds <- cds[as.character(genes_df$gene_id), colData(cds)$cluster_label %in% cls_ls]
  print(dim(cds))
  gene_fits <- fit_models(cds, model_formula_str=fm, cores = 5) #quasipoisson for large dataset but negbinomial is better, , expression_family="negbinomial"
  fit_coefs <- coefficient_table(gene_fits)
  # print(head(fit_coefs))
  # print(summary(as.factor(fit_coefs$term)))
  time_terms <- fit_coefs %>% filter(term != '(Intercept)')
  time_terms <- time_terms %>% filter (q_value < 0.05) %>%
    select(gene_short_name, term, q_value, estimate, gene_id)
  # time_terms$term[1:10]
  
  time_terms <- time_terms %>%
    dplyr::filter(gene_id %in% genes_df$gene_id)
  print(dim(time_terms))
  # time_terms <- time_terms[order(time_terms$estimate, decreasing=T),]
  data.table::fwrite(time_terms, paste0(output_dir,paste(cls_ls, collapse='_'),'_DE_genes.csv.gz'), quote=F)
  
  # signif_genes <- unique(time_terms$gene_short_name)
  # print(length(unique(time_terms$gene_short_name)))
  # View(time_terms[1:50,])
  # unique(time_terms$term)
  # time_terms$gene_short_name[1:4]
  res <- list(de_genes=time_terms,fit_coefs=fit_coefs)
  saveRDS(res, paste0(output_dir,paste(cls_ls, collapse='_'), "DE_genes_treated.rds"))
  # return(res)
}
# dim(cds_subset)
fit_regression_DE_genes <- function(cds, genes_df, output_dir, feature_use='drug_dose'){
  fm <- paste0('~',feature_use)
  print(fm)
  gene_fits <- fit_models(cds, model_formula_str=fm, cores = 5) #quasipoisson for large dataset but negbinomial is better, , expression_family="negbinomial"
  fit_coefs <- coefficient_table(gene_fits)
  print(head(fit_coefs))
  print(summary(as.factor(fit_coefs$term)))
  time_terms <- fit_coefs %>% filter(term != '(Intercept)')
  time_terms <- time_terms %>% filter (q_value < 0.05) %>%
    select(gene_short_name, term, q_value, estimate, gene_id)
  # time_terms$term[1:10]
  print(dim(time_terms))
  time_terms <- time_terms %>%
    dplyr::filter(gene_id %in% genes_df$gene_id)
  print(dim(time_terms))
  # colnames(time_terms)
  data.table::fwrite(time_terms, paste0(output_dir,'genes_signif_treated.csv.gz'), quote=F)
  dim(time_terms)
  # View(head(time_terms$term))
  summary(as.factor((time_terms$term)))
  # fit_coefs$p_value[1:10]
  dim(fit_coefs)
  dim(cds)
  time_terms <- time_terms[order(time_terms$estimate),]
  # signif_genes <- unique(time_terms$gene_short_name)
  # print(length(unique(time_terms$gene_short_name)))
  # View(time_terms[1:50,])
  # unique(time_terms$term)
  # time_terms$gene_short_name[1:4]
  p <- plot_genes_violin(cds[rowData(cds)$gene_short_name %in% time_terms$gene_short_name[1:6],], 
                         group_cells_by='treatmentSt', ncol=3) + #, normalize=T, log_scale=F
    theme(axis.text.x=element_text(angle=45, hjust=1))
  png(paste0(output_dir, "genes_violin_module2.png"), height =2*500, width= 2*800,res = 2*72)
  print(p)
  dev.off()
  res <- list(de_genes=time_terms,fit_coefs=fit_coefs)
  saveRDS(res, paste0(output_dir, "DE_genes_treated.rds"))
  return(res)
}

aggregate_gene_expression_v2 <- function(cds,
                                         gene_group_df = NULL,
                                         cell_group_df = NULL,
                                         norm_method=c("log", "binary",
                                                       "size_only"),
                                         pseudocount=1,
                                         scale_agg_values=TRUE,
                                         max_agg_value=3,
                                         min_agg_value=-3,
                                         exclude.na=TRUE){
  if (is.null(gene_group_df) && is.null(cell_group_df))
    stop("Error: one of either gene_group_df or cell_group_df must not be NULL")
  # agg_mat <- normalized_counts(cds, norm_method=norm_method,
  #                              pseudocount=pseudocount)
  agg_mat = SingleCellExperiment::counts(cds)
  # class(agg_mat)
  if (is.null(gene_group_df) == FALSE){
    gene_group_df <- as.data.frame(gene_group_df)
    gene_group_df <- gene_group_df[gene_group_df[,1] %in%
                                     fData(cds)$gene_short_name |
                                     gene_group_df[,1] %in%
                                     row.names(fData(cds)),,drop=FALSE]
    
    # Convert gene short names to rownames if necessary. The more
    # straightforward single call to recode took much longer.
    # Thanks to Christopher Johnstone who posted this on github.
    short_name_mask <- gene_group_df[[1]] %in% fData(cds)$gene_short_name
    if (any(short_name_mask)) {
      geneids <- as.character(gene_group_df[[1]])
      geneids[short_name_mask] <- row.names(fData(cds))[match(
        geneids[short_name_mask], fData(cds)$gene_short_name)]
      gene_group_df[[1]] <- geneids
    }
    
    # gene_group_df = gene_group_df[row.names(fData(cds)),]
    
    # FIXME: this should allow genes to be part of multiple groups. group_by
    # over the second column with a call to colSum should do it.
    agg_mat = as.matrix(monocle3:::my.aggregate.Matrix(agg_mat[gene_group_df[,1],],
                                                       as.factor(gene_group_df[,2]),
                                                       fun="mean"))
    if (scale_agg_values){
      agg_mat <- t(scale(t(agg_mat)))
      agg_mat[agg_mat < min_agg_value] <- min_agg_value
      agg_mat[agg_mat > max_agg_value] <- max_agg_value
    }
  }
  
  if (is.null(cell_group_df) == FALSE){
    
    cell_group_df <- as.data.frame(cell_group_df)
    cell_group_df <- cell_group_df[cell_group_df[,1] %in% row.names(pData(cds)),,
                                   drop=FALSE]
    agg_mat <- agg_mat[,cell_group_df[,1]]
    agg_mat <- monocle3:::my.aggregate.Matrix(Matrix::t(agg_mat),
                                              as.factor(cell_group_df[,2]),
                                              fun="mean")
    agg_mat <- Matrix::t(agg_mat)
  }
  
  if (exclude.na){
    agg_mat <- agg_mat[rownames(agg_mat) != "NA", colnames(agg_mat) != "NA",drop=FALSE]
  }
  return(agg_mat)
}


# gene_module_df[,2]
# rownames(preprocess_mat)[1:3]
aggregate_gene_expression_v3 <- function(preprocess_mat,
                                         gene_group_df = NULL,
                                         cell_group_df = NULL,
                                         norm_method=c("log", "binary",
                                                       "size_only"),
                                         pseudocount=1,
                                         scale_agg_values=TRUE,
                                         max_agg_value=3,
                                         min_agg_value=-3,
                                         exclude.na=TRUE){
  if (is.null(gene_group_df) && is.null(cell_group_df))
    stop("Error: one of either gene_group_df or cell_group_df must not be NULL")
  # agg_mat <- normalized_counts(cds, norm_method=norm_method,
  #                              pseudocount=pseudocount)
  agg_mat = preprocess_mat
  # class(agg_mat)
  if (is.null(gene_group_df) == FALSE){
    gene_group_df <- as.data.frame(gene_group_df)
    # gene_group_df <- gene_group_df[gene_group_df[,1] %in%
    #                                  fData(cds)$gene_short_name |
    #                                  gene_group_df[,1] %in%
    #                                  row.names(fData(cds)),,drop=FALSE]
    
    # Convert gene short names to rownames if necessary. The more
    # straightforward single call to recode took much longer.
    # Thanks to Christopher Johnstone who posted this on github.
    # short_name_mask <- gene_group_df[[1]] %in% fData(cds)$gene_short_name
    short_name_mask <- gene_group_df$id
    # if (any(short_name_mask)) {
    #   geneids <- as.character(gene_group_df[[1]])
    #   geneids[short_name_mask] <- row.names(fData(cds))[match(
    #     geneids[short_name_mask], fData(cds)$gene_short_name)]
    #   gene_group_df[[1]] <- geneids
    # }
    
    # gene_group_df = gene_group_df[row.names(fData(cds)),]
    
    # FIXME: this should allow genes to be part of multiple groups. group_by
    # over the second column with a call to colSum should do it.
    agg_mat = as.matrix(monocle3:::my.aggregate.Matrix(agg_mat[gene_group_df[,1],],
                                                       as.factor(gene_group_df[,2]),
                                                       fun="mean"))
    if (scale_agg_values){
      agg_mat <- t(scale(t(agg_mat)))
      print(max(agg_mat))
      print(min(agg_mat))
      agg_mat[agg_mat < min_agg_value] <- min_agg_value
      agg_mat[agg_mat > max_agg_value] <- max_agg_value
    }
  }
  
  if (is.null(cell_group_df) == FALSE){
    
    cell_group_df <- as.data.frame(cell_group_df)
    # cell_group_df <- cell_group_df[cell_group_df[,1] %in% row.names(pData(cds)),,
    #                                drop=FALSE]
    agg_mat <- agg_mat[,cell_group_df[,1]]
    agg_mat <- monocle3:::my.aggregate.Matrix(Matrix::t(agg_mat),
                                              as.factor(cell_group_df[,2]),
                                              fun="mean")
    agg_mat <- Matrix::t(agg_mat)
  }
  
  if (exclude.na){
    agg_mat <- agg_mat[rownames(agg_mat) != "NA", colnames(agg_mat) != "NA",drop=FALSE]
  }
  return(agg_mat)
}
# unique(tsce$treatmentSt)
# de_genes <- gene_module_df
# de_genes$ens_gene_id <- de_genes$id
# gmt_fn <- paste0(base_dir,"biodatabase/h.all.v7.0.symbols.gmt")
# gmt_fn <- paste0(base_dir,"biodatabase/c2.cp.kegg.v7.1.symbols.gmt")

get_logFC <- function(tsce, de_genes, genes_symb_df, gmt_fn){
  
  de_genes <- as.data.frame(de_genes)
  rownames(de_genes) <- de_genes$ens_gene_id
  sce_rx <- tsce[de_genes$ens_gene_id,grepl('*T',tsce$treatmentSt)]
  sce_unrx <- tsce[de_genes$ens_gene_id,grepl('*U',tsce$treatmentSt)]
  rxexp = apply(logcounts(sce_rx), 1, mean)
  unrxexp = apply(logcounts(sce_unrx), 1, mean)
  length(rxexp)
  logFC <- rxexp - unrxexp 
  max(logFC)
  min(logFC)
  logFC[1:3]
  deg_df <- data.frame(logFC=logFC,gene_symb=genes_symb_df[names(logFC),'gene_symb'])
  deg_df <- deg_df[order(abs(deg_df$logFC), decreasing=T),]
  deg_df <- deg_df[!duplicated(deg_df$gene_symb),]
  deg_stat <- deg_df$logFC
  names(deg_stat) <- deg_df$gene_symb   # TO DO: check duplicated genes, remove duplicated columns cause gene symbol are not unique in some case
  length(deg_stat)  
}

# h <- Hs.H
# ref_genes <- c()
# for(pt in names(h)){
#   ref_genes <- c(ref_genes, unlist(h[[pt]]))
# }
# sum(deg_df$gene_symb %in% ref_genes)
# dim(deg_df)

# data(examplePathways)
# data(exampleRanks)
# set.seed(42)
# View(head(exampleRanks))
# length(exampleRanks)
# fgseaRes <- fgsea(pathways = examplePathways, 
#                   stats    = exampleRanks,
#                   minSize  = 15,
#                   maxSize  = 500,
#                   nperm=1000)
# head(fgseaRes[order(pval), ])

get_significant_pathways <- function(obs_genes, save_dir, pathway_set='hallmark'){
  library(fgsea)
  ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  genome_genes_df <- read.csv(paste0(ref_dif, 'Symbol_ensembl.csv'), check.names = F, stringsAsFactors = F)  
  gmt_ls <- c("h.all.v7.0.symbols.gmt","c2.cp.kegg.v7.1.symbols.gmt","GO_c5.all.v7.1.symbols.gmt") #
  names(gmt_ls) <- c("hallmark","kegg","go") 
  gmt_fn <- paste0(ref_dif, gmt_ls[[pathway_set]])
  pathway_sets <- gmtPathways(gmt_fn) 
  
  stat <- tibble()
  for(pw in names(pathway_sets)){
    ci_tmp <- get_confident_interval_genes(obs_genes, pathway_sets[[pw]], pw, genome_genes_df)
    stat <- stat %>% bind_rows(ci_tmp)
  }
  dim(stat)
  
  
  stat <- stat %>%
    as.data.frame() %>%
    dplyr::filter(pval<0.05)
  stat <- stat[order(stat$pval, decreasing = F),]
  View(stat)
  data.table::fwrite(stat, paste0(save_dir,pathway_set,'_stats_downregulated.csv.gz'))
  
  # ref_genes <- ref_set
  # used_gene_set <- 'cisplatin_resistance'
  # used_gene_set <- 'oncogene_cosmic'
  # used_gene_set <- 'fitnessBF'
  # used_gene_set <- 'BroadSanger'
  # ci_tmp <- get_confident_interval_genes(obs_genes, ref_genes, used_gene_set, genome_genes_df)
  
}  
get_confident_interval_genes <- function(obs_genes, ref_genes, used_gene_set=' ', genome_genes_df=NULL){
  if(is.null(genome_genes_df)){
    ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
    genome_genes_df <- read.csv(paste0(ref_dif, 'Symbol_ensembl.csv'), check.names = F, stringsAsFactors = F)  
  }
  print(used_gene_set)
  # dim(genome_genes_df)
  genome_genes <- unique(genome_genes_df$Symbol)
  # print(length(genome_genes))
  ref_genes <- unique(ref_genes) # genes symbols
  print(length(ref_genes))
  CI_out <- get_bootstrap_stat(obs_genes, ref_genes, genome_genes,nsamples=10000) #nsamples=1000
  CI_out[['used_gene_set']] <- used_gene_set
  return(as.data.frame(CI_out))
}



get_bootstrap_stat <- function(our_obs, ref_genes, genome_genes, nsamples=1000){
  resamples <- lapply(1:nsamples, function(i) sample(genome_genes, size=length(our_obs), replace = T))
  count_occurance <- function(s){
    return(sum(s %in% ref_genes))
  }
  our <- sum(our_obs %in% ref_genes)
  print('Ours: ')
  print(our)
  occurs <- sapply(resamples, count_occurance)
  names(occurs) <- paste0('R',seq(1:length(occurs)))
  print('Randomized: ')
  print(summary(occurs))
  occurs['our'] <- our
  # length(occurs)
  r <- rank(occurs) #,ties.method = "max"
  pval <- (nsamples + 1 - r['our'])/nsamples
  return(list(CI=r['our'],pval=pval))
}

get_expr_ratio_across_time <- function(cds, exprs='counts', min_pct=0.1){
  # minimum counts per cell by gene
  # mnexprs<-10
  # # minimum cells with minimum counts by gene
  # mncells<-20
  # mnexprs_cell_by_gene <- apply(counts(cds), 1, function(v){v >= mnexprs})
  # mncells_by_gene <- apply(t(mnexprs_cell_by_gene), 1, function(v) {sum(v) >= mncells})
  # genes_with_mncells <- names(mncells_by_gene[mncells_by_gene])
  ts <- unique(colData(cds)$treatmentSt)
  genes_df <- data.frame(gene_id=rownames(cds))
  for(t in ts){
    zero_cbind <- DelayedArray::rowMeans(counts(cds[,colData(cds)$treatmentSt==t]) > 0)
    genes_df[,paste0('pct_expr_',t)] <- round(zero_cbind,2)
    if(t==ts[1]){
      genes_df$exp_quality <- round(zero_cbind,2) > min_pct
    }else{
      genes_df$exp_quality <- genes_df$exp_quality & (round(zero_cbind,2) > min_pct)
    }
  }
  # head(genes_df)
  print(summary(as.factor(genes_df$exp_quality)))
  # zero_cbind <- DelayedArray::rowMeans(assay(cds, exprs) > 0)
  # genes_use <- names(zero_cbind[zero_cbind >= min_pct])
  # genes_df <- data.frame(gene_id=genes_use, pct_expr=zero_cbind[genes_use])
  # genes_df <- genes_df[order(genes_df$pct_expr, decreasing=T),]
  return(genes_df)
}

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/encoder_trajectory/'
preprocess_mat <- data.table::fread(paste0(input_dir,'sctransform_SA535.csv.gz')) %>% as.data.frame()
dim(preprocess_mat)
rownames(preprocess_mat)[1:3]
colnames(preprocess_mat)[1:3]
preprocess_mat$V1 <- NULL
genes_df <- data.table::fread(paste0(input_dir, 'genes_ids_3000hvg.csv.gz'))
dim(genes_df)
rownames(preprocess_mat) <- genes_df$gene_id
# genes_use <- data.table::fread(paste0(input_dir, 'DE_genes_whole_dataset.csv')) %>% as.data.frame()
# branch_name = 'branches_8_5_10'
# genes_use <- data.table::fread(paste0(input_dir, 'DE_genes_DE_branches_8_5_10.csv')) %>% as.data.frame()
branch_name = 'branches_2_3_4'
genes_use <- data.table::fread(paste0(input_dir, 'DE_genes_DE_branches_2_3_4.csv')) %>% as.data.frame()
dim(genes_use)
sum(genes_use$beta_PDT<0)
genes_use <- genes_use %>%
  dplyr::filter(beta_PDT>0)


genes_use$gene_symb[1:5]
preprocess_mat <- preprocess_mat[as.character(genes_use$gene_symb),]
dim(preprocess_mat)


gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, resolution=0.04, max_components = 2)
data.table::fwrite(gene_module_df, paste0(input_dir,'gene_modules_',branch_name,'_positive.csv'))
summary(as.factor(gene_module_df$module))
gene_module_df <- gene_module_df %>%
  dplyr::filter(module %in% c(1,5))

genes_use$gene_symb
genes_use <- genes_use %>%
  dplyr::filter(gene_symb %in% gene_module_df$id)
# library(monocle3)
get_genes_modules <- function(preprocess_mat, reduction_method = c("UMAP"), max_components = 2, 
                              umap.metric = "cosine", umap.min_dist = 0.1, umap.n_neighbors = 15L, 
                              umap.fast_sgd = FALSE, umap.nn_method = "annoy", k = 20, 
                              leiden_iter = 1, partition_qval = 0.05, weight = FALSE, 
                              resolution = NULL, random_seed = 0L, cores = 5, verbose = F, 
                              ...) {
    method = "leiden"
    print('Run reduction features...')
    if (random_seed != 0L) 
      set.seed(random_seed)
    umap_res = uwot::umap(as.matrix(preprocess_mat), n_components = max_components, 
                          metric = umap.metric, min_dist = umap.min_dist, n_neighbors = umap.n_neighbors, 
                          fast_sgd = umap.fast_sgd, n_threads = cores, verbose = verbose, 
                          nn_method = umap.nn_method, ...)
    row.names(umap_res) <- row.names(preprocess_mat)
    colnames(umap_res) <- paste0("dim_", 1:ncol(umap_res))
    reduced_dim_res <- umap_res
    if (verbose) 
      message("Running leiden clustering algorithm ...")
    
    print('Running leiden clustering algorithm ...')
    # print(row.names(t(umap_res))[1:5])
    cluster_result <- monocle3:::leiden_clustering(data = reduced_dim_res, 
                                        pd = preprocess_mat, 
                                        k = k, weight = weight, num_iter = leiden_iter, resolution_parameter = resolution, 
                                        random_seed = random_seed, verbose = verbose, ...)
    print('Running partitions algorithm ...')
    cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g, 
                                            cluster_result$optim_res, partition_qval, verbose)
    partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
    names(partitions) <- row.names(reduced_dim_res)
    partitions <- as.factor(partitions)
    print('Get genes modules ...')
    gene_module_df <- tibble::tibble(id = row.names(preprocess_mat), 
                                     module = factor(igraph::membership(cluster_result$optim_res)), 
                                     supermodule = partitions)
    gene_module_df <- tibble::as_tibble(cbind(gene_module_df, 
                                              umap_res))
    return(gene_module_df)
  
  
}

preprocess_mat[1:5,1:5]
preprocess_mat$V1 <- NULL
max(preprocess_mat)
sum(preprocess_mat>3)
dim(preprocess_mat)

metacells <- data.table::fread(paste0(input_dir, 'grouping_SA535_v2.csv.gz')) %>% as.data.frame()
rownames(metacells) <- metacells$cell_id
metacells$treatmentSt
dim(preprocess_mat)
dim(gene_module_df)
dim(cell_group_df)
cell_group_df <- tibble::tibble(cell=colnames(preprocess_mat), 
                                cell_group=metacells[colnames(preprocess_mat),'treatmentSt']) #[colnames(cds)] #partitions(cds)
agg_mat <- aggregate_gene_expression_v3(preprocess_mat, gene_module_df, cell_group_df, max_agg_value = 4, min_agg_value = -4)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
# colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
print(dim(agg_mat))
# library(pheatmap)
p1 <- pheatmap::pheatmap(agg_mat, cluster_rows=T, cluster_cols=TRUE,
                         scale="column", clustering_method="ward.D2",
                         fontsize=12)
png(paste0(input_dir, "genes_modules_SA535_",branch_name,"_positive.png"), height =2*600, width= 2*800,res = 2*72)
print(p1)
dev.off()

