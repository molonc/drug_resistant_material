
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(RColorBrewer)
  # library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(edgeR)
})

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/TreeAlign_testing/'
# dir.create(output_dir)
datatag <- 'SA609'

fns <- list.files(input_dir)
fns <- fns[grepl('.csv',fns)]

obs_samples <- c('SA609X3XB01584','SA609X4XB03080','SA609X5XB03223',
                 
                 'SA609X6XB03447','SA609X7XB03554','SA609X4XB003083',
                 
                 'SA609X5XB03230','SA609X6XB03404','SA609X7XB03505',
                 
                 'SA609X5XB03231','SA609X6XB03401','SA609X7XB03510',
                 
                 'SA609X4XB03083','SA609X4XB003083') # note: SA609X4XB03083 and SA609X4XB003083
length(obs_samples)
obs_samples <- c('SA609X5XB03223','SA609X7XB03554')
fns <- fns[fns %in% paste0(obs_samples,'.csv')]
clonealign_df <- tibble::tibble()
for(f in fns){
  print("-------------------------------------")
  print(f)
  tmp <- data.table::fread(paste0(input_dir,f))
  print(summary(as.factor(tmp$clone)))
  clonealign_df <- dplyr::bind_rows(clonealign_df, tmp)
  print("")
}
dim(clonealign_df)

obs_fns <- c('SA609X5XB03223','SA609X7XB03554')
clonealign_df$id[1:3]
clonealign_df <- clonealign_df %>% 
  dplyr::filter(id %in% obs_fns)

dim(clonealign_df)
clonealign_df$library_id <- gsub('.cache/','',clonealign_df$Sample)
clonealign_df$library_id <- gsub('/filtered_feature_bc_matrix','',clonealign_df$library_id)
print(summary(as.factor(clonealign_df$library_id)))

clonealign_df <- clonealign_df %>% 
  dplyr::filter(id == 'SA609X7XB03554')
unique(clonealign_df$id)
colnames(clonealign_df)
sce <- readRDS(paste0(input_dir, obs_fns[1],'.rdata'))
dim(sce)

output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/TreeAlign_testing/TreeAlign_result_thrs08/"
clones_df1 <- data.table::fread(paste0(output_dir, "SA609X5XB03223_clone_assign.csv.gz"))
print(summary(as.factor(clones_df1$clone_id)))

## clonealign results
# [1] "SA609X5XB03223.csv"
# B    C    D    H 
# 237  157 1032 1347 

## TreeAlign results
# [1] "SA609X5XB03223.csv"
#      B   C   D   G   H 
# 26  76 821 547 702 601

# [1] "-------------------------------------"
## clonealign results
# [1] "SA609X7XB03554.csv"
# B    C    D    G    H 
# 9 1013   23  117  581 

## TreeAlign results
#  C    H 
# 1447  296 

clones_df2 <- data.table::fread(paste0(output_dir, "SA609X7XB03554_clone_assign.csv.gz"))
print(summary(as.factor(clones_df2$clone_id)))
clones_df2$cell_id[1]
dim(clonealign_df)
dim(clones_df2)
clonealign_df$cell_id[1]
clones_df2$cell_id <- paste0('SA609X7XB03554_',clones_df2$cell_id)
clonealign_df <- clonealign_df %>%
  left_join(clones_df2, by='cell_id')
table(clonealign_df$clone, clonealign_df$clone_id)
# cols: treealign, rows: clonealign
#      C    H
# B    8    1
# C 1007    6
# D   21    2
# G   17  100
# H  394  187





get_DE_genes <- function(sce, cut_off_overall=0.1){
  cells_clone_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cell_clones/SA609_cell_clones_metadata.csv'
  cnv_fn <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones_v3/mapped_wholedata_SA609_grch38_v2.csv.gz'
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/TreeAlign_testing/'
  input_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/"
  output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/TreeAlign_testing/TreeAlign_result_thrs08/"
  # sids <- c('SA609X5XB03223','SA609X7XB03554')
  # sids <- c('SA609X3XB01584','SA609X4XB003083','SA609X4XB03080',
  #           'SA609X5XB03230','SA609X5XB03231','SA609X6XB03447')
  # sid <- 'SA609X5XB03223'
  # sid <- 'SA609X7XB03554'
  sid <- 'SA609X6XB03447'
  sce_fn <- paste0(input_dir, sid,'.rdata')
  sce <- readRDS(sce_fn)
  dim(sce)
  assayNames(sce)
  cut_off_overall=0.1 # 10% expressed --> get DE genes, reduce execution time here
  if(cut_off_overall > 0){
    zero_cbind <- DelayedArray::rowMeans(assay(sce, 'counts') == 0)
    sce <- sce[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
  }
  
  # norm_sce <- sce_normalize_size_factors(sce, min_size=300, rlog=FALSE, exprs="counts")
  # dim(norm_sce)
  norm_sce <- sce
  clones_df <- data.table::fread(paste0(output_dir, sid,"_clone_assign.csv.gz"))
  print(summary(as.factor(clones_df$clone_id)))
  
  clones_df <- clones_df %>%
    dplyr::mutate(
      clone_id=case_when(
        clone_id=='' ~ 'None',
        TRUE ~ clone_id)
    )
  
  colnames(norm_sce) <- norm_sce$Barcode
  rownames(norm_sce) <- rowData(norm_sce)$ID
  
  norm_sce@colData <- as.data.frame(colData(norm_sce)) %>%
    left_join(clones_df, by=c('Barcode'='cell_id'))%>%
    DataFrame(check.names = FALSE)
  summary(as.factor(norm_sce$clone_id))
  
  
  obs_clones <- unique(clones_df$clone_id)
  obs_clones <- obs_clones[obs_clones != 'None']
  obs_clones
  
  # obs_clones <- c('C','D','G','H')
  norm_sce <- norm_sce[,norm_sce$clone_id %in% obs_clones]
  dim(norm_sce)
  
  
  de_genes <- edgeR_de(norm_sce, save_dir, use_normcounts=F)
  dim(de_genes)
  head(de_genes)
  de_genes$gene_id[1]
  sum(abs(de_genes$logFC)>=0.25)
  
  # map DE genes back to gene score
  gene_scores <- data.table::fread(paste0(output_dir, sid,"_gene_type_score.csv.gz"))
  dim(gene_scores)
  cnv <- data.table::fread(paste0(save_dir, sid,"_clones_cnv.csv.gz"))
  head(cnv)
  sum(cnv$ensembl_gene_id %in% de_genes$gene_id)
  cnv <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% de_genes$gene_id)
  dim(cnv)
  cnv1 <- cnv %>%
    dplyr::select(all_of(paste0('cell_',obs_clones)))
  colnames(cnv1)
  var_genes <- sapply(seq(1:dim(cnv1)[1]), function(x){
    var(as.numeric(cnv1[x,]))
  })
  names(var_genes) <- cnv$ensembl_gene_id
  var_genes <- var_genes[var_genes>0]
  print(sum(var_genes>0))
  cnv <- cnv %>%
    dplyr::filter(ensembl_gene_id %in% names(var_genes))
  dim(cnv)
  de_genes <- de_genes %>%
    dplyr::mutate(
      gene_type=case_when(
        gene_id %in% cnv$ensembl_gene_id ~ 'cis',
        TRUE ~ 'trans')
    )
  de_genes <- de_genes %>%
    dplyr::left_join(gene_scores, by=c('gene_id'='gene'))
  
  de_genes <- de_genes %>%
    dplyr::mutate(
      gene_type_score=case_when(
        is.na(gene_type_score) ~ 0,
        TRUE ~ gene_type_score)
    )
  
  gt <- summary(as.factor(de_genes$gene_type))
  gt
  
  
  trans_acc <- 100*sum(de_genes$gene_type=='trans' & de_genes$gene_type_score<0.6) / gt['trans'] #65% correct
  # sum(de_genes$gene_type=='trans' & de_genes$gene_type_score>=0.6) # 35% not correct
  cis_acc <- 100*sum(de_genes$gene_type=='cis' & de_genes$gene_type_score>=0.6)/gt['cis']
  desc_acc <- paste0('Prediction accuracy: Trans: ',round(trans_acc),'%', ',  Cis: ',round(cis_acc),'%')
  desc_acc
  dim(de_genes)
  
  de_genes$gene_type <- as.factor(de_genes$gene_type)
  
  p <- ggplot(de_genes, aes(x=logFC, y=gene_type_score, colour=gene_type)) + 
    geom_jitter(size=1) + 
    geom_hline(yintercept = 0.6) + 
    theme_bw() + 
    labs(title = paste0('TreeAlign - ',sid),
         subtitle = desc_acc) + 
    theme(title = element_text(size=14, color='black'),
          axis.text = element_text(size=13, color='black'))
  
  p <- ggExtra::ggMarginal(p, margins = "y", size = 3, type = "histogram", 
                           groupColour = TRUE, groupFill = TRUE)
  
  png(paste0(save_dir,sid,"_gene_scores_cis_trans.png"), height = 2*500, width=2*700,res = 2*72)
  print(p)
  dev.off()
  
}
sce_normalize_size_factors <- function(sce, min_size=100, rlog=FALSE, exprs="counts"){
  # library(scater)
  ## BiocManager::install('scran')
  # library("scran")
  
  print("Quick clustering")
  if(min_size < dim(sce)[2]){
    qclust <- scran::quickCluster(sce, min.size = min_size, assay.type=exprs)
    print("Compute sum factors")
    sce <- scran::computeSumFactors(sce, clusters = qclust, assay.type=exprs)
    sce$size_factor <- sizeFactors(sce)
    # print(sce$size_factor[1:5])
    print("Normalize data")
    # scater version 1.14.6
    
    # count --> normcounts
    # normcounts --> logcounts 
    # String containing an assay name for storing the output normalized values. 
    # Defaults to "logcounts" when log=TRUE and "normcounts" otherwise.
    sce_normalized <- scater::logNormCounts(sce, log=rlog, exprs_values=exprs, size_factors=sce$size_factor)
    
  } else{
    print("Small nb cells in this library")
    sce_normalized <- logNormCounts(sce, log=rlog, exprs_values=exprs, size_factors=NULL)
  }
  # ?logNormCounts
  # scater version < 1.14.6
  # sce_normalized <- normalize(sce, return_log=rlog, exprs_values=exprs)
  print(assayNames(sce_normalized))
  return(sce_normalized)
}

edgeR_de <- function(sce_de, save_dir, use_normcounts=T){
  
  print("Filtering data...")
  # sce_de <- sce_de[, sce_de$total_features_by_counts > 1500]
  # rs <- rowSums(as.matrix(counts(sce_de)))
  # qplot(rs, log='x') + geom_vline(xintercept = 100)
  
  # sce_de <- sce_de[rowSums(as.matrix(counts(sce_de))) > 100, ]
  print(dim(sce_de))
  if(use_normcounts){
    mycounts <- as.matrix(normcounts(sce_de))    # zeros are good
    print(mycounts[1:3,1:3])
  }else{
    mycounts <- as.matrix(counts(sce_de))    # zeros are good  
  }
  
  print("Create DGE edgeR object...")
  # dge <- DGEList(counts=mycounts, group=sce_de$clone)
  dge <- DGEList(counts=mycounts, group=sce_de$clone_id)
  
  print("DE Analysis...")
  # which vs which ???
  # design <- model.matrix(~ clone, data = colData(sce_de))
  design <- model.matrix(~ clone_id, data = colData(sce_de))
  # This describes the edgeR user manual
  # http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  
  dge <- edgeR::estimateDisp(dge, design = design)
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit)
  tt <- edgeR::topTags(qlf, n = Inf)
  
  print("Generate output...")
  tt <- as.data.frame(tt) %>% 
    rownames_to_column("gene_id")
  # tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$ID
  tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$Symbol
  
  # Saving the logFC file
  # tt <- tt[tt$FDR<0.01 & tt$PValue<0.05 & abs(tt$logFC)>0.25,] # 
  # write.csv(tt, file=paste0(save_dir, 'edgeR_significant_genes.csv'), row.names=FALSE, quote=FALSE)
  
  data.table::fwrite(tt, paste0(save_dir, 'edgeR_significant_genes.csv.gz'))
  print("With FDR<0.01 and PValue<0.05 and abs(tt$logFC)>0.25, number of significant genes is: ")
  tt <- tt %>%
    dplyr::filter(abs(logFC)>0.25 & PValue<=0.05 & FDR<=0.01)
  print(dim(tt))
  # print("With threshold logFC>0.25, number of significant genes is: ")
  # print(nrow(tt[abs(tt$logFC)>0.25,]))
  return(tt)
}



load_cnv_TreeAlign <- function(cnv_fn, sid, input_dir, save_dir, min_cells_per_clone=15){
  cnv_fn <- 'mapped_wholedata_SA535_v2.csv.gz'
  cnv <- data.table::fread(cnv_fn)
}  
# [1] "-------------------------------------"
# [1] "SA609X3XB01584.csv"
# B   C   D 
# 15  99 160 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X4XB003083.csv"
# B   C   D   R 
# 185   5   6 404 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X4XB03080.csv"
# B        C_D          G          H unassigned 
# 317       3977          1        572          2 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X5XB03223.csv"
# B    C    D    H 
# 237  157 1032 1347 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X5XB03230.csv"
# B    R 
# 39 2380 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X5XB03231.csv"
# B    R 
# 2472  328 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X6XB03401.csv"
# R 
# 1844 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X6XB03404.csv"
# R 
# 3047 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X6XB03447.csv"
# B    C    D    H 
# 708 1132  104  506 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X7XB03505.csv"
# R 
# 1934 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X7XB03510.csv"
# R 
# 2372 
# [1] ""
# [1] "-------------------------------------"
# [1] "SA609X7XB03554.csv"
# B    C    D    G    H 
# 9 1013   23  117  581 
# [1] ""
#  

plot_umap <- function(umap_df, cols_use, datatag, output_dir, plottitle='', plotlegend=F){
  xl <- c(min(umap_df$UMAP_1),max(umap_df$UMAP_1))
  yl <- c(min(umap_df$UMAP_2),max(umap_df$UMAP_2))
  
  my_font <- "Helvetica"
  
  
  umap_df <- umap_df %>%
    dplyr::mutate(clone=replace(clone, is.na(clone),'None'))
  # umap_df$treatmentSt <- umap_df$treat
  tmp <- umap_df
  
  if(dim(tmp)[1]==0){
    print('Double check input passage info, exit')
    print('No cell at this passage')
    return(NULL)
  }
  
  obs_clones <- sort(unique(tmp$clone))
  unassign_clones <- c('unassigned','Unassigned','Un','un','None')
  obs_clones <- obs_clones[!obs_clones %in% unassign_clones]
  # print(obs_clones)
  # if(is.null(obs_clones)){
  #   
  # }
  cols_use <- cols_use[obs_clones]
  other_clones <- '#F0F0F0'
  names(other_clones) <- 'Others'
  cols_use <- c(cols_use,other_clones)
  # print(cols_use)
  # umap_df <- umap_df %>%
  #   dplyr::mutate(clone=replace(clone, !clone %in% obs_clones,'Others'))
  
  rownames(umap_df) <- umap_df$cell_id
  
  
  umap_df$clone <- factor(umap_df$clone, levels = sort(unique(umap_df$clone)))
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) + 
    geom_point(color='#e0e0e0') +  # grey color for all cells landscape displaying in background
    geom_point(data=tmp, aes(color=clone), size=0.7) + 
    scale_color_manual(values=cols_use, name='') + 
    labs(title = plottitle) +
    xlim(xl[1], xl[2]) + 
    ylim(yl[1], yl[2])
  
  p <- p + thesis_theme
  # lg <- cowplot::get_legend(p + guides(color = guide_legend(nrow = 1, 
  #                                                           title.position = "left", 
  #                                                           override.aes = list(size=2))))
  # plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  if(!plotlegend){
    lg_pos <- "none"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }else{
    lg_pos <- "right"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }
  
  results <- list(p=p, cols_use=cols_use, df=tmp)
  # basename <- paste0(datatag,"_", gsub(' ','',plottitle))
  # saveRDS(results, paste0(output_dir, basename, ".rds"))
  # 
  # png(paste0(output_dir,basename,".png"), height = 2*350, width=2*500,res = 2*72)
  # print(p)
  # dev.off()
  return(results)
  
}

plot_SA609 <- function(output_dir){
  source('/home/htran/Projects/farhia_project/rnaseq/drug_manuscript/viz_umap_figs/viz_umaps.R')
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/figs_rna_testTreeAlign/'
  # dir.create(output_dir)
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
  
  datatag <- 'SA609'
  basename <- datatag
  clonealign_stat <- data.table::fread(paste0(dirname(output_dir),'/clonealign_labels.csv'))
  clonealign_stat <- clonealign_df
  # unique(clonealign_stat$unique_clone)
  # clonealign_stat <- clonealign_stat %>%
  #   dplyr::filter(datatag==basename)%>%
  #   dplyr::select(cell_id, unique_clone)%>%
  #   dplyr::rename(clone=unique_clone)
  # colnames(clonealign_stat)
  # rm(umap_df)
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv.gz')) %>% as.data.frame()
  dim(umap_df)
  umap_df$clone <- NULL
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  summary(as.factor(umap_df$clone))
  summary(as.factor(umap_df$clone_id))
  umap_df$clone_clonealign <- umap_df$clone
  cols_use <- get_color_clones(datatag)
  res <- plot_umap(umap_df, cols_use, datatag, output_dir, plottitle='CloneAlign SA609X7XB03554', plotlegend=T)
  res$p
  umap_df$clone <- umap_df$clone_id
  res_ta <- plot_umap(umap_df, cols_use, datatag, output_dir, plottitle='TreeAlign SA609X7XB03554', plotlegend=T)
  res_ta$p
  p_total <- cowplot::plot_grid(res$p, res_ta$p)  
  dir.create(output_dir)
  png(paste0(output_dir,datatag,".png"), height = 2*400, width=2*1000,res = 2*72)
  print(p_total)
  dev.off()
}