suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  # require("argparse")
  require("tidyr")
  require("SingleCellExperiment")
  require("stringr")
  require("scran")
  require("Seurat")
  require("data.table")
  require("clusterProfiler")
  require("org.Hs.eg.db")
  require("fgsea")
  require("DOSE")
  # require("enrichplot")
  require("ggplot2")
  require("ComplexHeatmap")
  require("annotables")
  require("dplyr")
  require("edgeR")
  # require("grid")
  
})
script_dir <- '/home/htran/Projects/farhia_project/rscript/pipeline'
source(paste0(script_dir, "/utils/deg_utils.R"))


get_treatment_status <- function(status) {
  labels <- sapply(strsplit(status, "-"), function(x) {
    return(x[3])
  })
  return(as.character(labels))
}

# datatag: 'SA535_CX', 'SA535_CY','SA1035', 'SA609'
edgeR_DE_analysis_by_clones_samples <- function(pair_groups, 
                                                sce_fn, output_file, datatag='SA', 
                                                min_features=1500, min_pct=0.1){
  # pair_groups <- read.csv(pair_groups_fn, header=T, check.names=F)
  # pair_groups$desc <- paste0(pair_groups$datatag,'_',
  #                            pair_groups$sample1,'-',pair_groups$clone1,'_',
  #                            pair_groups$sample2,'-',pair_groups$clone2)
  # pair_groups$desc <- paste0(pair_groups$datatag,'_',pair_groups$clone1,'-',pair_groups$clone2)
  
  observed_pair_groups <- pair_groups[pair_groups$datatag==datatag,'desc']
  rownames(pair_groups) <- pair_groups$desc
  print("Number of pair comparison: ")
  print(length(observed_pair_groups))
  if(length(observed_pair_groups)==0){
    stop("DEBUG: No comparison to do, double check input setting")
  }
  
  save_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir, recursive=T)
  }
  
  sce <- readRDS(sce_fn)
  print("Initialized sce")
  print(dim(sce))
  print(summary(as.factor(sce$clone)))
  print(summary(as.factor(sce$sample)))
  if(grepl('SA535',datatag)){
    meta_df <- read.csv('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/snakemake/metasample_SA535.csv', check.names = F)
    rownames(meta_df) <- meta_df$library_id
    sce$pdx <- meta_df[sce$library_id,'PDX']
    sce$treatmentSt <- meta_df[sce$library_id,'treatmentSt']
    if(datatag=='SA535_CY'){
      sce <- sce[,sce$pdx!='SA535_CX']
    }else{
      sce <- sce[,sce$pdx!='SA535_CY']
    }
  }
  # table(sce$sample,sce$clone)
  print(dim(sce))
  
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  sum(mito_genes==TRUE)
  
  ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
  sum(ribo_genes==TRUE)
  sce <- sce[(!mito_genes) & (!ribo_genes), ]
  print("Observed sce: ")
  print(paste0('Removing mito and ribo genes: ',dim(sce)[1],'_',dim(sce)[2]))
  print(sce$series[1:3])
  sce$treatmentSt <- get_treatment_status(sce$series)
  sce$treatmentSt[1:3]
  meta_data <- as.data.frame(colData(sce))
  # summary(as.factor(meta_data$clone))
  
  for(de in observed_pair_groups){
    cl1 <- unlist(strsplit(as.character(pair_groups[de,'clone1']), "_"))
    cl2 <- unlist(strsplit(as.character(pair_groups[de,'clone2']), "_"))
    # cells_use_g1 <- rownames(meta_data)[meta_data$clone %in% cl1 & meta_data$sample==pair_groups[de,'sample1']]
    # cells_use_g2 <- rownames(meta_data)[meta_data$clone %in% cl2 & meta_data$sample==pair_groups[de,'sample2']]
    cells_use_g1 <- rownames(meta_data)[meta_data$treatmentSt %in% cl1]
    cells_use_g2 <- rownames(meta_data)[meta_data$treatmentSt %in% cl2]
    
    if(length(cells_use_g1) < 100 || length(cells_use_g2) < 100){
      print("There are no cells or small nb cells 
            only which satisfy the input condition ")
      print(paste0("\n Observed clones: ",pair_groups[de,'clone1'],' vs ',pair_groups[de,'clone2'], 
                   "  nb cells in ",pair_groups[de,'clone1'], ":",length(cells_use_g1),
                   "  nb cells in ",pair_groups[de,'clone2'], ":",length(cells_use_g2)))
      
    } else{
      save_dir_pw <- paste0(save_dir,de,"/")
      
      if (!file.exists(save_dir_pw)){
        dir.create(save_dir_pw, recursive = T)
      }
      cells_use <- c(cells_use_g1, cells_use_g2)
      # cells_use <- rownames(meta_data)
      cat("\n DE analysis ", file = paste0(save_dir_pw,"de_analysis_log.txt"))
      print(length(cells_use))
      cat(paste0("\n Observed clones: ",pair_groups[de,'clone1'],' vs ',pair_groups[de,'clone2'], 
                 "  nb cells in ",pair_groups[de,'clone1'], ":",length(cells_use_g1),
                 "  nb cells in ",pair_groups[de,'clone2'], ":",length(cells_use_g2)), 
          file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      
      sce_tmp <- sce[,cells_use]
      print(dim(sce_tmp))
      pair_groups[de,'nbcells_g1'] <- length(cells_use_g1)
      pair_groups[de,'nbcells_g2'] <- length(cells_use_g2)
      
      
      # Similar to findMarkers func in Seurat, only test genes that are detected in a minimum fraction 
      # of min.pct cells in either of the two populations. Meant to speed up the function
      # by not testing genes that are very infrequently expressed. Default is 0.1
      print("Filtering by pct minimum fraction genes in each population")
      # zero_g1 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$clone %in% cl1], "counts") == 0)
      zero_g1 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$treatmentSt %in% cl1], "counts") == 0)
      
      obs_genes1 <- names(zero_g1[zero_g1 <= (1 - min_pct)])
      print(length(obs_genes1))
      # zero_g2 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$clone %in% cl2], "counts") == 0)
      zero_g2 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$treatmentSt %in% cl2], "counts") == 0)
      
      obs_genes2 <- names(zero_g2[zero_g2 <= (1 - min_pct)])
      print(length(obs_genes2))
      sce_tmp <- sce_tmp[intersect(obs_genes1, obs_genes2),]
      cat("\n Filtering by pct minimum fraction genes:\n ",file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      cat(dim(sce_tmp),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      sce_tmp <- sce_tmp[, sce_tmp$total_features_by_counts > min_features]
      cat("\n Filtering by total_features_by_counts greater than 1500:\n ",file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      cat(dim(sce_tmp),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      # sce_tmp$clone <- ifelse(sce_tmp$clone %in% cl1, paste0("2_",pair_groups[de,'clone1']),paste0("1_",pair_groups[de,'clone2']))
      sce_tmp$treatmentSt <- ifelse(sce_tmp$treatmentSt %in% cl1, paste0("2_",pair_groups[de,'clone1']),paste0("1_",pair_groups[de,'clone2']))
      cat(summary(as.factor(sce_tmp$treatmentSt)),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      
      minLogFC <- 0.25
      pValueThrs <- 0.05
      de_genes <- edgeR_de(sce_tmp, save_dir_pw)
      
      # de_genes <- de_genes[abs(de_genes$logFC)>minLogFC,]
      # cat(nrow(de_genes),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      # # Plot DE genes
      # plttitle <- paste0(pair_groups[de,'datatag'],":  ",pair_groups[de,'clone1']," versus ",pair_groups[de,'clone2'])
      # nbtopup <- 30
      # nbtopdown <- 30
      # markers_ls_upreg <- de_genes[de_genes$logFC>minLogFC,]
      # markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$logFC,decreasing = T),] 
      # cat(paste0('\n ',de), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      # cat(paste0('\n Nb resistant genes: ', nrow(markers_ls_upreg)),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      # pair_groups[de,'resistant_genes'] <- nrow(markers_ls_upreg)
      # 
      # # dim(markers_ls_upreg)
      # if(nrow(markers_ls_upreg) < nbtopup){
      #   nbtopup <- nrow(markers_ls_upreg)
      # }
      # markers_ls_upreg <- markers_ls_upreg[1:nbtopup,]
      # 
      # # Get top down-regulated genes
      # markers_ls_downreg <- de_genes[de_genes$logFC<(-minLogFC),]
      # markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$logFC,decreasing = F),] 
      # cat(paste0('\n Nb sensitive genes: ', nrow(markers_ls_downreg)),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      # pair_groups[de,'sensitive_genes'] <- nrow(markers_ls_downreg)
      # # dim(markers_ls_downreg)
      # if(nrow(markers_ls_downreg) < nbtopdown){
      #   nbtopdown <- nrow(markers_ls_downreg)
      # }
      # markers_ls_downreg <- markers_ls_downreg[1:nbtopdown,]
      # topGenes <- c(as.character(markers_ls_upreg$gene_symbol),as.character(markers_ls_downreg$gene_symbol))
      # # rownames(markers_ls_tmp) <- markers_ls_tmp$gene_symb
      # plot_DE_genes_edgeR(de_genes, topGenes, nrow(de_genes), 
      #                     0.01, 0.25, plttitle, save_dir_pw, legendVisible=F)
    }  
  }
  
  
  # write.csv(pair_groups,file=pair_groups_fn, quote=F, row.names=F)
}



datatag <- 'SA609'
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6'
output_file <- paste0(input_dir,'/SA609_rna/deg_analysis/SA609_deg_pathway.rds')
save_dir <- paste0(input_dir,'/',datatag,'_rna/deg_analysis/')
pair_groups_fn <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/deg_analysis/pair_groups_SA1035_SA609.csv')
sce_fn <- paste0(results_10x_dir,'/total_sce_clones.rds')
pair_groups <- data.frame(clone1=c('UTTTT','UTTT','UTT'),
                          clone2=c('UTTT','UTT','UT'))
pair_groups$datatag <- rep(datatag, nrow(pair_groups))
pair_groups$desc <- paste0(pair_groups$datatag,'_',pair_groups$clone1,'-',pair_groups$clone2)
rownames(pair_groups) <- pair_groups$desc
pair_groups$time <- c('X8','X7','X6')
dim(pair_groups)

pair_groups1 <- data.frame(clone1=c('UUUUU','UUUU','UUU','UU'),
                          clone2=c('UUUU','UUU','UU','U'))
pair_groups1$datatag <- rep(datatag, nrow(pair_groups1))
pair_groups1$desc <- paste0(pair_groups1$datatag,'_',pair_groups1$clone1,'-',pair_groups1$clone2)
rownames(pair_groups1) <- pair_groups1$desc
pair_groups1$time <- c('X8','X7','X6','X5')
dim(pair_groups1)
edgeR_DE_analysis_by_clones_samples(pair_groups1, sce_fn, output_file, datatag, 1500, 0.1)

edgeR_DE_analysis_by_clones_samples(pair_groups, sce_fn, output_file, datatag, 1500, 0.1)


base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
de_x3 = read.csv(paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTT_R_UUUU_H/signif_genes.csv'))
de_x4 = read.csv(paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTTT_R_UUUUU_H/signif_genes.csv'))
observed_genes <- intersect(de_x3$ensembl_gene_id, de_x4$ensembl_gene_id)
length(observed_genes)
treated_gene_df <- genes_regression(pair_groups, observed_genes, save_dir)
untreated_gene_df <- genes_regression(pair_groups1, observed_genes, save_dir)

dim(treated_gene_df1)
dim(untreated_gene_df1)
length(unique(treated_gene_df$gene_id))
treated_gene_df1 <- treated_gene_df[!duplicated(treated_gene_df$gene_id),]
untreated_gene_df1 <- untreated_gene_df[!duplicated(untreated_gene_df$gene_id),]
colnames(untreated_gene_df1)
untreated_gene_df1 <- untreated_gene_df1 %>%
                      dplyr::select(gene_id, gene_type, slope)

untreated_gene_df1 <- untreated_gene_df1 %>%
  dplyr::rename(gene_type_untreated=gene_type, slope_untreated=slope)

treated_gene_df1 <- treated_gene_df1 %>% inner_join(untreated_gene_df1, by=c("gene_id"))
summary(as.factor(treated_gene_df1$gene_type))
table(treated_gene_df1$gene_type, treated_gene_df1$gene_type_untreated)


input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
datatag <- 'SA1035'
results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA1035-v6'
output_file <- paste0(input_dir,'/',datatag,'_rna/deg_analysis/',datatag,'_deg_pathway.rds')
sce_fn <- paste0(results_10x_dir,'/total_sce_clones.rds')


pair_groups <- data.frame(clone1=c('UTTTT','UTTT','UTT'),
                          clone2=c('UTTT','UTT','UT'))
pair_groups$datatag <- rep(datatag, nrow(pair_groups))
pair_groups$desc <- paste0(pair_groups$datatag,'_',pair_groups$clone1,'-',pair_groups$clone2)
rownames(pair_groups) <- pair_groups$desc
pair_groups$time <- c('X8','X7','X6')
dim(pair_groups)
edgeR_DE_analysis_by_clones_samples(pair_groups, sce_fn, output_file, datatag, 1500, 0.1)



# Read csv files of DE genes
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
de_x3 = read.csv(paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/SA1035_UTTT_H_UUUU_E/signif_genes.csv'))
de_x4 = read.csv(paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/SA1035_UTTTT_H_UUUUU_E/signif_genes.csv'))
observed_genes <- intersect(de_x3$ensembl_gene_id, de_x4$ensembl_gene_id)
length(observed_genes)
# observed_genes[1:3]
#Have logFC values, need to load data, and do linear regression
save_dir <- paste0(input_dir,'/',datatag,'_rna/deg_analysis/')


genes_regression <- function(pair_groups, observed_genes, save_dir){
  print(length(observed_genes))
  de_ls <- list()
  for(de in pair_groups$desc){
    de_genes <- read.csv(paste0(save_dir, de, '/edgeR_significant_genes.csv'))
    
    observed_genes <- intersect(observed_genes, de_genes$gene_id)
    de_ls[[as.character(de)]] <- de_genes
  }
  print(length(observed_genes))
  
  for(de in pair_groups$desc){
    de_genes <- de_ls[[as.character(de)]]
    de_genes <- de_genes %>%
      dplyr::filter(gene_id %in% observed_genes)
    de_genes$desc <- de
    de_genes$time <- pair_groups[de,'time']
    de_ls[[as.character(de)]] <- de_genes
  }
  
  genes_df <- do.call(rbind, de_ls)
  class(genes_df)
  genes_df$mean <- genes_df$logFC
  # rownames(genes_df) <- genes_df$desc
  dim(genes_df)
  # View(head(genes_df))
  genes_df$gene_type <- ''
  for(g in unique(genes_df$gene_id)) {
    # print(gene)
    ## untreated
    # genelm <- genedf[genedf$name==gene & genedf$condition %in% c("UUU","UUUU","UUUUU"),]
    # linearMod <- lm(mean ~ time, data=genelm)
    # genedf[genedf$name==gene & genedf$condition %in% c("UUU","UUUU","UUUUU"),"slope"] <- linearMod$coefficients[2]
    ## treated
    # genelm <- genedf[genedf$name==gene & genedf$condition %in% c("UTT","UTTT","UTTTT"),]  
    genelm <- genes_df[genes_df$gene_id==g,]  
    linearMod <- lm(mean ~ time, data=genelm)
    genes_df[genes_df$gene_id==g,"slope"] <- linearMod$coefficients[2]  
    genes_df[genes_df$gene_id==g,"gene_type"] <- detect_gene_trend(genelm, linearMod)
  }
  return(genes_df)
}
# View(head(genes_df))
detect_gene_trend <- function(genelm, linearMod){
  # gene_cluster <- as.data.frame(gene_cluster)
  # exp_mean_df <- as.data.frame(exp_mean_df)
  epsilon <- 0.1 # error rate between different batches
  
  coeffs <- coef(linearMod)
  # Get means in each treatment
  # ts_med <- tapply(tmp$exp_scaled, tmp$treatment_status, mean)
  ts_med <- c(genelm[genelm$time=='X6','mean'],genelm[genelm$time=='X7','mean'],genelm[genelm$time=='X8','mean'])
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
  return(gene_type)
}






