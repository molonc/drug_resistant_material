edgeR_de <- function(sce_de, save_dir){
  
  print("Filtering data...")
  sce_de <- sce_de[, sce_de$total_features_by_counts > 1500]
  # rs <- rowSums(as.matrix(counts(sce_de)))
  # qplot(rs, log='x') + geom_vline(xintercept = 100)
  
  sce_de <- sce_de[rowSums(as.matrix(counts(sce_de))) > 100, ]
  print(dim(sce_de))
  mycounts <- as.matrix(counts(sce_de))    # zeros are good
  print("Create DGE edgeR object...")
  # dge <- DGEList(counts=mycounts, group=sce_de$clone)
  dge <- DGEList(counts=mycounts, group=sce_de$treatmentSt)
  dge <- calcNormFactors(dge, method = "TMM")
  print("DE Analysis...")
  # which vs which ???
  # design <- model.matrix(~ clone, data = colData(sce_de))
  design <- model.matrix(~ treatmentSt, data = colData(sce_de))
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
  write.csv(tt, file=paste0(save_dir, 'edgeR_significant_genes.csv'), row.names=FALSE, quote=FALSE)
  print("With FDR<0.01 and PValue<0.05 and abs(tt$logFC)>0.25, number of significant genes is: ")
  print(dim(tt))
  # print("With threshold logFC>0.25, number of significant genes is: ")
  # print(nrow(tt[abs(tt$logFC)>0.25,]))
  return(tt)
}



run_edgeR <- function(sce, pair_groups_fn, save_dir,
                      datatag='SA', 
                      min_features=1500, min_pct=0.1){
  
  pair_groups <- read.csv(pair_groups_fn, header=T, check.names=F)
  pair_groups$desc <- paste0(pair_groups$datatag,'_',
                             pair_groups$sample1,'-',pair_groups$clone1,'_',
                             pair_groups$sample2,'-',pair_groups$clone2)
  observed_pair_groups <- pair_groups[pair_groups$datatag==datatag,'desc']
  rownames(pair_groups) <- pair_groups$desc
  print("Number of pair comparison: ")
  print(length(observed_pair_groups))
  if(length(observed_pair_groups)==0){
    stop("DEBUG: No comparison to do, double check input setting")
  }
  
  
  print(dim(sce))
  
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  sum(mito_genes==TRUE)
  
  ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
  sum(ribo_genes==TRUE)
  sce <- sce[(!mito_genes) & (!ribo_genes), ]
  print("Observed sce: ")
  print(paste0('Removing mito and ribo genes: ',dim(sce)[1],'_',dim(sce)[2]))
  
  
  meta_data <- as.data.frame(colData(sce))
  # summary(as.factor(meta_data$clone))
  
  for(de in observed_pair_groups){
    cl1 <- unlist(strsplit(as.character(pair_groups[de,'clone1']), "_"))
    cl2 <- unlist(strsplit(as.character(pair_groups[de,'clone2']), "_"))
    cells_use_g1 <- rownames(meta_data)[meta_data$clone %in% cl1 & meta_data$sample==pair_groups[de,'sample1']]
    cells_use_g2 <- rownames(meta_data)[meta_data$clone %in% cl2 & meta_data$sample==pair_groups[de,'sample2']]
    
    
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
      # by not testing genes that are very infrequently expressed. Default is 0.1, means 10%
      print("Filtering by pct minimum fraction genes in each population")
      zero_g1 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$clone %in% cl1], "counts") == 0)
      obs_genes1 <- names(zero_g1[zero_g1 <= (1 - min_pct)])
      print(length(obs_genes1))
      
      zero_g2 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$clone %in% cl2], "counts") == 0)
      obs_genes2 <- names(zero_g2[zero_g2 <= (1 - min_pct)])
      print(length(obs_genes2))
      
      sce_tmp <- sce_tmp[intersect(obs_genes1, obs_genes2),]
      
      sce_tmp <- sce_tmp[, sce_tmp$total_features_by_counts > min_features]
      sce_tmp$clone <- ifelse(sce_tmp$clone %in% cl1, paste0("2_",pair_groups[de,'clone1']),paste0("1_",pair_groups[de,'clone2']))
      
      minLogFC <- 0.25  # consider to change this threshold to 0.5, get only genes that are significant difference between 2 groups
      pValueThrs <- 0.05
      de_genes <- edgeR_de(sce_tmp, save_dir_pw)
      
      de_genes <- de_genes %>%
        dplyr::filter(abs(logFC)>minLogFC)
    }
  }  
      
}