suppressMessages({
  # library("sleuth")
  library("dplyr")
  library("tximport")
  library("tximportData")
  library("TxDb.Hsapiens.UCSC.hg19.knownGene")
  library("edgeR")
  library("csaw")
  library(tidyverse)
  library("ggplot2")
})

load_kallisto_data <- function(meta_df){
  tx2gene <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene.csv', header=T) %>% as.data.frame()
  tx2gene$V1 <- NULL
  head(tx2gene)
  dim(tx2gene)
  files <- paste0(input_dir, meta_df$bulk_output_id,'/abundance.tsv')
  for(f in files){
    if(!file.exists(f)){
      print(f)
      stop('File do not exist')
      
    }
  }
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  head(txi$counts)
  
  colnames(txi$counts) <- meta_df$sample_id
  
  # Get genes symbols here
  chrs <- c(as.character(1:22), "X")
  t <- annotables::grch38
  colnames(t)
  annots <- annotables::grch38 %>%
    dplyr::select(ensembl_gene_id = ensgene, gene_symbol=symbol,chr) %>%
    dplyr::filter(chr %in% chrs)
  annots <- annots[!duplicated(annots$ensembl_gene_id),]
  colnames(annots)
  dim(txi$counts)
  txi$counts_values <- as.data.frame(txi$counts)
  
  gn <- sapply(strsplit(rownames(txi$counts_values),'\\.'),function(x){
    return(x[1])
  })
  
  txi$counts_values$ensembl_gene_id <- as.character(gn)
  # head(txi$counts_values)
  txi$counts_values <- txi$counts_values %>% left_join(annots, by='ensembl_gene_id')
  # sum(!is.na(txi$counts_values$gene_symbol))
  # df <- txi$counts_values
  # View(head(df))
  # df2 <- df %>%
  #   filter(gene_symbol %in% grep('MAEL',df$gene_symbol, value=T))
  # df2
  # df$gene_symbol[50:150]
  # sum(grepl('HIST',df$gene_symbol))
  # class(df$gene_symbol)
  return(txi)
  
}

# Normalizing data using edgeR
# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normalize_by_size_factor <- function(txi, groups_use){
  cts <- txi$counts
  normMat <- txi$length
  
  # Obtaining per-observation scaling factors for length, adjusted to avoid
  # changing the magnitude of the counts.
  normMat <- normMat/exp(rowMeans(log(normMat)))
  normCts <- cts/normMat
  
  # Computing effective library sizes from scaled counts, to account for
  # composition biases between samples.
  # library(edgeR)
  eff.lib <- edgeR::calcNormFactors(normCts) * colSums(normCts)
  
  # Combining effective library sizes with the length factors, and calculating
  # offsets for a log-link GLM.
  normMat <- sweep(normMat, 2, eff.lib, "*")
  normMat <- log(normMat)
  
  # Creating a DGEList object for use in edgeR.
  # class(cts)
  # dim(cts)
  # head(cts)
  dge <- edgeR::DGEList(cts, group = groups_use)
  dge <- edgeR::scaleOffset(dge, normMat)
  # filtering
  print(dim(dge))
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, ]
  print(dim(dge))
  class(dge)
  cpms <- edgeR::cpm(dge, offset = dge$offset, log = FALSE)
  
  # y is now ready for estimate dispersion functions see edgeR User's Guide
  return(list(dge=dge,cpm=cpms))
}

load_dge_counts <- function(cts, groups_use){
  # Creating a DGEList object for use in edgeR.
  # class(cts)
  # dim(cts)
  # head(cts)
  dge <- edgeR::DGEList(cts, group = groups_use)
  
  # dge <- edgeR::calcNormFactors(cts)
  # dge$genes <- meta_genes
  
  # normMat <- edgeR::cpm(cts, log=TRUE) 
  # dge <- edgeR::scaleOffset(dge, normMat)
  # filtering
  print(dim(dge))
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, ]
  print(dim(dge))
  print(class(dge))
  dge$samples$lib.size <- colSums(dge$counts)
  dge <- calcNormFactors(dge, method = "TMM")
  
  # lcpm_dge <- cpm(dge, log=TRUE)
  # dim(lcpm_dge)
  # View(head(lcpm_dge))
  # y is now ready for estimate dispersion functions see edgeR User's Guide
  return(dge)
}


edgeR_de <- function(dge, meta_df, save_dir, tag=NULL){
  
  # print("Filtering data...")
  # sce_de <- sce_de[, sce_de$total_features_by_counts > 1500]
  # # rs <- rowSums(as.matrix(counts(sce_de)))
  # # qplot(rs, log='x') + geom_vline(xintercept = 100)
  # 
  # sce_de <- sce_de[rowSums(as.matrix(counts(sce_de))) > 100, ]
  # print(dim(sce_de))
  # mycounts <- as.matrix(counts(sce_de))    # zeros are good
  # print("Create DGE edgeR object...")
  
  # dge <- DGEList(counts=mycounts, group=sce_de$clone)
  # dge <- DGEList(counts=mycounts, group=sce_de$treatmentSt)
  
  print("DE Analysis...")
  # which vs which ???
  # design <- model.matrix(~ clone, data = colData(sce_de))
  design <- model.matrix(~condition, data = meta_df)
  if(is.null(rownames(design))){
    rownames(design) <- meta_df$sample
  }
  
  # 
  # design <- model.matrix(~ treatmentSt, data = colData(sce_de))
  # This describes the edgeR user manual
  # http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  dge <- calcNormFactors(dge, method = "TMM")
  # print(dge$samples)
  dge <- edgeR::estimateDisp(dge, design = design)
  # dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
  
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit)
  tt <- edgeR::topTags(qlf, n = Inf)
  
  print("Generate output...")
  tt <- as.data.frame(tt) %>% 
    tibble::rownames_to_column("gene_id")
  # tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$ID
  # tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$Symbol
  # View(head(tt))
  print(dim(tt))
  # tt$ens_gene_id <- get_ens_gene(tt$gene_id)
  tt$gene_symbol <- get_gene_symbol(tt$gene_id)
  print(length(unique(tt$gene_symbol)))
  # print(summary(tt$logFC))
  # print(summary(tt$PValue))
  # Saving the logFC file
  # tt <- tt[abs(tt$logFC)>0.25,] # tt$FDR<0.01 & tt$PValue<0.05 
  print(dim(tt))
  if(is.null(tag)){
    tag <- paste(unique(meta_df$condition),collapse='_')
  }
  write.csv(tt, file=paste0(save_dir, 'edgeR_significant_genes_',tag,'.csv'), row.names=FALSE, quote=FALSE)
  print("With FDR<0.01 and PValue<0.05 and abs(tt$logFC)>0.25, number of significant genes is: ")
  print(dim(tt))
  # print("With threshold logFC>0.25, number of significant genes is: ")
  # print(nrow(tt[abs(tt$logFC)>0.25,]))
  return(tt)
}

get_gene_symbol <- function(ens_gene_ids){
  
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  genes_symb_df <- read.csv(paste0(base_dir,'biodatabase/meta_genes.csv'), check.names = F, stringsAsFactors = F)
  # dim(genes_symb_df)
  rownames(genes_symb_df) <- genes_symb_df$gene_ens
  return(genes_symb_df[ens_gene_ids,'gene_symb'])
}

edgeR_de_v2 <- function(dge, meta_df, save_dir, save_dge=T, tag=NULL){
  
  # print("Filtering data...")
  # sce_de <- sce_de[, sce_de$total_features_by_counts > 1500]
  # # rs <- rowSums(as.matrix(counts(sce_de)))
  # # qplot(rs, log='x') + geom_vline(xintercept = 100)
  # 
  # sce_de <- sce_de[rowSums(as.matrix(counts(sce_de))) > 100, ]
  # print(dim(sce_de))
  # mycounts <- as.matrix(counts(sce_de))    # zeros are good
  # print("Create DGE edgeR object...")
  
  # dge <- DGEList(counts=mycounts, group=sce_de$clone)
  # dge <- DGEList(counts=mycounts, group=sce_de$treatmentSt)
  
  print("DE Analysis...")
  # which vs which ???
  # design <- model.matrix(~ clone, data = colData(sce_de))
  design <- model.matrix(~condition, data = meta_df)
  if(is.null(rownames(design))){
    rownames(design) <- meta_df$sample
  }
  # 
  # design <- model.matrix(~ treatmentSt, data = colData(sce_de))
  # This describes the edgeR user manual
  # http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- edgeR::estimateDisp(dge, design = design)
  # dge$common.dispersion
  # dge <- edgeR::estimateCommonDisp(dge, design = design)
  # dge <- edgeR::estimateTagwiseDisp(dge)
  # dge <- edgeR::estimateGLMCommonDisp(dge, design = design)
  # de.com <- exactTest(dge, dispersion=predefined_disp)
  
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit)
  tt <- edgeR::topTags(qlf, n = Inf)
  
  print(summary(decideTests(qlf)))
  png(paste0(save_dir,tag,"_DE_genes.png"), height = 2*350, width=2*500, res = 2*72)
  plotMD(qlf)
  abline(h=c(-1, 1), col="blue")
  dev.off()
  
  
  print("Generate output...")
  tt <- as.data.frame(tt) %>% 
    tibble::rownames_to_column("gene_symbol")
  # tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$ID
  # tt$gene_symbol <- rowData(sce_de[tt$gene_id,])$Symbol
  # View(head(tt))
  print(dim(tt))
  # tt$ens_gene_id <- get_ens_gene(tt$gene_id)
  tt$gene_symbol <- get_gene_symbol(tt$gene_id)
  # print(length(unique(tt$gene_symbol)))
  print(summary(tt$logFC))
  print(summary(tt$PValue))
  # Saving the logFC file
  tt <- tt[abs(tt$logFC)>1,] # tt$FDR<0.01 & tt$PValue<0.05 
  print(head(tt))
  if(is.null(tag)){
    tag <- paste(unique(meta_df$condition),collapse='_')
  }
  write.csv(tt, file=paste0(save_dir, 'edgeR_significant_genes_',tag,'.csv'), row.names=FALSE, quote=FALSE)
  # print("With threshold logFC>0.25, number of significant genes is: ")
  # print(nrow(tt[abs(tt$logFC)>0.25,]))
  if(save_dge){
    saveRDS(dge,paste0(save_dir, 'edgeR_dge',tag,'.rds'))  
  }
  
  return(tt)
}

get_gene_symbol_from_transcript <- function(target_ids){
  genes <- mclapply(strsplit(target_ids,'\\|'), function(x){
    return(x[6])
  }, mc.cores =3)
  return(as.character(genes))
}
get_ens_gene_from_transcript <- function(target_ids){
  genes <- mclapply(strsplit(target_ids,'\\|'), function(x){
    return(x[2])
  }, mc.cores =3)
  genes <- as.character(genes)
  # genes <- genes[grepl('^ENSG',genes)]
  return(genes)
}


viz_heatmap_rawdata <- function(counts_df, meta_df, save_dir, datatag){
  # library(ComplexHeatmap)
  lcpm_counts <- edgeR::cpm(counts_df, log=T)
  samples_use <- colnames(lcpm_counts)
  predefined_cols <- c("#666666","#66A61E")
  names(predefined_cols) <- c('2_Rx','1_UnRx')
  # cols_use <- predefined_cols[meta_df[gsub(paste0(datatag,'_'),samples_use),'condition']]
  print(samples_use)
  cols_use <- predefined_cols[meta_df[samples_use,'condition']]
  # cols_use=c("#66A61E", "#66A61E", "#666666") # green and grey
  names(cols_use) <- samples_use
  print(cols_use)
  left_anno = rowAnnotation(Sample = factor(samples_use),
                            col = list(Sample=cols_use))
  p <- ComplexHeatmap::Heatmap(as.matrix(t(lcpm_counts)), na_col = "white",
                               show_column_names=F,
                               show_row_names = T,
                               cluster_rows=T,
                               cluster_columns=F,
                               name = "LogCPMExp", 
                               # row_order = sort(rownames(test)),
                               # row_split= samples_use,
                               row_title_rot = 0,
                               row_gap = unit(2, "mm"),
                               # column_split = genes_type$gt,
                               column_title = paste0("Raw Data Clustering ",datatag),
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 10),
                               row_names_gp = grid::gpar(fontsize = 10),
                               show_heatmap_legend = T,
                               # top_annotation=top_anno,
                               left_annotation = left_anno,
                               # cell_fun = cell_func,
                               row_dend_reorder=F)
  
  png(paste0(save_dir,"hm_rawdata.png"), height = 2*280, width=2*700, res = 2*72)
  print(p)
  dev.off()
  
}

viz_heatmap_DE_genes <- function(counts_df, meta_df, save_dir, datatag, de_genes_fn=NULL){
  # library(ComplexHeatmap)
  if(is.null(de_genes_fn)){
    de_genes_fn <- paste0(save_dir, 'edgeR_significant_genes_1_UnRx_2_Rx.csv')
    de_genes <- data.table::fread(de_genes_fn) %>% as.data.frame()
    print(dim(de_genes))
    
  }
  counts_df <- counts_df[de_genes$gene_symbol,]
  print(dim(counts_df))
  lcpm_counts <- edgeR::cpm(counts_df, log=T)
  samples_use <- colnames(lcpm_counts)
  predefined_cols <- c("#666666","#66A61E")
  names(predefined_cols) <- c('2_Rx','1_UnRx')
  cols_use <- predefined_cols[meta_df[samples_use,'condition']]
  # cols_use=c("#66A61E", "#66A61E", "#666666") # green and grey
  names(cols_use) <- samples_use
  cols_use
  
  left_anno = rowAnnotation(Sample = factor(samples_use),
                            col = list(Sample=cols_use))
  genes_type <- data.frame(gene_symbol=rownames(lcpm_counts), stringsAsFactors=F)
  genes_type <- genes_type %>% left_join(de_genes, by=c('gene_symbol'))
  # dim(genes_type)
  genes_type$gt <- ifelse(genes_type$logFC>=1 & genes_type$PValue<=0.05,'up_regulated',
                          ifelse(genes_type$logFC<=-1 & genes_type$PValue<=0.05,'down_regulated','NotSig'))
  
  genes_type$gt <- factor(genes_type$gt, levels = c('up_regulated','down_regulated','NotSig'))
  top_anno = HeatmapAnnotation(Gene_Type = factor(genes_type$gt),
                               col = list(Gene_Type = c(up_regulated = "#0000FF", 
                                                        down_regulated = "#D95F02",
                                                        NotSig='#AAB2B1')))
  lcpm_counts <- lcpm_counts[genes_type$gene_symbol,]
  
  p <- ComplexHeatmap::Heatmap(as.matrix(t(lcpm_counts)), na_col = "white",
                               show_column_names=F,
                               show_row_names = T,
                               cluster_rows=T,
                               cluster_columns=F,
                               name = "LogCPMExp", 
                               # row_order = sort(rownames(test)),
                               # row_split= samples_use,
                               row_title_rot = 0,
                               row_gap = unit(2, "mm"),
                               column_split = genes_type$gt,
                               column_title = paste0("Raw Data Clustering ",datatag),
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 10),
                               row_names_gp = grid::gpar(fontsize = 10),
                               show_heatmap_legend = T,
                               top_annotation=top_anno,
                               left_annotation = left_anno,
                               # cell_fun = cell_func,
                               row_dend_reorder=F)
  
  png(paste0(save_dir,"hm_DE_genes.png"), height = 2*280, width=2*700, res = 2*72)
  print(p)
  dev.off()
  
}

# suffix='-R.counts.genes'
#suffix='_R2.counts.genes' in some cases Farhia's project
load_counts_data <- function(samples, input_dir, datatag, suffix='-R.counts.genes'){ #
  counts <- list()
  for(sample in samples){
    sample_fn <- paste0(input_dir, sample,'/counts/',sample,suffix)
    suffix2='_R.counts.genes'
    sample_fn2 <- paste0(input_dir, sample,'/counts/',sample,suffix2)
    if(file.exists(sample_fn)){
      scounts <- data.table::fread(sample_fn, header=F,sep='\t')  
      print(sample_fn)
    }else if(file.exists(sample_fn2)){
      scounts <- data.table::fread(sample_fn2, header=F,sep='\t')  
      print(sample_fn2)
    }else{
      stop('Double check input data')
    }
    
    
    dim(scounts)
    colnames(scounts) <- c('gene_symbol',sample)
    # head(scounts)
    print(colnames(scounts))
    
    # cts <- scounts
    
    scounts <- scounts %>% remove_rownames %>% column_to_rownames(var="gene_symbol")
    counts[[sample]] <- scounts
  }
  counts_df <- dplyr::bind_cols(counts)
  print(dim(counts_df))
  # View(head(counts_df))  
  return(counts_df)
}

run_edgeR <- function(datatag, meta_df, input_dir, base_dir){
  # library(stringr)
  groups_use <- c('Rx','UnRx')
  meta_df$condition <- ifelse(grepl(groups_use[2], meta_df$treatmentSt),paste0("1_",groups_use[2]),
                              paste0("2_",groups_use[1]))
  print(meta_df)
  meta_df$condition <- factor(meta_df$condition,levels = c("1_UnRx","2_Rx"))
  # meta_df$sample <- c('3617','3616','3664')
  meta_df$sample <- stringr::str_sub(meta_df$sample, str_length(meta_df$sample)-3, str_length(meta_df$sample))
  print(meta_df$sample)
  rownames(meta_df) <- meta_df$sample
  counts_df <- load_counts_data(meta_df$sample, input_dir, datatag)
  print(colnames(counts_df))
  # rm(dge)
  dge <- load_dge_counts(counts_df, meta_df$condition)
  # predefined_disp <- 0.01342844
  save_dir_de <- paste0(base_dir,paste(meta_df$sample,collapse = '_'),'/')
  dir.create(save_dir_de, recursive = T)
  de_genes <- edgeR_de_v2(dge, meta_df, save_dir_de, tag=NULL)
  res <- list(counts_df=counts_df,de_genes=de_genes, 
              save_dir=save_dir_de, datatag=datatag,
              meta_df=meta_df)
  return(res)
  
}


load_txt2genes <- function(){
  # txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  # k <- keys(txdb, keytype = "TXNAME")
  # tx2gene <- select(txdb, k, "GENEID", "TXNAME")
  # tx2gene <- select(txdb, k, "GENEID", "TXNAME")
}

getLogFC <- function(gene_df, save_dir, tag=''){
  tag <- paste(colnames(gene_df), collapse = '_')
  gene_df$logFC <- round(log2(gene_df[,1]) - log2(gene_df[,2]), 2)
  gene_df$bulk_ens_gene_id <- rownames(gene_df)
  gene_df$ens_gene_id <- sapply(gene_df$bulk_ens_gene_id,function(x){
    return(as.character(strsplit(as.character(x), "[.]")[[1]][1]))
  })
  ref <- annotables::grch37 %>%
    # dplyr::select(gene_symbol = symbol, chr, start, end) %>%
    dplyr::select(ens_gene_id = ensgene, gene_symbol=symbol, chr)
  gene_df <- gene_df %>% left_join(ref, by=c('ens_gene_id'))
  print(sum(is.na(gene_df$gene_symbol)))
  # gene_df$gene_symbol <- get_gene_symbol(gene_df$ens_gene_id)
  write.csv(gene_df, file=paste0(save_dir, 'log2FC_',tag,'.csv'), row.names=FALSE, quote=FALSE)
}
## Reading txt files
# files <- file.path(dir, "kallisto", samples$run, "abundance.tsv.gz")
# names(files) <- paste0("sample", 1:6)
# txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
# head(txi.kallisto.tsv$counts)

viz_cross_correlation <- function(de_genes_fn, de_genes_fn_bulk, save_dir, 
                                  de_desc_sc='',de_desc_bulk='',
                                  tag='',datatag=''){
  de_genes <- data.table::fread(de_genes_fn) %>% as.data.frame()
  print(dim(de_genes))
  # head(de_genes)
  de_genes <- de_genes %>%
    dplyr::filter(abs(avg_log2FC)>0.25 & p_val_adj <= 0.05)%>%
    dplyr::rename(singlecell_logFC=avg_log2FC, gene_symbol=gene_symb)
  
  # edgeR
  # de_genes <- de_genes %>%
  #   dplyr::filter(abs(logFC)>0.25 & FDR<0.01 & PValue < 0.05)%>%
  #   dplyr::rename(singlecell_logFC=logFC)
  print(dim(de_genes))
  
  de_genes_bulk <- data.table::fread(de_genes_fn_bulk) %>% as.data.frame()
  print(dim(de_genes_bulk))
  head(de_genes_bulk)
  
  de_genes_bulk$gene_symbol
  # View(de_genes_bulk)
  de_genes_bulk <- de_genes_bulk %>%
    dplyr::filter(abs(logFC)>0.25) %>%
    # dplyr::select(gene_symbol,logFC)%>%
    dplyr::rename(bulk_logFC=logFC)
  # de_genes_bulk$bulk_logFC[1]
  de_genes <- de_genes %>% inner_join(de_genes_bulk, by=c('gene_symbol'))
  dim(de_genes)
  
  # de_genes$bulk_logFC <- ifelse(!is.na(de_genes$bulk_logFC),de_genes$bulk_logFC, 0)
  save_dir <- paste0(save_dir, tag, '_')
  
  cr <- cor(de_genes$singlecell_logFC, de_genes$bulk_logFC, method = "spearman")
  pltttitle <- paste0(datatag,': SingleCell-Bulk \n Correlation: ',round(cr,2))
  de_genes$DE_gene <- ifelse(de_genes$singlecell_logFC>0,'Up-regulated','Down-regulated')
  library(ggplot2)
  # de_desc_sc <- 'Res(X8:U) - Sen(X9:J)'
  # de_desc_bulk <- 'Res(X8:Rx) - Sen(X8,X6:UnRx)'
  # dim(de_genes)
  p <- ggplot(de_genes, aes(singlecell_logFC, bulk_logFC, color=DE_gene, shape=DE_gene)) + 
    geom_point(size=2.3, alpha=1) +
    scale_shape_manual(values=c('Up-regulated'=1, 'Down-regulated'=2)) +
    scale_color_manual(values=c('Up-regulated'='#E69F00', 'Down-regulated'='#56B4E9')) + 
    geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.3) + 
    geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) 
  p <- p + labs(x=paste0('SingleCell log2FC ',de_desc_sc),y=paste0('Bulk log2FC ',de_desc_bulk), title = pltttitle) + 
    theme(plot.title = element_text(color="black", size=13, hjust = 0.5, face = "bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.y = element_text(color="black", size=12),
          axis.text.x = element_text(color="black", size=12),
          axis.title = element_text(color="black", size=12),
          legend.title = element_text(color="black", size=12), 
          legend.text = element_text(color="black", size=12))
  # p
  # saveRDS(p, paste0(output_dir,"edgeR_corr_eval_",de_desc,".rds"))
  png(paste0(save_dir,"singlecell_bulk_corr_eval.png"), height = 2*350, width=2*500,res = 2*72)
  print(p)
  dev.off()
  return(p)
  
}