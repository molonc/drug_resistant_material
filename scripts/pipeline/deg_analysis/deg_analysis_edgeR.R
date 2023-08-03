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




# datatag: 'SA535_CX', 'SA535_CY','SA1035', 'SA609'
edgeR_DE_analysis_by_clones_samples <- function(pair_groups_fn, 
                                                sce_fn, output_file, datatag='SA', 
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
  # sce <- sce[(!mito_genes) & (!ribo_genes), ]
  sce <- sce[!mito_genes, ] # regress out mito genes
  print("Observed sce: ")
  print(paste0('Removing mito genes, new sce file: ',dim(sce)[1],'_',dim(sce)[2]))
  
  
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
      # by not testing genes that are very infrequently expressed. Default is 0.1
      print("Filtering by pct minimum fraction genes in each population")
      zero_g1 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$clone %in% cl1], "counts") == 0)
      obs_genes1 <- names(zero_g1[zero_g1 <= (1 - min_pct)])
      print(length(obs_genes1))
      zero_g2 <- DelayedArray::rowMeans(assay(sce_tmp[,sce_tmp$clone %in% cl2], "counts") == 0)
      obs_genes2 <- names(zero_g2[zero_g2 <= (1 - min_pct)])
      print(length(obs_genes2))
      sce_tmp <- sce_tmp[intersect(obs_genes1, obs_genes2),]
      cat("\n Filtering by pct minimum fraction genes:\n ",file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      cat(dim(sce_tmp),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      sce_tmp <- sce_tmp[, sce_tmp$total_features_by_counts > min_features]
      cat("\n Filtering by total_features_by_counts greater than 1500:\n ",file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      cat(dim(sce_tmp),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      sce_tmp$clone <- ifelse(sce_tmp$clone %in% cl1, paste0("2_",pair_groups[de,'clone1']),paste0("1_",pair_groups[de,'clone2']))
      cat(summary(as.factor(sce_tmp$clone)),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      
      minLogFC <- 0.25
      pValueThrs <- 0.05
      de_genes <- edgeR_de(sce_tmp, save_dir_pw)
      
      de_genes <- de_genes[abs(de_genes$logFC)>minLogFC,]
      cat(nrow(de_genes),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      # Plot DE genes
      plttitle <- paste0(pair_groups[de,'datatag'],":  ",pair_groups[de,'clone1']," versus ",pair_groups[de,'clone2'])
      nbtopup <- 30
      nbtopdown <- 30
      markers_ls_upreg <- de_genes[de_genes$logFC>minLogFC,]
      markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$logFC,decreasing = T),] 
      cat(paste0('\n ',de), file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      cat(paste0('\n Nb resistant genes: ', nrow(markers_ls_upreg)),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      pair_groups[de,'resistant_genes'] <- nrow(markers_ls_upreg)
     
      # dim(markers_ls_upreg)
      if(nrow(markers_ls_upreg) < nbtopup){
        nbtopup <- nrow(markers_ls_upreg)
      }
      markers_ls_upreg <- markers_ls_upreg[1:nbtopup,]
      
      # Get top down-regulated genes
      markers_ls_downreg <- de_genes[de_genes$logFC<(-minLogFC),]
      markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$logFC,decreasing = F),] 
      cat(paste0('\n Nb sensitive genes: ', nrow(markers_ls_downreg)),file = paste0(save_dir_pw,"de_analysis_log.txt"), append = TRUE)
      pair_groups[de,'sensitive_genes'] <- nrow(markers_ls_downreg)
      # dim(markers_ls_downreg)
      if(nrow(markers_ls_downreg) < nbtopdown){
        nbtopdown <- nrow(markers_ls_downreg)
      }
      markers_ls_downreg <- markers_ls_downreg[1:nbtopdown,]
      topGenes <- c(as.character(markers_ls_upreg$gene_symbol),as.character(markers_ls_downreg$gene_symbol))
      # rownames(markers_ls_tmp) <- markers_ls_tmp$gene_symb
      plot_DE_genes_edgeR(de_genes, topGenes, nrow(de_genes), 
                          0.01, 0.25, plttitle, save_dir_pw, legendVisible=F)
    }  
  }
  
  
  write.csv(pair_groups,file=pair_groups_fn, quote=F, row.names=F)
}

datatag <- 'SA609'
# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
# results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6'
# output_file <- paste0(input_dir,'/SA609_rna/deg_analysis/SA609_deg_pathway.rds')
pair_groups_fn <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/deg_analysis/pair_groups_SA1035_SA609.csv')
# sce_fn <- paste0(results_10x_dir,'/total_sce_clones.rds')
# edgeR_DE_analysis_by_clones_samples(pair_groups_fn, sce_fn, output_file, datatag, 1500, 0.1)
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results','/',datatag,'_rna/deg_analysis/SA609-v6/')
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/added_segments/clonealign/whole_data/'
cnv_mat <- readRDS(paste0(dlp_dir,datatag,'_cnv_mat.rds'))
cnv_mat <- cnv_mat[,obs_clones]
dim(cnv_mat)
cnv_mat$ensembl_gene_id <- rownames(cnv_mat)
cnv_mat <- cnv_mat %>%
  pivot_longer(!ensembl_gene_id, names_to = "cluster", values_to = "cnv")

df_cnv <- cnv_mat
dim(df_cnv)
cnv_mat <- read.csv(paste0(input_dir,'SA609_UTTT_R_UUUU_H/cnv_mat_total.csv'), check.names = F, stringsAsFactors = F)
head(cnv_mat_clones)
get_gene_type_edgeR(pair_groups_fn, cnv_mat, 
                    datatag, input_dir, 
                    cancer_ref_genes_fn=NULL, outlier_FC_thrs=4.25)
p <- ggplot(stat1, aes(fill=gene_type, y=pct_genes, x=desc)) + 
  geom_bar(position="fill", stat="identity")+
  facet_grid(. ~ datatag, scales="free", space='free')
p

obs_clones=c('R','H')
View(head(cnv_mat))
input_deg_fn <- paste0(input_dir,'SA609_UTTT_R_UUUU_H/signif_genes.csv')



input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
datatag <- 'SA1035'
results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA1035-v6'
output_file <- paste0(input_dir,'/',datatag,'_rna/deg_analysis/',datatag,'_deg_pathway.rds')
pair_groups_fn <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/deg_analysis/pair_groups_SA1035_SA609.csv')
sce_fn <- paste0(results_10x_dir,'/total_sce_clones.rds')
# edgeR_DE_analysis_by_clones_samples(pair_groups_fn, sce_fn, output_file, datatag, 1500, 0.1)

#Get intrans, incis genes 
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results','/',datatag,'_rna/deg_analysis/')
dlp_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/clonealign/whole_data/'
cnv_mat <- readRDS(paste0(dlp_dir,datatag,'_cnv_mat.rds'))
get_gene_type_edgeR(pair_groups_fn, cnv_mat, 
                    datatag, input_dir, 
                    cancer_ref_genes_fn=NULL, outlier_FC_thrs=4.25)
p <- ggplot(stat1, aes(fill=gene_type, y=pct_genes, x=desc)) + 
  geom_bar(position="fill", stat="identity")+
  facet_grid(. ~ datatag, scales="free", space='free')


input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
datatag <- 'SA535_CY'
results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA535-v6'
output_file <- paste0(input_dir,'/SA535_total_rna_v2/deg_analysis_',datatag,'/',datatag,'_deg_pathway.rds')
pair_groups_fn <- paste0(input_dir,'/SA609_rna/deg_analysis/pair_groups_SA535.csv')
sce_fn <- paste0(results_10x_dir,'/total_sce_clones.rds')
edgeR_DE_analysis_by_clones_samples(pair_groups_fn, sce_fn, output_file, datatag, 1500, 0.1)


input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
datatag <- 'SA535_CX'
results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA535-v6'
output_file <- paste0(input_dir,'/SA535_total_rna_v2/deg_analysis_',datatag,'/',datatag,'_deg_pathway.rds')
pair_groups_fn <- paste0(input_dir,'/SA609_rna/deg_analysis/pair_groups_SA535.csv')
sce_fn <- paste0(results_10x_dir,'/total_sce_clones.rds')
edgeR_DE_analysis_by_clones_samples(pair_groups_fn, sce_fn, output_file, datatag, 1500, 0.1)



# pair_groups_1 <- paste0(input_dir,'/SA609_rna/deg_analysis/pair_groups_SA535.csv')
# pair_groups1 <- read.csv(pair_groups_1, header=T, check.names=F)
# pair_groups_2 <- paste0(input_dir,'/SA609_rna/deg_analysis/pair_groups_SA1035_SA609.csv')
# pair_groups2 <- read.csv(pair_groups_2, header=T, check.names=F)
# pair_groups <- rbind(pair_groups1,pair_groups2)
# pair_groups_fn <- paste0(input_dir,'/rnaseq_v6/DE_analysis.csv')
# dim(pair_groups)
# write.csv(pair_groups,file=pair_groups_fn, quote=F, row.names=F)
