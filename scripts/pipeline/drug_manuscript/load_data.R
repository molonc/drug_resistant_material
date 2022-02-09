suppressPackageStartupMessages({
  # library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/encoder_trajectory/'
meta_df <- data.table::fread(paste0(output_dir,'grouping_SA535.csv.gz')) %>% as.data.frame()
dim(meta_df)
table(meta_df$clone, meta_df$treatmentSt)
colnames(metadata2)
metadata3 <- meta_df %>%
  dplyr::filter(grepl('A', clone) & treatmentSt!='UUUU' & !grepl('TU$', treatmentSt))
dim(metadata3)  
metadata3 <- meta_df %>%
  dplyr::filter(grepl('G', clone) & !grepl('TU$', treatmentSt))
dim(metadata3)  
table(metadata3$clone, metadata3$treatmentSt)
table(metadata3$mouse_id, metadata3$treatmentSt)

# Option 1: UT - SA535X6XB03101, UTT-SA535X7XB03304, UTTT - SA535X8XB03431, UTTTT - SA535X9XB03617 (main clone A - resistance to drug)
# analysis: drug treatment across time
# Option 2: UTT-SA535X7XB03304, UTTT - SA535X8XB03431, UTTTT - SA535X9XB03617, UUUU-SA535X8XB03664 or UUUUU-SA535X9XB03776, X8 seems better quality(main clone A - resistance to drug, G - sensitive to drug, late time points)
# analysis: drug treatment across time + validate DE analysis resistance A vs sensitive G
metadata2 <- metadata2 %>%
  select(-PDX)
metadata <- data.table::fread(paste0(base_dir,'out/SA535_sctransform_umap.csv')) %>% as.data.frame()
metadata2 <- data.table::fread(paste0(base_dir,'out/metasample_SA535_cisplatin_untreated.csv')) %>% as.data.frame()
meta_df <- meta_df %>% left_join(metadata, by=c('cell_id','clone'))
meta_df <- meta_df %>% left_join(metadata2, by=c('mouse_id','library_id'))
dim(meta_df)
colnames(meta_df)
data.table::fwrite(meta_df, paste0(output_dir,'grouping_SA535_v2.csv.gz'))
  
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v7/SA535_v7/'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/monocle_v2/'
dir.create(output_dir)
datatag <- 'SA535'
input_dir <- paste0(base_dir,'clonealign/')
sce_dir <- paste0(base_dir,'sce/')
output_dir <- base_dir
load_metadata <- function(sce_dir, input_dir){
  fns <- list.files(input_dir)
  fns <- fns[grepl('.csv',fns)]
  clonealign_df <- tibble::tibble()
  for(f in fns){
    tmp <- data.table::fread(paste0(input_dir,f))
    clonealign_df <- dplyr::bind_rows(clonealign_df, tmp)
  }
  dim(clonealign_df)
  # unique(clonealign_df$clone)
  clonealign_df$library_id <- gsub('(.cache/)','',clonealign_df$Sample)
  clonealign_df$library_id <- gsub('(/filtered_feature_bc_matrix)','',clonealign_df$library_id)
  colnames(clonealign_df)
  View(head(clonealign_df))
  clonealign_df <- clonealign_df %>%
    dplyr::select(-Sample)%>%
    dplyr::rename(mouse_id=id)
  data.table::fwrite(clonealign_df, paste0(input_dir,'clonealign_fitness_cisplatin.csv'))
  clonealign_df <- as.data.frame(clonealign_df)
  rownames(clonealign_df) <- paste0(clonealign_df$library_id,'_',clonealign_df$Barcode)
  # clonealign_df$cell_id[1]
  # [1]
  sce_combine <- load_data(clonealign_df, sce_dir,tag='SA535_cisplatin')
  
  
  rowData(sce_combine)$Symbol[1]
  
  # st <- DelayedArray::rowSums(assay(sce_combine, 'counts') != 0)
  dim(sce_combine)
  datatag <- 'SA535'
  norm_sce <- normalize_SCTransform(sce_combine, base_dir, datatag, return_data=T, output_fn=NULL)
  sce <- sce_combine
  umap_df <- umap_df %>% inner_join(clonealign_df, by=c("cell_id"))
  pca_df <- pca_df %>% inner_join(clonealign_df, by=c("cell_id"))
  dim(umap_df)
  write.csv(pca_df, paste0(output_dir, datatag, "_sctransform_pca.csv"), row.names = F, quote = F)
  write.csv(umap_df, paste0(output_dir, datatag, "_sctransform_umap.csv"), row.names = F, quote = F)
  
  var_genes <- srt@assays[['RNA']]@var.features
  var_genes[1]
  write.csv(data.frame(var_genes), paste0(output_dir, datatag, "_var_genes.csv"), row.names = F, quote = F)
  sce <- readRDS(paste0(base_dir,'sce/combined_SA535_cisplatin.rds'))
  dim(sce)
  cut_off_overall <- 0.025
  if(cut_off_overall > 0){
    zero_cbind <- DelayedArray::rowMeans(assay(sce, 'counts') == 0)
    sce <- sce[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
  }
  dim(sce)
  saveRDS(sce, paste0(base_dir,'sce/combined_SA535_cisplatin_filtered_zeros.rds'))
  norm_sce <- sce_normalize_size_factors(sce, min_size=300, rlog=FALSE, exprs="counts")
  dim(norm_sce)
  assayNames(norm_sce)
  genes_use <- rownames(norm_sce)[rowData(norm_sce)$Symbol %in% var_genes]
  length(genes_use)
  genes_use[1]
  saveRDS(norm_sce, paste0(base_dir,'sce/scran normcounts_SA535_cisplatin.rds'))
  scran_df <- as.data.frame(normcounts(norm_sce))
  scran_df <- scran_df[var_genes,]
  
  dim(scran_df)
  View(scran_df[1:3,1:3])
  data.table::fwrite(scran_df, paste0(base_dir,'encoder_trajectory/scran_normcounts_SA535.csv.gz'))
  unique(srt$seurat_clusters)
  dim(clonealign_df)
  grouping <- data.frame(cell_id=colnames(scran_df))
  grouping$clone <- 'unassigned'
  cells_use <- intersect(grouping$cell_id, clonealign_df$cell_id)  
  cells_use <- intersect(grouping$cell_id,cells_use)
  rownames(grouping) <- grouping$cell_id  
  grouping[cells_use,'clone'] <- clonealign_df[cells_use,'clone']
  summary(as.factor(grouping$clone))
  cells_use <- clonealign_df %>%
    dplyr::filter(clone !='unassigned')%>%
    dplyr::pull(cell_id)
  scran_df <- scran_df[,cells_use]
  dim(scran_df)
  rownames(scran_df)[1]
  grouping <- grouping[colnames(scran_df),]
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/encoder_trajectory/'
  data.table::fwrite(scran_df, paste0(save_dir,'scran_normcounts_SA535.csv.gz'), row.names=T)
  data.table::fwrite(grouping, paste0(base_dir,'encoder_trajectory/grouping_SA535.csv.gz'))
  ?data.table::fwrite
  genes_ids <- data.frame(gene_id = rownames(norm_df))
  data.table::fwrite(genes_ids, paste0(save_dir,'genes_ids_3000hvg.csv.gz'))
  
  
  norm_sce <- readRDS(paste0(save_dir,'SA535_sctransform_normalized.rds'))
  rownames(norm_sce) <- rowData(norm_sce)$Symbol
  norm_sce <- norm_sce[rowData(norm_sce)$Symbol %in% var_genes,cells_use]
  norm_df <- logcounts(norm_sce)
  View(as.matrix(norm_df[1:5,1:5]))
  data.table::fwrite(as.matrix(norm_df), paste0(save_dir,'sctransform_SA535.csv.gz'), row.names=T)
  rownames(norm_df)[778]
  
}



base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/SA609_same_chip/'

fns <- list.files(base_dir)
fns <- fns[grepl('*.csv',fns)]
fns <- fns[grepl('XB',fns)]

clones_df <- tibble::tibble()
for(f in fns){
  tmp <- data.table::fread(paste0(base_dir, f))
  clones_df <- dplyr::bind_rows(clones_df, tmp)
}
dim(clones_df)
dim(sce_combine)
sum(colnames(sce_combine) %in% clones_df$cell_id)
clones_df$cell_id[1]
sce_combine$id[1]
colnames(sce_combine) <- paste0(sce_combine$id,'_',sce_combine$Barcode)
sum(t %in% clones_df$cell_id)
colData(sce_combine)[clones_df$cell_id,'clone'] <- clones_df$clone
data.table::fwrite(clones_df, paste0(base_dir, 'SA609_10x_samechip_clonelabels.csv'))
umap_df$cell_id[1]
umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
umap_df$scell_id[1]
colnames(clones_df)
clones_df <- clones_df %>%
  dplyr::select(clone, cell_id)
umap_df <- umap_df %>% left_join(clones_df, by=c('scell_id'='cell_id'))
data.table::fwrite(umap_df, paste0(input_dir,datatag,'_norm_umap.csv'))

fns <- list.files(base_dir)
fns <- fns[grepl('.rds',fns)]
fns <- gsub('.rds','',fns)
meta_cells <- data.frame(mouse_id=fns)
meta_cells$timepoint <- stringr::str_sub(meta_cells$mouse_id, 6, 7)
data.table::fwrite(meta_cells, paste0(base_dir, 'SA609_10x_samechip.csv'))
sce_combine <- load_data(meta_cells, base_dir, tag='SA609')
dim(sce_combine)
saveRDS(sce_combine, paste0(base_dir, 'SA609_samechip_rawdata.rds'))
cut_off_overall <- 0.025
if(cut_off_overall > 0){
  zero_cbind <- DelayedArray::rowMeans(assay(sce_combine, 'counts') == 0)
  sce_combine <- sce_combine[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
}
dim(sce_combine)
saveRDS(sce_combine, paste0(base_dir, 'SA609_samechip_filtered_genes.rds'))

datatag <- 'SA609_samechip'
rownames(sce_combine) <- rowData(sce_combine)$ID
colnames(sce_combine)[1]
rownames(sce_combine)[1]

res <- normalize_SCTransform(sce_combine, base_dir, datatag, return_data=T, output_fn=NULL)
rownames(res$sce)[1:10] 
rownames(res$sctransform)[1:5]
dim(rowData(res$sctransform))
colnames(rowData(res$sctransform))
sum(rownames(res$sce)==rownames(res$sctransform))
length(unique(rownames(res$sce)))
rownames(res$sce)[!rownames(res$sce) %in% rownames(res$sctransform)]
rownames(res$sctransform)[!rownames(res$sctransform) %in% rownames(res$sce)]

rownames(res$sctransform) <- rownames(res$sce) # need to confirm

norm_data <- res$sctransform
rowData(norm_data) <- dplyr::bind_cols(as.data.frame(rowData(res$sce)),
                                       as.data.frame(rowData(norm_data)))
?dplyr::bind_cols

dim(rowData(norm_data))
output_fn <- paste0(base_dir, datatag,'_sctransform_normalized.rds')
saveRDS(norm_data, output_fn)
dim(norm_data)



# SA604
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA604_rna/raw/'

fns <- list.files(base_dir)
fns <- fns[grepl('.rds',fns)]
fns <- gsub('.rds','',fns)
meta_cells <- data.frame(mouse_id=fns)
data.table::fwrite(meta_cells, paste0(base_dir, 'metasample_10x.csv'))
sce_combine <- load_data(meta_cells, base_dir, tag='SA604')
dim(sce_combine)

cut_off_overall <- 0.025
if(cut_off_overall > 0){
  zero_cbind <- DelayedArray::rowMeans(assay(sce_combine, 'counts') == 0)
  sce_combine <- sce_combine[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
}
dim(sce_combine)
datatag <- 'SA604'
res <- normalize_SCTransform(sce_combine, base_dir, datatag, return_data=T, output_fn=NULL)
rownames(res$sce)[1:10] 
rownames(res$sctransform)[1:5]
dim(rowData(res$sctransform))
colnames(rowData(res$sctransform))
sum(rownames(res$sce)==rownames(res$sctransform))
length(unique(rownames(res$sce)))
rownames(res$sce)[!rownames(res$sce) %in% rownames(res$sctransform)]
rownames(res$sctransform)[!rownames(res$sctransform) %in% rownames(res$sce)]

rownames(res$sctransform) <- rownames(res$sce) # need to confirm

norm_data <- res$sctransform
rowData(norm_data) <- dplyr::bind_cols(as.data.frame(rowData(res$sce)),
                                       as.data.frame(rowData(norm_data)))
?dplyr::bind_cols

dim(rowData(norm_data))
output_fn <- paste0(base_dir, datatag,'_sctransform_normalized.rds')
saveRDS(norm_data, output_fn)
dim(norm_data)
dim(colData(norm_data))
print(assayNames(norm_data))
output_fn <- paste0("/home/htran/storage/datasets/drug_resistance/rna_results/SA604_rna/normalized/", datatag,'_sctransform_normalized.rds')

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA604_rna/normalized/'
norm_sce <- readRDS(paste0(base_dir, 'SA604_sctransform_normalized.rds'))
dim(norm_sce)
t <- colData(norm_sce)
head(t)
library(SingleCellExperiment)
norm_sce2 <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(counts(norm_sce)), 
                                                                    logcounts=as.matrix(logcounts(norm_sce))),
                                                         colData = as.data.frame(colData(norm_sce)), 
                                                         rowData = as.data.frame(rowData(norm_sce)))

saveRDS(norm_sce2, paste0(base_dir, 'SA604_sctransform_normalized_v2.rds'))
data.table::fwrite(as.matrix(counts(norm_sce)),file = paste0(base_dir, 'SA604_sctransform_normalized_counts_mtx.csv.gz'),row.names = T, quote = F)
data.table::fwrite(as.matrix(logcounts(norm_sce)),file = paste0(base_dir, 'SA604_sctransform_normalized_logcounts_mtx.csv.gz'),row.names = T, quote = F)
data.table::fwrite(as.data.frame(colData(norm_sce)),file = paste0(base_dir, 'SA604_sctransform_colData.csv.gz'),row.names = T, quote = F)
data.table::fwrite(as.data.frame(rowData(norm_sce)),file = paste0(base_dir, 'SA604_sctransform_rowData.csv.gz'),row.names = T, quote = F)


?SingleCellExperiment::SingleCellExperiment

df <- data.frame(c1=c(1,3,4,7), c2=c(4,2,6,10))
p <- ggplot(df) +
  geom_point(aes(x = c1, y = c2),size = 7, alpha=0.6)
p <- p + annotate(geom="text", x=df$c1, y=df$c2, label="*",color="red", size=10)
p


# TO DO 
# Load data for 1 DE analysis SA609 
# Get incis, intrans genes 
# Get significant custom pathways 
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/cis_trans_scran/'

de_genes <- data.table::fread(paste0(base_dir,'signif_genes_FDR0.01.csv')) %>% as.data.frame()
dim(de_genes)

unique(de_genes$classified_gene)
de_genes$gene_symbol
deg_df <- de_genes
# deg_df <- deg_df %>%
#   dplyr::filter(classified_gene=='in_cis')
dim(deg_df)
datatag <- 'SA609'
# [1]           
# [2] 
# [3]         
# [4] 
# [5] 
gsea_out <- get_custom_pathway_results(deg_df,                      # named vector of statistical significance 
                                       desc, base_name=datatag,   
                                       pathway_name=c('custom_pathways'), #'cosmic' or 'cisplatin_resistance', or 'metastasis'
                                       groups_use=c("UTT-R","UUU-H"),    # vector of 2 elements, 2 group name used for DE analysis
                                       save_dir=base_dir,
                                       reference_genes_set="CISPLATIN_RESISTANCE")
gsea_out$pathway
rm(gsea_out)
gsea_out$reference_genes_set
gsea_out$datatag
class(gsea_out)




signf_ls <- tibble::tibble()
fgs <- c("COSMIC","CISPLATIN_RESISTANCE","CORE_FITNESS","BS_ESSENTIAL_CANCER") #,"METASTASIS"
for(reference_genes_set in fgs){
  for(fn in pair_groups$result_fn){
    # fn <- gsub('_logfc_results.csv','',fn)
    de_genes <- data.table::fread(paste0(base_dir, fn,'/','signif_genes_FDR0.01.csv')) %>% as.data.frame()
    print(dim(de_genes))
    gsea_out <- get_custom_pathway_results(de_genes,                      # named vector of statistical significance 
                                           fn, base_name=datatag,   
                                           pathway_name='custom_pathways', 
                                           groups_use=NULL,    # vector of 2 elements, 2 group name used for DE analysis
                                           save_dir=paste0(base_dir, fn,'/'),
                                           reference_genes_set)
    if(!is.null(gsea_out)){
      signf_ls <- dplyr::bind_rows(signf_ls, gsea_out)  
    }
    
    
  }
}

signf_ls <- as.data.frame(signf_ls)
signf_ls <- signf_ls %>%
  select(pathway, padj, desc, datatag, everything())
data.table::fwrite(signf_ls, paste0(base_dir,datatag,'_reference_sets.csv'))
colnames(signf_ls)



View(signf_ls)
"/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/cis_trans_scran/scrande_SA609_3_SA609_UTTT_R_UUUU_H_logfc_results.csv/signif_genes_FDR0.01.csv"

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA1035-v6/'
samples <- list.files(input_dir, pattern = '*.rdata')
samples <- gsub('.rdata','',samples)
cells_use <- colnames(sce)
cells_use[1]

load_data_testing <- function(samples, cells_use, input_dir, tag='CIS_Rx_RxH'){
  sce_list <- list()
  c <- 0
  for (f in unique(samples)){
    norm_fn <- paste0(input_dir,f,".rdata")
    if(file.exists(norm_fn)){
      sce_normalized <- readRDS(norm_fn)
      print(dim(sce_normalized))
      # sce_normalized <- sce_normalized[,sce_normalized$Barcode %in% cells_use]
      print(dim(sce_normalized))
      c <- c + 1
      if(c==1){
        cols_use <- colnames(colData(sce_normalized))
      }else{
        cols_use <- intersect(cols_use, colnames(colData(sce_normalized)))
      }
    } else{
      warning('Error loading files')
      print(f)
    }  
    print(paste0('DEBUG: ',f,'  ',dim(sce_normalized)[1],' ',dim(sce_normalized)[2]))
    sce_list[[f]] <- sce_normalized
    
  }  
  length(cols_use)
  
  sce_combine <- sce_cbind_func_v2(sce_list, cut_off_overall = 0, exprs = c("counts"), 
                                   colData_names = cols_use, save_raw=F, save_dir=input_dir, tag) 
  sce_combine <- sce_cbind_func_v2(sce_list, cut_off_overall = 0.1, exprs = c("counts"), 
                                   colData_names = cols_use, save_raw=F, save_dir=input_dir, tag='SA1035') 
  print(dim(sce_combine))
  # sce_combine <- readRDS(paste0(output_dir,'combined_total_genes_filtered.rds'))
  # Twice scater normalization
  # output are saved in logcounts exp values
  print(("Normalizing total data..."))
  
  ## cell name = library_id + barcode
  sce_combine$library_id <- gsub('.cache/','',sce_combine$Sample)
  sce_combine$library_id <- gsub('/filtered_feature_bc_matrix','',sce_combine$library_id)
  print(unique(sce_combine$library_id))
  # sum(unique(sce_combine$library_id) %in% meta_data$library_id)
  colnames(sce_combine) <- paste0(sce_combine$sample,'_',sce_combine$Barcode)
  # sce_combine <- sce_combine[,colnames(sce_combine) %in% meta_data$cell_id]
  # rownames(sce_combine)
  # print(rowData(sce_combine_raw)$Symbol[1:3])
  
  saveRDS(sce_combine, file=paste0(input_dir,"combined_",tag,".rds"))
  
  return(sce_combine)
}

