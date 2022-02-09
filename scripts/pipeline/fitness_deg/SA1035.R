
# Resistant clone: G, H
# Sensitive clone: E
save_dir <- "/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/deg_analysis/"
results_10x_dir <- "/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/"
datatag <- 'SA1035'
input_file <- paste0(results_10x_dir,'normalized/','SA1035_normalized_output.rds')
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results'
input_dir=paste0(input_dir,'/')
sce <- readRDS(input_file)
print("Initialized sce")
print(dim(sce))

mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
sum(mito_genes==TRUE)

ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
sum(ribo_genes==TRUE)
sce <- sce[(!mito_genes) & (!ribo_genes), ]
print(dim(sce))
sce_backup <- sce


sce$clone_id <- 'unassigned'
clone_df <- read.csv(paste0(results_10x_dir,'clonealign/','total_clones.csv'), check.names=F, stringsAsFactors=F)
View(head(clone_df))
summary(as.factor(clone_df$clone_id))
rownames(clone_df) <- clone_df$cell_id
colnames(colData(sce))


cells_use <- intersect(clone_df$cell_id, colnames(sce))  #16122 in 17473: output of clonealign, TO DO: change threshold of QC
length(cells_use)
colData(sce)[cells_use,'clone_id'] <- clone_df[cells_use,'clone_id']
summary(as.factor(sce$clone_id))
dim(sce)
saveRDS(sce, file=paste0(results_10x_dir,'normalized/','SA1035_normalized_clones.rds'))

genes_map_symb <- read.csv(paste0(input_dir, "biodatabase/meta_genes.csv"), header=T, stringsAsFactors=F)
rownames(genes_map_symb) <- genes_map_symb$gene_ens
dim(genes_map_symb)

srt <- as.Seurat(sce, counts = "counts", data = "logcounts") 
dim(srt)
meta_data <- srt@meta.data
colnames(meta_data)

groups_use_ls <- list(T1=c('G','E'))  #T1=c('G','E'), T2=c('H','E')

treatment_st <- list(E=c('UUUU','UUUUU'),  #"U","UU",'UUU',
                     G=c('UTTT','UTTTT'), #"UT","UTT",,'UTU',,'UTTU','UTTTU'
                     H=c('UTTT','UTTTT'))  #"UT","UTT",'UTU','UTTU','UTTTU'



for(groups_use in groups_use_ls){
  cells_use_g1 <- rownames(meta_data)[meta_data$clone_id==groups_use[1] & meta_data$treatmentSt %in% treatment_st[[groups_use[1]]]]
  cells_use_g2 <- rownames(meta_data)[meta_data$clone_id==groups_use[2] & meta_data$treatmentSt %in% treatment_st[[groups_use[2]]]]
  if(length(cells_use_g1) < 100 || length(cells_use_g2) < 100){
    print(groups_use)
    print("There are no cells or small nb cells 
            only which satisfy the input condition ")
    print(paste0("\n Observed clones: ",groups_use[1],' vs ',groups_use[2], 
                 "  nb cells in ",groups_use[1], ":",length(cells_use_g1),
                 "  nb cells in ",groups_use[2], ":",length(cells_use_g2)))
    
  } else{
    
    pathway_ls <- calculate_DE_analysis_v2(base_name=datatag, meta_data=meta_data, 
                                           cells_use_g1, cells_use_g2,
                                           feature_use="clone_id",
                                           srt, genes_map_symb, groups_use,
                                           test_use, save_dir, input_dir,
                                           pAdjustThrs=pAdjust, minLogFC=minlfc,
                                           nbtopup = 30, nbtopdown = 30, 
                                           save_data=T, viz=F)
    
    
  }  
}