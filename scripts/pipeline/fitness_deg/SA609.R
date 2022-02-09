#SA609


results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609_rna/'
sample_ids <- c('SA609X10XB02454','SA609X6XB03404','SA609X7XB03505',
                'SA609X4XB003083','SA609X5XB03230','',
                '','')
datatag <- 'SA609'

input_sce <- paste0(results_10x_dir,'clonealign/observed_sce.rds')
data_dir <- paste0(results_10x_dir,'clonealign/SA609_clonealign_X10/tenxqc/')
f <- 'SA609X10XB02454'
sce_normalized$cell_id[1:3]
rowData(sce_normalized)$ID[1:3]
rowData(sce_normalized)$Symbol[1:3]
sce_normalized$timepoint[1:3]
sce_normalized$treatment[1:3]
sce_normalized$lib



colnames(colData(sce_normalized))[which(colnames(colData(sce_normalized)) == "clone")] <- "clone_id"
sce_normalized$treatmentSt[1:3]


# Get clone info
clone_df_1 <- read.table(file =paste0(results_10x_dir,'cnv/fitness_cell_assignment_feb07_2020.tsv'), sep = ',', header = TRUE)
View(head(clone_df))
rownames(clone_df) <- clone_df$single_cell_id
clone_df <- clone_df[grep('SA609', clone_df$single_cell_id, value = T),]
dim(clone_df)
colnames(clone_df)[which(colnames(clone_df) == "single_cell_id")] <- "cell_id"
colnames(clone_df)[which(colnames(clone_df) == "letters")] <- "clone_id"
clone_df <- clone_df[,!(colnames(clone_df) %in% c('datatag','V1'))]
write.csv(clone_df, file=paste0(results_10x_dir,'clonealign/Mirela_output/total_clones.csv'), quote=F, row.names=F)
sample_ids <- unique(clone_df$id)

clone_df <- read.csv(paste0(results_10x_dir,'clonealign/Mirela_output/total_clones.csv'), check.names=F, stringsAsFactors=F)
dim(clone_df)
summary(as.factor(clone_df$clone_id))
colnames(clone_df)[which(colnames(clone_df) == "clone")] <- "clone_id"
ext <- clone_df$cell_id[!clone_df$cell_id %in% cells_use]
length(ext)
clone_df1 <- clone_df[ext,]
clone_df1 <- clone_df1[clone_df1$clone_id=='R',]
clone_df1$cell_id

dim(clone_df)
colnames(sce_normalized)[1:3]
clone_df$cell_id[1:3]
cells_use <- intersect(clone_df$cell_id, colnames(sce_normalized))  #16122 in 17473: output of clonealign, TO DO: change threshold of QC
length(cells_use)
rownames(clone_df) <- clone_df$cell_id
summary(as.factor(sce_normalized$clone))
sce_normalized$clone_id <- 'unassigned'
colData(sce_normalized)[cells_use,'clone_id'] <- clone_df[cells_use,'clone_id']
sce_normalized <- sce_normalized[,sce_normalized$clone_id != 'unassigned']
dim(sce_normalized)
saveRDS(sce_normalized, file=paste0(results_10x_dir,'normalized/sce_normalized_DE.rds'))

saveRDS(sce_normalized, file=paste0(results_10x_dir,'normalized/sce_normalized_total.rds'))
length(sce_normalized$cell_id)

sce_normalized <- readRDS(paste0(results_10x_dir,'clonealign/observed_sce.rds'))

mouse_id <- c()
passage <- c()
treatmentSt <- c()
sce_list <- list()
c <- 0
sample_ids <- sample_ids[sample_ids !="SA609X4XB003083"]
for(s in sample_ids){
  norm_fn <- paste0(data_dir,s,".rds")
  if(file.exists(norm_fn)){
    print(s)
    sce_normalized <- readRDS(norm_fn)
    print(dim(sce_normalized))
    sce_normalized <- sce_normalized[,sce_normalized$cell_id %in% clone_df$cell_id]
    print(dim(sce_normalized))
    colnames(sce_normalized) <- sce_normalized$cell_id
    
    rownames(sce_normalized) <- rowData(sce_normalized)$ID
    print(rowData(sce_normalized)$Symbol[1:3])
    c <- c + 1
    sce_list[c] <- sce_normalized
    print(assayNames(sce_normalized))
    print(dim(sce_normalized))
    mouse_id <- c(mouse_id, rep(f,dim(sce_normalized)[2]))
    passage <- c(passage, sce_normalized$timepoint)
    treatmentSt <- c(treatmentSt, sce_normalized$treatment)
  }
}
 

output_file <- paste0(results_10x_dir,'deg_analysis/de_output.png')
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results'
input_file <- paste0(results_10x_dir,'clonealign/observed_sce.rds')

execute_seurat_DEG_between_clones <- function(sample_id, clonealign_fn, input_dir, input_file, output_file, cluster_rm=NULL,
                                              test_use = 'wilcox',
                                              pAdjust = 0.05,   #pAdjust = 0.0125, minlfc = 0.25
                                              minlfc = 0.25, datatag='SA'){
  input_dir=paste0(input_dir,'/')
  save_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  # save_dir_deg <- paste0(save_dir,'deg/')
  # if (!file.exists(save_dir_deg)){
  #   dir.create(save_dir_deg)
  # }
  # base_name <- basename(output_file)
  # base_name <- gsub('_pathway.rds','',base_name)
  # print(paste0("base_name is: ", base_name))
  sce <- readRDS(input_file)
  print("Initialized sce")
  print(dim(sce))
  # sce$clone <- 'unassigned'
  # sce <- sce[,sce$cluster_label!=cluster_rm & sce$Grouping=="Primary"]
  
  # rowData(sce)$Symbol <- genes_map_symb[rownames(sce),'gene_symb']
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  sum(mito_genes==TRUE)
  
  ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
  sum(ribo_genes==TRUE)
  sce <- sce[(!mito_genes) & (!ribo_genes), ]
  print(dim(sce))
  
  # sce_backup <- sce
  # print(summary(as.factor(clone_df$clone_id)))
  colnames(colData(sce))[which(colnames(colData(sce)) == "clone")] <- "clone_id"
  print(summary(as.factor(sce$clone_id)))
  sce1 <- sce[,sce$treatmentSt=='UT']
  print(summary(as.factor(sce1$clone_id)))
  genes_map_symb <- read.csv(paste0(input_dir, "biodatabase/meta_genes.csv"), header=T, stringsAsFactors=F)
  rownames(genes_map_symb) <- genes_map_symb$gene_ens
  dim(genes_map_symb)
  
  # sce <- sce[,!sce$clone_id %in% c('unassigned')]
  srt <- as.Seurat(sce, counts = "counts", data = "logcounts") 
  meta_data <- srt@meta.data
  # sce1 <- sce[,sce$clone_id == 'H']
  # summary(sce1)
  # sce$clone_id <- ifelse(sce$clone_id=="D-E","DE",sce$clone_id)
  dim(sce)
  
  # Resistant clone: A, D-E
  # Sensitive clone: G, B
  #SA609 
  
  groups_use_ls <- list(T1=c('R','H'))
  groups_use_ls <- list(T1=c('R','E'))
  treatment_st <- list(R=c("UTTTT"),E=c("UUUUUUUU"))
  # treatment_st <- list(R=c("UTTT"),H=c("UUUUUUUU"))
  # treatment_st <- list(R=c("UTT"),H=c("UUUUUUUU"))
  treatment_st <- list(R=c("UTTTT"),H=c("UUUUUUUU"))
  # treatment_st <- list(R=c("UTTTT"),E=c("UUUUUUUU"))
  # treatment_st <- list(R=c("UT"),E=c("UUUUUUUU"))  srt <- as.Seurat(sce, counts = "counts", data = "logcounts") 
  dim(srt)
  
  colnames(meta_data)
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
      
      pathway_ls <- calculate_DE_analysis_v2(base_name=paste0(datatag,'_',treatment_st[[groups_use[1]]],'_',treatment_st[[groups_use[2]]]), 
                                             meta_data=meta_data, 
                                             cells_use_g1, cells_use_g2,
                                             feature_use="clone_id",
                                             srt, genes_map_symb, groups_use,
                                             test_use, save_dir, input_dir,
                                             pAdjustThrs=pAdjust, minLogFC=minlfc,
                                             nbtopup = 30, nbtopdown = 30, 
                                             save_data=T, viz=F)
      
      
    }  
  }
  
  
  # Resistant clone: A, D-E
  # Sensitive clone: G, B
  treatment_st <- list(T1=c("UTT","UT"),
                       T2=c("UTTT","UTT"),
                       T3=c("UTTTT","UTTT"),
                       T4=c("UTTT","UT"),
                       T5=c("UTTTT","UTT"),
                       T6=c("UTTTT","UT"))
  
  groups_use_ls <- c('R')
  
  
  
  for(groups_use in groups_use_ls){
    for(ts in treatment_st){
      print(paste0('groups_use: ',groups_use,' ts1: ',ts[1], ' ts2: ', ts[2]))
      cells_use_g1 <- rownames(meta_data)[meta_data$clone_id==groups_use & meta_data$treatmentSt %in% ts[1]]
      cells_use_g2 <- rownames(meta_data)[meta_data$clone_id==groups_use & meta_data$treatmentSt %in% ts[2]]
      if(length(cells_use_g1) < 100 || length(cells_use_g2) < 100){
        print(groups_use)
        print("There are no cells or small nb cells
            only which satisfy the input condition ")
        print(paste0("\n Observed clones: ",groups_use,
                     "  nb cells in ts: ",ts[1], ":",length(cells_use_g1),
                     "  nb cells in ts: ",ts[2], ":",length(cells_use_g2)))
        
      } else{
        
        pathway_ls <- calculate_DE_analysis_v2(base_name=paste0(datatag,'_',groups_use), meta_data=meta_data,
                                               cells_use_g1, cells_use_g2,
                                               feature_use="treatmentSt",
                                               srt, genes_map_symb, ts,
                                               test_use, save_dir, input_dir,
                                               pAdjustThrs=pAdjust, minLogFC=minlfc,
                                               nbtopup = 30, nbtopdown = 30,
                                               save_data=T, viz=F)
        
        
      }
    }
    
  }
  
  
  
  treatment_st <- list(T1=c("UTTTT","UT"))
  for(ts in treatment_st){
    print(paste0('groups_use:  ts1: ',ts[1], ' ts2: ', ts[2]))
    cells_use_g1 <- rownames(meta_data)[meta_data$treatmentSt == ts[1]]
    cells_use_g2 <- rownames(meta_data)[meta_data$treatmentSt == ts[2]]
    if(length(cells_use_g1) < 100 || length(cells_use_g2) < 100){
      print("There are no cells or small nb cells
            only which satisfy the input condition ")
      print(paste0("\n Observed clones: ",
                   "  nb cells in ts: ",ts[1], ":",length(cells_use_g1),
                   "  nb cells in ts: ",ts[2], ":",length(cells_use_g2)))
      
    } else{
      
      pathway_ls <- calculate_DE_analysis_v2(base_name=datatag, meta_data=meta_data,
                                             cells_use_g1, cells_use_g2,
                                             feature_use="treatmentSt",
                                             srt, genes_map_symb, ts,
                                             test_use, save_dir, input_dir,
                                             pAdjustThrs=pAdjust, minLogFC=minlfc,
                                             nbtopup = 30, nbtopdown = 30,
                                             save_data=T, viz=F)
      
      
    }
  }
  
  
  
  
}
