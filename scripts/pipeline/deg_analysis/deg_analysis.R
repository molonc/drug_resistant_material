
# BiocManager::install("enrichplot")
suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  # require("argparse")
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
  # require("grid")
  
})

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(paste0(script.basename, "/utils/deg_utils.R"))

option_list <- list(make_option(c("-d", "--input_dir"), type="character", default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-i", "--input_file"), type="character", default=NULL, help="input_file", metavar="character"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-p", "--p_adjust"), type="double", default=0.0125, help="p_adjust"),
                    make_option(c("-m", "--min_logfc"), type="double", default=0.25, help="min_logfc"),
                    make_option(c("-c", "--cluster_rm"), type="character", default=NULL, help="exception_cluster"),
                    make_option(c("-b", "--datatag"), type="character", default='SA', help="basename", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



# input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/'
# input_file <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/clonealign/SA535X9XB03616/SA535X9XB03616_sce_filtered.rds'
# output_file <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/clonealign/SA535X9XB03616/SA535X9XB03616_pathway.rds'
# clonealign_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/clonealign/SA535X9XB03616/SA535X9XB03616_clones_045.csv'
# execute_seurat_DEG_between_clones(clonealign_fn, input_dir, input_file, output_file, cluster_rm,
#                                               test_use,pAdjust, minlfc)
# obs_clones <- c('E','H')
# input_deg_fn
# sample_id <- 'SA535X9XB03616'
# plot_DE_genes_chr(input_deg_fn=input_deg_fn, clones=obs_clones, sample_name=sample_id,
#                               additional_genes = NULL, n_genes_to_annotate = 25)
  

# output_file <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/pathway.rds'
execute_seurat_DEG_between_clones <- function(input_dir, input_file, output_file, cluster_rm=NULL,
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
  print("Initialized sce...")
  sce <- readRDS(input_file)
  print(dim(sce))
  print("Removing replicates...")
  replicates <- c('SA609X4XB03084','SA609X5XB03235','SA609X7XB03573')
  sce <- sce[,!sce$sample %in% replicates]
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
  
  print(summary(as.factor(sce$clone)))
  print(summary(as.factor(sce$treatmentSt)))
  # sce$clone <- ifelse(sce$clone=="A_C",'A',
  #              ifelse(sce$clone=="B_F",'B',sce$clone))
  
  # For SA1035
  # sce$clone <- ifelse(sce$clone=="C_D",'D',sce$clone)
  
  # For SA535
  # sce$clone <- ifelse(sce$clone=="C_D",'D',sce$clone)
  # sce$clone <- ifelse(sce$clone=="G",'GF',
  #                     ifelse(sce$clone=="F",'GF',sce$clone))
  
  genes_map_symb <- read.csv(paste0(input_dir, "biodatabase/meta_genes.csv"), header=T, stringsAsFactors=F)
  rownames(genes_map_symb) <- genes_map_symb$gene_ens
  
  
  sce <- sce[,!sce$clone %in% c('unassigned','Un','un')]
  # sce <- sce[,sce$clone %in% c('B','E')]
  
  # groups_use_ls <- list(X34=c("X3","X4"),X37=c("X3","X7"),X47=c("X4","X7"))
  # groups_use_ls <- list()
  # 
  # clones <- unique(meta_data$clone)
  # clones <- clones[clones != 'unassigned']
  # for(c1 in rep(1:(length(clones)-1), 1)){
  #   for(c2 in rep(2:(length(clones)), 1)){
  #     if(clones[c1]!=clones[c2]){
  #       lb <- paste0(clones[c1],clones[c2])
  #       groups_use_ls[[lb]] <- c(clones[c1],clones[c2])
  #     }
  #   }  
  # }
  # SA1035
  # groups_use_ls <- list(T1=c('B','D'),T2=c('B','A'),T3=c('E','D'),T4=c('E','A'))
  # groups_use_ls <- list(T1=c('B','D'),T2=c('E','D'))
  # treatment_st <- list(B=c("UT","UTT",'UTTT','UTTTT'),  #,'UTU','UTTU','UTTTU'
  #                      E=c("UT","UTT",'UTTT','UTTTT'),  #,'UTU','UTTU','UTTTU'
  #                      D=c("U","UU",'UUU','UUUU','UUUUU'),
  #                      A=c("U","UU",'UUU','UUUU','UUUUU'),
  #                      C=c("UTU",'UTTU','UTTTU'))
  # treatment_st <- list(B=c('UTU','UTTU','UTTTU'),  #,
  #                      E=c('UTU','UTTU','UTTTU'),  #,'UTU','UTTU','UTTTU'
  #                      D=c("U","UU",'UUU','UUUU','UUUUU'),
  #                      A=c("U","UU",'UUU','UUUU','UUUUU'),
  #                      C=c("UTU",'UTTU','UTTTU'))
  # groups_use_ls <- list(T1=c('B','D'),T3=c('E','D'))
  
  
  # SA535
  # groups_use_ls <- list(T1=c('GF','D'), T2=c('E','D'),T3=c('GF','C'),T4=c('E','C'))
  # treatment_st <- list(D=c("U","UU",'UUU','UUUU','UUUUU'),C=c("U","UU",'UUU','UUUU','UUUUU'),  #,'UTU','UTTU','UTTTU'
  #                      GF=c("UT","UTT",'UTTT','UTTTT'), #'UTU','UTTU','UTTTU'
  #                      E=c("UT","UTT",'UTTT','UTTTT'))  #,'UTU','UTTU','UTTTU'
  # treatment_st <- list(D=c("U","UU",'UUU','UUUU','UUUUU'),C=c("U","UU",'UUU','UUUU','UUUUU'),  #,'UTU','UTTU','UTTTU'
  #                      GF=c('UTU','UTTU','UTTTU'), #
  #                      E=c('UTU','UTTU','UTTTU'))  #,'UTU','UTTU','UTTTU'
  
  
  #SA609 
  groups_use_ls <- list(T1=c('R','C'))
  treatment_st <- list(R=c("UTTT"),C=c("UUUU"))
  # groups_use_ls <- list(T1=c('R','C'))
  # treatment_st <- list(R=c("UTTT"),C=c("U","UUU"))
  # 
  # treatment_st <- list(R=c("UTTT"),C=c("UUUU"))
  # treatment_st <- list(R=c("UTTT"),C=c("UTT"))
  # treatment_st <- list(R=c("UTTT"),C=c("UT"))
  # treatment_st <- list(R=c("UT"),E=c("UUUUUUUU"))
  
  # groups_use <- groups_use_ls[[2]]
  
  # treatment_ls <- list(T1=c('UTTT','UT'), T2=c('UTTTT','UT'))
  
  srt <- as.Seurat(sce, counts = "counts", data = "logcounts") 
  dim(srt)
  meta_data <- srt@meta.data
  for(groups_use in groups_use_ls){
    cells_use_g1 <- rownames(meta_data)[meta_data$clone==groups_use[1] & meta_data$treatmentSt %in% treatment_st[[groups_use[1]]]]
    cells_use_g2 <- rownames(meta_data)[meta_data$clone==groups_use[2] & meta_data$treatmentSt %in% treatment_st[[groups_use[2]]]]
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
                                             feature_use="clone",
                                             srt, genes_map_symb, groups_use,
                                             test_use, save_dir, input_dir,
                                             pAdjustThrs=pAdjust, minLogFC=minlfc,
                                             nbtopup = 30, nbtopdown = 30, 
                                             save_data=T, viz=F)
      
      
    }  
  }
  
  
  
  # for(groups_use in treatment_ls){
  #   cells_use_g1 <- rownames(meta_data)[meta_data$treatmentSt==groups_use[1]]
  #   cells_use_g2 <- rownames(meta_data)[meta_data$treatmentSt==groups_use[2]]
  #   if(length(cells_use_g1) < 100 || length(cells_use_g2) < 100){
  #     print(groups_use)
  #     print("There are no cells or small nb cells 
  #           only which satisfy the input condition ")
  #     print(paste0("\n Observed clones: ",groups_use[1],' vs ',groups_use[2], 
  #                  "  nb cells in ",groups_use[1], ":",length(cells_use_g1),
  #                  "  nb cells in ",groups_use[2], ":",length(cells_use_g2)))
  #     
  #   } else{
  #     
  #     pathway_ls <- calculate_DE_analysis_v2(base_name=datatag, meta_data=meta_data, 
  #                                            cells_use_g1, cells_use_g2,
  #                                            feature_use="treatmentSt",
  #                                            srt, genes_map_symb, groups_use,
  #                                            test_use, save_dir, input_dir,
  #                                            pAdjustThrs=pAdjust, minLogFC=minlfc,
  #                                            nbtopup = 30, nbtopdown = 30, 
  #                                            save_data=T, viz=F)
  #     
  #     
  #   }  
  # }
    

  # saveRDS(pathway_ls, file = output_file)
  
}


execute_seurat_DEG_between_treatment_status <- function(clonealign_fn, input_dir, 
                                                        input_file, output_file, 
                                                        cluster_rm=NULL, test_use = 'wilcox',
                                                        pAdjust = 0.0125, minlfc = 0.5, 
                                                        datatag='SA'){
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
  print("Raw sce: ")
  print(dim(sce))
  print(summary(as.factor(sce$clone)))
  
  genes_map <- read.csv(paste0(input_dir, "biodatabase/meta_genes.csv"), header=T, stringsAsFactors=F)
  rownames(genes_map) <- genes_map$gene_ens
  dim(genes_map)
  # groups_use_ls <- list(X34=c("X3","X4"),X37=c("X3","X7"),X47=c("X4","X7"))
  
  obs_clones <- unique(sce$clone)
  obs_clones <- obs_clones[obs_clones != 'unassigned']
  cluster_rm <- c('C','D') 
  obs_clones <- obs_clones[!obs_clones %in% cluster_rm]
  print(obs_clones)
  sce <- sce[,sce$clone_id %in% obs_clones]
  print("Observed sce: ")
  print(dim(sce))
  mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")
  sum(mito_genes==TRUE)
  
  ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")  # or ^RP[L|S]?
  sum(ribo_genes==TRUE)
  sce <- sce[(!mito_genes) & (!ribo_genes), ]
  print("Observed sce: ")
  print(paste0('Removing mito and ribo genes: ',dim(sce1)[1],'_',dim(sce1)[2]))
  
  srt <- as.Seurat(sce, counts = "counts", data = "logcounts") 
  meta_data <- srt@meta.data
  # ts_ls <- unique(sce$treatmentSt)
  group_use_dh <- list(T23=c("UT","UTU"),T33=c("UTT","UTU"),T44=c("UTTT","UTTU"), T55=c("UTTTT","UTTTU"),
                       T53=c("UTTTT","UTU"),T43=c("UTTT","UTU"))
  group_use_ts <- list(T23=c("UT","UTT"),T34=c("UTT","UTTT"),T24=c("UT","UTTT"), 
                             T35=c("UTT","UTTTT"),T25=c("UT","UTTTT"),T45=c("UTTT","UTTTT"))
  
  
  group_use_dh <- list(T53=c("UTTTT","UTU"),T43=c("UTTT","UTU"))
  # for(c1 in rep(1:(length(clones)-1), 1)){
  #   for(c2 in rep(2:(length(clones)), 1)){
  #     if(clones[c1]!=clones[c2]){
  #       lb <- paste0(clones[c1],clones[c2])
  #       groups_use_ls[[lb]] <- c(clones[c1],clones[c2])
  #     }
  #   }  
  # }
  # groups_use_ls <- do.call(c, list(group_use_ts, group_use_dh))
  
  for(groups_use in group_use_ts){
    meta_data_tmp <- meta_data
    print(dim(meta_data_tmp))
    cells_use_g1 <- rownames(meta_data_tmp)[meta_data_tmp$treatmentSt==groups_use[1]]
    cells_use_g2 <- rownames(meta_data_tmp)[meta_data_tmp$treatmentSt==groups_use[2]]
    if(length(cells_use_g1) < 100 & length(cells_use_g2) < 100){
      print("There are no cells which satisfy the input condition or small nb cells only")
    } else{
      pathway_ls <- calculate_DE_analysis_v2(datatag, meta_data_tmp, 
                                             cells_use_g1, cells_use_g2,
                                             feature_use="treatmentSt",
                                             srt, genes_map, groups_use,
                                             test_use, save_dir, input_dir,
                                             pAdjust, minlfc)
    }
    
  }
  
  for(groups_use in group_use_dh){
    meta_data_tmp <- meta_data
    print(dim(meta_data_tmp))
    cells_use_g1 <- rownames(meta_data_tmp)[meta_data_tmp$treatmentSt==groups_use[1]]
    cells_use_g2 <- rownames(meta_data_tmp)[meta_data_tmp$treatmentSt==groups_use[2]]
    if(length(cells_use_g1) < 100 & length(cells_use_g2) < 100){
      print("There are no cells which satisfy the input condition or small nb cells only")
    } else{
      pathway_ls <- calculate_DE_analysis_v2(datatag, meta_data_tmp, cells_use_g1, cells_use_g2,
                                             feature_use="treatmentSt",
                                             srt, genes_map, groups_use,
                                             test_use, save_dir, input_dir,
                                             pAdjust, minlfc)
    }
    
  }
  
  
  saveRDS(pathway_ls, file = output_file)
  
}


print(opt$input_file)
print(opt$output_file)
print(opt$input_dir)
print(opt$p_adjust)
print(opt$min_logfc)
# print(opt$cluster_rm)
print(opt$datatag)
# cluster_rm <- c('C','D')
# execute_seurat_DEG_between_treatment_status(opt$input_dir, opt$input_file, opt$output_file, cluster_rm,
#                                   'wilcox', opt$p_adjust,
#                                   opt$min_logfc, opt$datatag)

# execute_seurat_DEG_between_primary(opt$input_dir, opt$input_file, opt$output_file, opt$cluster_rm,
#                                    test_use = 'wilcox', pAdjust = opt$p_adjust,
#                                    minlfc = opt$min_logfc)
# results_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA919'
# input_file <- paste0(results_10x_dir,'/clonealign/normalized_clones.rds')



input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6'
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna'
output_file <- paste0(save_dir,'/deg_analysis/SA609_deg_pathway.rds')

input_file <- paste0(results_10x_dir,'/total_sce_clones.rds')
datatag <- 'SA609'
execute_seurat_DEG_between_clones(input_dir, input_file, output_file, 
                                  cluster_rm=NULL,
                                  'wilcox', 0.05, 0.25, datatag)  #pAdjust = 0.0125, minlfc = 0.25
                                   

# edgeR
pair_groups_fn <- paste0(input_dir,'/SA609_rna/deg_analysis/pair_groups.csv')

View(pair_groups)

pair_groups$desc <- paste0(pair_groups$datatag,'_',
                           pair_groups$sample1,'-',pair_groups$clone1,'_',
                           pair_groups$sample2,'-',pair_groups$clone2)
pair_groups$nbcells_g1 <- NA
pair_groups$nbcells_g2 <- NA

# results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna'
# input_file <- paste0(results_10x_dir,'/clonealign/Mirela_output_v1/normalized_clones.rds')
# output_file <- paste0(results_10x_dir,'/deg_analysis/SA1035_deg_pathway.rds')
# input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results'
# datatag <- 'SA1035'

# results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total'
# input_file <- paste0(results_10x_dir,'/clonealign/normalized_clones_SA535_cispl.rds')
# output_file <- paste0(results_10x_dir,'/deg_analysis/SA1035_deg_pathway.rds')
# input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results'
# datatag <- 'SA535'


