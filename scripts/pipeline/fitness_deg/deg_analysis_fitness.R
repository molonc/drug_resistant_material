
# BiocManager::install("enrichplot")
suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  # require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
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


# Load sce data 
# Get sample id 
# DE analysis
results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_Tyler_v3/'
datatag <- 'SA535'
library_grouping <- read.csv(paste0(results_dir,'library_groupings.csv'), check.names=F, stringsAsFactors=F)
# View(head(library_grouping))
dim(library_grouping)
length(unique(library_grouping$sample))
sample_ids <- unique(library_grouping$sample)
length(sample_ids)

results_10x_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/'
output_file <- paste0(results_10x_dir,'deg_analysis/de_output.png')
input_dir <- paste0(results_10x_dir,'clonealign')
input_file <- paste0(input_dir,'/normalized_sce.rds')
sample_id <- ''
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results'

norm_data <- logcounts(sce)
class(norm_data)
norm_data <- as.data.frame(as.matrix(norm_data))
dim(norm_data)
z_score <- scale(t(norm_data))
dim(z_score)

summary(as.double(z_score[1,]))
View(z_score[1:10,1:10])
z_score <- t(z_score)
min(z_score)
z_score_backup <- z_score
rownames(z_score)[1:3]



input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/'

input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/deg_analysis/'
tag <- 'SA535_A_B'
convert_de_csv(tag, save_dir)

tag <- 'SA535_A_G'
convert_de_csv(tag, save_dir)

tag <- 'SA535_DE_B'
convert_de_csv(tag, save_dir)

tag <- 'SA535_DE_G'
convert_de_csv(tag, save_dir)

tag_ls <- c('SA535_A_UTTT_UTTTT','SA535_A_UT_UTTT',
            'SA535_B_UUU_UUUU','SA535_B_UU_UUU','SA535_B_UU_UUUU',
            'SA535_DE_UTTT_UTTTT','SA535_DE_UT_UTTT',
            'SA535_G_UUU_UUUU','SA535_G_UU_UUU','SA535_G_UU_UUUU')
for(tag in tag_ls){
  convert_de_csv(tag, save_dir)
}


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
  sce_backup <- sce
  print(summary(as.factor(sce$clone_id)))
  print(summary(as.factor(sce$treatmentSt)))
  genes_map_symb <- read.csv(paste0(input_dir, "biodatabase/meta_genes.csv"), header=T, stringsAsFactors=F)
  rownames(genes_map_symb) <- genes_map_symb$gene_ens
  dim(genes_map_symb)
  
  sce <- sce[,!sce$clone_id %in% c('unassigned')]
  # sce <- sce[,sce$clone_id %in% c('A','B')]
  
  sce$clone_id <- ifelse(sce$clone_id=="D-E","DE",sce$clone_id)
  srt <- as.Seurat(sce, counts = "counts", data = "logcounts") 
  dim(srt)
  meta_data <- srt@meta.data
  
  # Resistant clone: A, D-E
  # Sensitive clone: G, B
  groups_use_ls <- list(T1=c('A','G'), T2=c('A','B'),
                        T3=c('DE','G'), T4=c('DE','B'))
  
  treatment_st <- list(G=c("U","UU",'UUU','UUUU','UUUUU'),
                       B=c("U","UU",'UUU','UUUU','UUUUU'),
                       A=c("UT","UTT",'UTTT','UTTTT','UTU','UTTU','UTTTU'), 
                       DE=c("UT","UTT",'UTTT','UTTTT','UTU','UTTU','UTTTU'))
  
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
  # Resistant clone: A, D-E
  # Sensitive clone: G, B
  # treatment_st <- list(T1=c("UT","UTT"),
  #                      T2=c("UTT","UTTT"),
  #                      T3=c("UTTT","UTTTT"),
  #                      T4=c("UT","UTTT"),
  #                      T5=c("UTT","UTTTT"))
  
  treatment_st <- list(T1=c("UTT","UT"),
                       T2=c("UTTT","UTT"),
                       T3=c("UTTTT","UTTT"),
                       T4=c("UTTT","UT"),
                       T5=c("UTTTT","UTT"),
                       T6=c("UTTTT","UT"))
  
  groups_use_ls <- c('A','DE')
  
  # treatment_st <- list(T1=c("UU","UUU"),
  #                      T2=c("UUU","UUUU"),
  #                      T3=c("UUUU","UUUUU"),
  #                      T4=c("UU","UUUU"),
  #                      T5=c("UUU","UUUUU"))
  # 
  # groups_use_ls <- c('G','B')
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
  
}
  