suppressPackageStartupMessages({
  require("optparse")
  # require("scater")
  # require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  # require("scran")
  require("Seurat")
  
})

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
print(script.basename)
source(paste0(dirname(script.basename), "/utils/normalize_utils.R"))

option_list <- list(make_option(c("-l", "--library_ids"), type="character", default=NULL, help="library_ids", metavar="character"),
                    make_option(c("-i", "--input_dir"), type="character", default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-s", "--suffix"), type="character", default='_f_dt_m.rds', help="suffix of filtered data", metavar="character"),
                    make_option(c("-d", "--datatag"), type="character", default='SA', help="basename", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt$library_ids)
print(opt$input_dir)
print(opt$output_file)
print(opt$suffix)
print(opt$datatag)


## suffix: f: filtered - quality control, dt: doublet removed, m: mutation detection and removed.
SCTransform_normalize <- function(library_ids_ls, input_dir, output_file, 
                                  suffix='_f_dt_m.rds', datatag='SA'){
  library_ids = as.character(unlist(strsplit(library_ids_ls, ",")))
  print(library_ids)
  output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  # base_name <- basename(output_file)
  # base_name <- gsub('_normalized_output.rds','',base_name)
  # print(paste0("base_name is: ", base_name))
  
  sce_list <- list()
  for (f in library_ids){
    # fn <- str_sub(f, 1, str_length(f)-19)  #-16
    
    filtered_fn <- paste0(input_dir, '/', f, suffix)
    print(paste0("Processing file:  ",filtered_fn))
    if(file.exists(filtered_fn)){
      fsce <- readRDS(filtered_fn)
      print(dim(fsce))
      if(dim(fsce)[2]>0){  
        obs_features <- c('seqnames',"ranges","strand","start","end","width","element")
        obs_features <- intersect(obs_features, colnames(rowData(fsce)))
        rowData(fsce)[,obs_features] <- NULL
        
        sce_list[[f]] <- fsce
      } else{
        stop(paste0("Have no filtered data, please double check library id: ",f))
      }  
    } 
    
  }
  # 'strand' %in% colnames(rowData(fsce))
  print(length(sce_list))
  # Combine all data into 1 total data 
  # From Nick: gene should express in at least 10% of population
  print(paste0("Create filtered sce list with length: ",length(sce_list)))
  # sce_combine_10percent <- sce_cbind_func(sce_list, cut_off_overall = 0.1, exprs = c("counts", "normcounts"),
  #                                         colData_names = NULL, meta_data)
  
  # Raw data
  print("Combining...")
  
  ## Genes should expressed in at least 1% of total cells
  sce_combine <- sce_cbind_func_v2(sce_list, cut_off_overall = 0.025, exprs = c("counts"), 
                   colData_names = NULL, save_raw=T, save_dir=input_dir, tag=datatag) 
    
  print(dim(sce_combine))
  # sce_combine <- readRDS(paste0(output_dir,'combined_total_genes_filtered.rds'))
  # Twice scater normalization
  # output are saved in logcounts exp values
  print(("Normalizing total data..."))
  
  ## cell name = library_id + barcode
  sce_combine$library_id <- gsub('.cache/','',sce_combine$Sample)
  sce_combine$library_id <- gsub('/filtered_feature_bc_matrix','',sce_combine$library_id)
  print(unique(sce_combine$library_id))
  colnames(sce_combine) <- paste0(sce_combine$library_id,'_',sce_combine$Barcode)
  if('ID' %in% colnames(rowData(sce_combine))){
    saveRDS(sce_combine, paste0(input_dir,'/sce_combine.rds'))
    rownames(sce_combine) <- rowData(sce_combine)$ID # using ens gene id as index  
  }else{
    saveRDS(sce_combine, paste0(input_dir,'/sce_combine.rds'))  
    stop('Please double check ens gene id in rowData!!!')
  }
  
  normalize_SCTransform(sce_combine, output_dir, datatag, return_data=F, output_file)
  sce <- readRDS(output_file)
  dim(sce)
  dim(counts(sce))
  dim(logcounts(sce))
  
  # meta_cells <- data.table::fread('/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/snakemake_10x_v2/SA535-10x_updated_jan_2022.csv') %>% as.data.frame()
  # dim(meta_cells)
  # colnames(meta_cells)
  # meta_cells$Noted <- NULL
  # meta_cells$lid <- gsub('SCRNA10X_SA_CHIP0','L',meta_cells$library_id)
  # viz_umap(norm_sce, output_dir, datatag, return_seurat_obj=F)
  # srt@meta.data <- srt@meta.data %>% left_join(meta_cells, by=c('library_id'))
  # srt@meta.data$lid <- paste0('C',srt@meta.data$lid)
  # srt@meta.data$lid <- as.factor(srt@meta.data$lid)
  # srt@meta.data$Site_origin <- as.factor(srt@meta.data$Site_origin)
  # sce$cell_id <- colnames(sce)
  # cells_features <- as.data.frame(colData(sce)) %>%
  #   left_join(meta_cells, by=c('library_id'))%>%
  #   DataFrame(check.names = FALSE)
  # rownames(cells_features) <- cells_features$cell_id
  # sce@colData <- cells_features
  # sce$pdxid[1]
}



SCTransform_normalize(opt$library_ids, opt$input_dir, opt$output_file, 
                      opt$suffix, opt$datatag)




# library_ids = c("SCRNA10X_SA_CHIP0077_004","SCRNA10X_SA_CHIP0176_001","SCRNA10X_SA_CHIP0213_001"
#                 ,"SCRNA10X_SA_CHIP0213_002","SCRNA10X_SA_CHIP0220_001")
# input_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered'
# output_dir <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/normalized/'
# datatag <- 'SA535'
# output_file <- paste0(output_dir,'SA535_sctransform_normalized.rds')
# suffix <- '_f_dt_m.rds'

# norm1 <- readRDS(paste0(output_dir,'SA535_sctransform_normalized_v1.rds'))
# dim(norm1)
# library_ids = c("SCRNA10X_SA_CHIP0077_001","SCRNA10X_SA_CHIP0077_002")
# input_dir = '/home/htran/storage/datasets/metastasis_results/rnaseq/filtered'
# output_file = "/home/htran/storage/datasets/metastasis_results/rnaseq/filtered/normalized_output.rds"

## Process meta data
# meta_fn <- '/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/snakemake_10x_v2/SA535-10x_updated_march.csv'
# meta_data <- read.csv(meta_fn, stringsAsFactors=F, check.names=F)
# dim(meta_data)
# View(meta_data)
# meta_data <- meta_data %>%
#   dplyr::rename(sample_id=mouse_id)
# 
# meta_data <- meta_data %>%
#   dplyr::select(-Noted)
# dim(norm_sce)
# 
# norm_sce <- readRDS(output_file)
# saveRDS(norm_sce, output_file)
# norm_sce$cell_id <- colnames(norm_sce)
# norm_sce$cell_id[1]
# dim(meta_cells)
# class(colData(norm_sce))
# class(meta_cells)
# dim(colData(norm_sce))
# norm_sce@colData <- as.data.frame(colData(norm_sce)) %>%
#   inner_join(meta_data, by=c('library_id'))%>%
#   DataFrame(check.names = FALSE)
# 
# colnames(norm_sce) <- norm_sce$cell_id
# rownames(norm_sce)[1]
# dim(rowData(fsce))
# s <- colnames(rowData(fsce))
# s1 <- colnames(rowData(norm_sce))
# s[!s %in% s1]
# dim(colData(norm_sce))
# colnames(rowData(norm_sce))
# rownames(norm_sce)[1]
