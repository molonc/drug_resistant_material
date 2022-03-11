suppressPackageStartupMessages({
  require("scater")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  require("scran")
  
})

script_dir <- '/home/htran/Projects/farhia_project/rnaseq/pipeline/utils/'
source(paste0(script_dir, "normalize_utils.R"))

# rlog: log=TRUE, FALSE in logNormCounts function, Logical scalar indicating 
# whether normalized values should be log2-transformed.if rlog=FALSE: save data to normcounts, rlog=TRUE, save to logcounts
sce_normalize_size_factors <- function(sce, min_size=100, rlog=FALSE, exprs="counts", name){
  # library(scater)
  print("Quick clustering")
  if(min_size < dim(sce)[2]){
    qclust <- quickCluster(sce, min.size = min_size, assay.type=exprs)
    print("Compute sum factors")
    sce <- computeSumFactors(sce, clusters = qclust, assay.type=exprs)
    sce$size_factor <- sizeFactors(sce)
    print("Normalize data")
    # scater version 1.14.6
    # count --> normcounts
    # normcounts --> logcounts 
    
    # String containing an assay name for storing the output normalized values. 
    # Defaults to "logcounts" when log=TRUE and "normcounts" otherwise.
    sce_normalized <- logNormCounts(sce, log=rlog, exprs_values=exprs, name=name, size_factors=sce$size_factor)
    
  } else{
    print("Smaller nb cells than min size threshold in this library")
    sce_normalized <- logNormCounts(sce, log=rlog, exprs_values=exprs, name=name, size_factors=NULL)
  }
  
  # scater version < 1.14.6
  # sce_normalized <- normalize(sce, return_log=rlog, exprs_values=exprs)
  print(assayNames(sce_normalized))
  return(sce_normalized)
}

# sce_cbind_func <- function(sce_list, exprs = c("counts", "logcounts","normcounts"), 
#                               colData_names = NULL) {
#   # Combining expression
#   assay_list <- list()
#   for (i in seq_len(length(exprs))) {
#     assay_list[[i]] <- do.call(cbind, lapply(sce_list, 
#                                              function(y) assay(y, exprs[i])))
#   }
#   names(assay_list) <- exprs
#   
#   # Combining metadata
#   colData_list <- do.call(DelayedArray::rbind, 
#                           lapply(sce_list, function(y) colData(y)[, colData_names, drop = FALSE]))
#   sce_combine <- SingleCellExperiment::SingleCellExperiment(assay = assay_list, 
#                                                             colData = colData_list)
#   
#   print(paste0("Dim sce combine: ",dim(sce_combine)[1],' ',dim(sce_combine)[2]))
#   # print(colnames(colData(sce_combine)))
#   print(paste0("sce combine assay name: ",assayNames(sce_combine)))
#   return(sce_combine)
# }


twice_scran_normalize <- function(library_ids, input_dir, output_dir, datatag='SA'){
  #library_ids = strsplit(library_ids_ls, ",")[[1]]
  print(library_ids)
  #output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  sce_list <- list()
  saved_rowdata <- list()
  for (f in library_ids){
    print(paste0("Processing file:  ",f))
    
    sce <- readRDS(paste0(input_dir, '/', f,'.rds'))
    saved_rowdata <- c(saved_rowdata, rowData(sce))
    print(dim(sce))
    if(dim(sce)[2]>0){
      # First time normalize data using scran, input is counts exp, normalized output are kept in normcounts exp values 
      sce$sample <- f   # adding the name so I can save the individual files later
      # rlog=FALSE means the normalized counts are NOT log transformed
      sce_normalized <- sce_normalize_size_factors(sce, min_size=100, rlog=FALSE, exprs="counts", name="normcounts")  ##  rlog=FALSE: count --> normcounts
      if(!is.null(sce_normalized)){  
        sce_list <- c(sce_list, sce_normalized)
      } else{
        print(paste0("Please double check library id: ",f))
      }
    }
  }
  # Combining all sces to 1 total sce
  # Take out the tenx column
  sce_normalized$tenx <- NULL
  sce_combine <- sce_cbind_func(sce_list, exprs = c("counts", "normcounts"),
                                colData_names = colnames(colData(sce_normalized)))
  
  # Twice scater normalization
  # normalized output are kept in logcounts exp values
  print(("Normalizing total data..."))
  sce_normalized_total <- sce_normalize_size_factors(sce_combine, min_size=300, rlog=TRUE, exprs="normcounts", name="logcounts") ##  rlog=TRUE: normcounts --> logcounts
  saveRDS(sce_normalized_total, file = paste0(output_dir,"/allsce_normalized.rds"))
  # Saving individual files
  j <- 1
  for (sample in library_ids) {
    print(sample)
    sample_level <- sce_normalized_total[,sce_normalized_total$sample==sample]
    rowData(sample_level) <- saved_rowdata[[j]]
    print(sample_level)
    saveRDS(sample_level, file=paste0(output_dir,"/",sample,".rdata"))
    j <- j+1
  }
}  

twice_scran_normalize_v2 <- function(sce, input_dir, output_dir, datatag, return_data=F){
  #library_ids = strsplit(library_ids_ls, ",")[[1]]
  # print(library_ids)
  #output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  # sce <- sce_combine
  print(dim(sce))
  # colnames(colData(sce))
  # sce$sample[1:3]
  sce_list <- list()
  for (s in unique(sce$sample)){
    print(paste0("Processing sample:  ",s))
    
    # sce <- readRDS(paste0(input_dir, '/', f,'.rds'))
    sce_tmp <- sce[,sce$sample==s] 
    # saved_rowdata <- c(saved_rowdata, rowData(sce))
    print(dim(sce_tmp))
    if(dim(sce_tmp)[2]>0){
      # First time normalize data using scran, input is counts exp, normalized output are kept in normcounts exp values 
      # rlog=FALSE means the normalized counts are NOT log transformed
      sce_normalized <- sce_normalize_size_factors(sce_tmp, min_size=100, rlog=FALSE, exprs="counts", name="normcounts")  ##  rlog=FALSE: count --> normcounts
      if(!is.null(sce_normalized)){  
        sce_list <- c(sce_list, sce_normalized)
      } else{
        print(paste0("Please double check sample id: ",s))
      }
    }
  }
  # Combining all sces to 1 total sce
  # Take out the tenx column
  sce_normalized$tenx <- NULL
  sce_combine <- sce_cbind_func(sce_list, rowData(sce), cut_off_overall = 0, exprs = c("counts", "normcounts"),
                                colData_names = colnames(colData(sce_normalized)), meta_data=NULL)
  # colData(sce_combine)[1:3,1:3]
  # Twice scater normalization
  # normalized output are kept in logcounts exp values
  print(("Normalizing total data..."))
  sce_normalized_total <- sce_normalize_size_factors(sce_combine, min_size=300, rlog=TRUE, exprs="normcounts", name="logcounts") ##  rlog=TRUE: normcounts --> logcounts
  saveRDS(sce_normalized_total, file = paste0(output_dir,datatag,"_twice_scran_normalized_v3.rds"))
  if(return_data){
    return(sce_normalized_total)
  }  
}  




# Mirela
# twice_scran_normalize(gsub(".rds","", list.files("SA609_v6/sce_annotated/")), input_dir="SA609_v6/sce_annotated/", output_dir="SA609_v6/sce_twice_scran/")
# twice_scran_normalize(gsub(".rds","", list.files("SA1035_v6/sce_annotated/")), input_dir="SA1035_v6/sce_annotated/", output_dir="SA1035_v6/sce_twice_scran/")
# twice_scran_normalize(gsub(".rds","", list.files("SA535_v6/sce_annotated/")), input_dir="SA535_v6/sce_annotated/", output_dir="SA535_v6/sce_twice_scran/")


