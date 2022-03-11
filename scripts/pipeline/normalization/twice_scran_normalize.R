suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  # require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  require("scran")
  require("Seurat")
  
})

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(paste0(script.basename, "/utils/normalize_utils.R"))

option_list <- list(make_option(c("-l", "--library_ids"), type="character", default=NULL, help="library_ids", metavar="character"),
                    make_option(c("-i", "--input_dir"), type="character", default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="output_file", metavar="character"),
                    make_option(c("-d", "--datatag"), type="character", default='SA', help="basename", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt$library_ids)
print(opt$input_dir)
print(opt$output_file)
print(opt$datatag)



twice_scran_normalize <- function(library_ids_ls, input_dir, output_file, datatag='SA'){
  library_ids = strsplit(library_ids_ls, ",")[[1]]
  print(library_ids)
  output_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(output_dir)){
    dir.create(output_dir)
  }
  # base_name <- basename(output_file)
  # base_name <- gsub('_normalized_output.rds','',base_name)
  # print(paste0("base_name is: ", base_name))
  
  mouse_id <- c()
  pdxid <- c()
  passage <- c()
  batch_info <- c()
  library_label <- c()
  treatmentSt <- c()
  cphases <- c()
  sce_list <- list()
  c <- 0
  
  # f <- 'TENX069'
  # download_10x_dir <- '/home/htran/storage/datasets/drug_resistance_RNAseq/SA609_human/'
  # sce_raw <- readRDS(paste0(download_10x_dir, '/', f,'/',f,'.rdata'))
  # colnames(sce_raw)[1:3]
  # rownames(sce_raw)[1:3]
  # dim(sce_raw)
  # View(head(rowData(sce_raw)))
  # length(rowData(sce_raw)$ID)
  # colnames(sce_raw) <- sce_raw$Barcode
  # rownames(sce_raw) <- rowData(sce_raw)$ID
  # length(unique(rowData(sce_raw)$ID))
  # rowData(sce_raw)$ID[1:5]
  # # sce_raw$Barcode[1:3]
  # # colnames(sce_raw) <- sce_raw$Barcode
  # # rownames(sce_raw)
  # # dim(sce_raw)
  # saveRDS(sce_raw, file = paste0(download_10x_dir, '/', f,'/',f,'.rdata'))
  # input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609_rna/filtered'
  # # output_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609_rna/normalized/'
  # sce_filtered <- readRDS(paste0(input_dir,'/',f,'_filtered.rds'))
  # dim(sce_filtered)
  # colnames(sce_filtered)[1:3]
  # rownames(sce_filtered)[1:3]
  # rowData(sce_filtered)$ID[1:5]
  # colnames(sce_filtered) <- sce_filtered$Barcode
  # rownames(sce_filtered) <- rowData(sce_filtered)$ID
  # saveRDS(sce_filtered, file = paste0(input_dir,'/',f,'_filtered.rds'))
  
  for (f in library_ids){
    # fn <- str_sub(f, 1, str_length(f)-19)  #-16
    print(paste0("Processing file:  ",f))
    norm_fn <- paste0(output_dir,f,"_normalized.rds")
    if(file.exists(norm_fn)){
      sce_normalized <- readRDS(norm_fn)
    } else{
      sce <- readRDS(paste0(input_dir, '/', f,'_filtered.rds'))
      print(dim(sce))
      if(dim(sce)[2]>0){
        # colnames(sce) <- paste0(f,"_",colnames(sce))
        colnames(sce) <- paste0(sce$mouse_id,colnames(sce))
        # print("Normalizing data")
        sce_normalized <- sce_normalize_size_factors(sce, min_size=80, rlog=FALSE, exprs="counts")
        print(dim(sce_normalized))
        saveRDS(sce_normalized, file = norm_fn)
      }
    }  
    if(dim(sce_normalized)[2]>0){  
      mouse_id <- c(mouse_id, sce_normalized$mouse_id)
      pdxid <- c(pdxid, sce_normalized$pdxid)
      passage <- c(passage, sce_normalized$passage)
      treatmentSt <- c(treatmentSt, sce_normalized$treatmentSt)
      print(paste0('DEBUG: ',f, ': cell cycle: ',length(sce_normalized$cell_cycle_phases)))
      cphases <- c(cphases, sce_normalized$cell_cycle_phases)
      library_label <- c(library_label, sce_normalized$library_label)
      batch_info <- c(batch_info, sce_normalized$batch_info)
      c <- c + 1
      sce_list[c] <- sce_normalized
    } else{
      print(paste0("Please double check library id: ",f))
    }
    print(paste0('DEBUG: ',f,'  ',dim(sce_normalized)[1],' ',dim(sce_normalized)[2]))
  }
  # print("mouse_id \n")
  # print(length(mouse_id))
  # print("pdxid \n")
  # print(length(pdxid))
  # print("passage \n")
  # print(length(passage))
  # print("treatmentSt \n")
  # print(length(treatmentSt))
  # print("library_label \n")
  # print(length(library_label))
  # print("batch_info \n")
  # print(length(batch_info))
  # print("cphases \n")
  # print(length(cphases))
  meta_data <- list(mouse_id=mouse_id,pdxid=pdxid,passage=passage,treatmentSt=treatmentSt,
                    library_label=library_label, batch_info=batch_info, cphases=cphases)
  print(length(sce_list))
  print("DEBUG")
  print(length(meta_data))
  
 
  meta_genes <- data.frame(gene_ens=rownames(rowData(sce_list[[1]])),gene_symb=rowData(sce_list[[1]])$Symbol, 
                           row.names = rownames(rowData(sce_list[[1]])), stringsAsFactors = F)
  write.csv(meta_genes, file = paste0(output_dir,'meta_genes.csv'), quote = F, row.names = F)
  print("DEBUG")
  print(dim(meta_genes))
  # Combine all data into 1 total data 
  # From Nick: gene should express in at least 10% of population
  print(paste0("Create normalized sce list with length: ",length(sce_list)))
  # sce_combine_10percent <- sce_cbind_func(sce_list, cut_off_overall = 0.1, exprs = c("counts", "normcounts"),
  #                                         colData_names = NULL, meta_data)
  
  # Raw data
  print("Combining...")
  sce_combine_raw <- sce_cbind_func(sce_list, cut_off_overall = 0, exprs = c("counts", "normcounts"),
  colData_names = NULL, meta_data)
  print(dim(sce_combine_raw))
  rowData(sce_combine_raw)$Symbol <- meta_genes[rownames(rowData(sce_combine_raw)),'gene_symb']
  # print(rowData(sce_combine_raw)$Symbol[1:3])
  saveRDS(sce_combine_raw, file = paste0(output_dir, datatag, "_raw.rds"))
  
  sce_combine <- sce_cbind_func(sce_list, cut_off_overall = 0.1, exprs = c("counts", "normcounts"),
                                colData_names = NULL, meta_data)
  print(dim(sce_combine))
  print(("Save combined data"))
  print(assayNames(sce_combine))
  rowData(sce_combine)$Symbol <- meta_genes[rownames(rowData(sce_combine)),'gene_symb']
  saveRDS(sce_combine, file = paste0(output_dir, datatag, "_normcounts_normalized.rds"))
  
  
  # Twice scater normalization
  # output are saved in logcounts exp values
  print(("Normalizing total data..."))
  sce_normalized_total <- sce_normalize_size_factors(sce_combine, min_size=300, rlog=TRUE, exprs="normcounts")
  saveRDS(sce_normalized_total, file = output_file)
  if(file.exists(output_file)){
    return(TRUE)
  } else{
    return(FALSE)
  }
  # return(sce_normalized_total)
  
}

# library_ids = c("SCRNA10X_SA_CHIP0077_001","SCRNA10X_SA_CHIP0077_002")
# input_dir = '/home/htran/storage/datasets/metastasis_results/rnaseq/filtered'
# output_file = "/home/htran/storage/datasets/metastasis_results/rnaseq/filtered/normalized_output.rds"



twice_scran_normalize(opt$library_ids, opt$input_dir, opt$output_file, opt$datatag)





# testing sth
# libs <- c("SCRNA10X_SA_CHIP0071_000","TENX071")
# download_dir <- "/home/htran/storage/datasets/drug_resistance_RNAseq/human"
# sce1 <- readRDS(paste0(download_dir, '/', libs[1],'/',libs[1],'.rdata'))
# sce2 <- readRDS(paste0(download_dir, '/', libs[2],'/',libs[2],'.rdata'))
# inter_cells <- intersect(colnames(sce1), colnames(sce2))
# length(inter_cells)
# dim(sce1)
# dim(sce2)