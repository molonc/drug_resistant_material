
suppressPackageStartupMessages({
  require("optparse")
  require("data.table")
  require("feather")
  require("dplyr")
})
option_list <- list(make_option(c("-l", "--library_ids"), type="character", 
                                default=NULL, help="library_ids", metavar="character"),
                    make_option(c("-s", "--samples"), type="character",
                                default=NULL, help="sample_id", metavar="character"),
                    make_option(c("-i", "--input_dir"), type="character", 
                                default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-t", "--lib_ticket_fn"), type="character", 
                                default=NULL, help="lib_ticket_filename", metavar="character"),
                    make_option(c("-o", "--output_file"), type="character", 
                                default=NULL, help="output_file", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
# print(opt$library_ids) 
print(opt$lib_ticket_fn)
print(opt$input_dir)
print(opt$output_file)
print('Samples are: ')
# print(opt$samples)

# library(parallel)
# cores_use <- 8
# nbcores <- detectCores()
# if(cores_use > nbcores){
#   cores_use <- nbcores
# }


# initial.options <- commandArgs(trailingOnly = FALSE)
# file.arg.name <- "--file="
# script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
# script.basename <- dirname(script.name)
# print(script.basename)
# source(paste0(script.basename, "/utils.R"))
get_allele_counts_v2 <- function(samples, library_ids, input_dir, output_file,
                                 lib_ticket_fn, prefix='total'){
  save_dir <- paste0(dirname(output_file),'/')
  input_dir <- paste0(input_dir,'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  
  library_ids = strsplit(library_ids, ",")[[1]]
  samples = strsplit(samples, ",")[[1]]
  print('DEBUG_allele_load data \n')
  print(library_ids)
  print(samples)
  
  meta_df <- data.frame(library_id=library_ids, sample=samples,
                        stringsAsFactors=F, check.names=F)
  
  lib_ticket_df <- read.csv(lib_ticket_fn, stringsAsFactors=F, check.names=F)
  # View(lib_ticket_df)
  meta_df <- meta_df %>% inner_join(lib_ticket_df, by = 'library_id')
  rownames(meta_df) <- meta_df$library_id
  print(colnames(meta_df))
  print(dim(meta_df))
  
  
  allele_counts_ls <- list()
  # keep_features <- c("allele_id", "chromosome", "end","hap_label","readcount","start","cell_id")
  for(library_id in meta_df$library_id){
    search_dir = paste0(input_dir, meta_df[library_id,'jira_ticket'],'/results/count_haps/sample_',
                        meta_df[library_id,'sample'],'/')
    allele_fn_tsv <- paste0(search_dir, 'allele_counts.tsv') # old standard
    allele_fn_csv <- paste0(search_dir, 'allele_counts.csv') # new standard?
    
    # In case of gz file
    # allele_fn_gz <- paste0(search_dir, 'allele_counts.csv.gz')
    # zz=gzfile(allele_fn_gz,'rt')  
    # dat=read.csv(zz,header=F) 
    if(dir.exists(search_dir)){
      if(file.exists(allele_fn_tsv)){
        data <- as.data.frame(data.table::fread(allele_fn_tsv))
      } else if(file.exists(allele_fn_csv)){
        data <- as.data.frame(data.table::fread(allele_fn_csv, sep=','))
      } else{
        data <- NULL
        print(paste0("**** ERROR, check allele output for lib id: ",library_id, ' sample id: ',meta_df[library_id,'sample']))
      }
      
      # sample_id <- gsub('sample_','',dirname(f))
      # data$library_id <- rep(library_id, dim(data)[1])
      # data$sample_id <- rep(sample_id, dim(data)[1])
      
      if(!is.null(data)){
        print(dim(data))
        # data <- data[,colnames(data) %in% keep_features]
        allele_counts_ls[[library_id]] <- data
      }
      
      
    }
    else{
      print(paste0("**** ERROR, check allele output for lib id: ",library_id, ' sample id: ',meta_df[library_id,'sample']))
    }
   
   
  }
  
  total_allele_counts <- dplyr::bind_rows(allele_counts_ls)
  print(dim(total_allele_counts))
  colnames(total_allele_counts)[which(colnames(total_allele_counts) == "chromosome")] <- "chr"
  # total_allele_counts <- total_allele_counts[total_allele_counts$cell_id %in% cells_use,]
  # print(dim(total_allele_counts))
  saveRDS(total_allele_counts, file = output_file)  #paste0(save_dir,prefix,'_allele_data.rds')
  print("Completed!")
  
}


get_allele_counts_v2(opt$samples, opt$library_ids, opt$input_dir, opt$output_file,
                     opt$lib_ticket_fn, 'total')



# filenames <- list.files(searchpath, pattern=fn, recursive = TRUE)
# if(!is.null(filenames)){
#   for(f in filenames){
#     
#     if(grep(meta_df[library_id,'sample'], f, value = F)==1){
#       print(f)
#       print(meta_df[library_id,'sample'])
#       print(paste0(searchpath, f))
#      
#     }
#    
#   }

# input_dir <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535/pseudobulk/SC-3712/results/count_haps/sample_SA535X8XB03431/'
# df <- data.table::fread(paste0(input_dir,'allele_counts.csv'), sep=',')
# dim(df)
# View(head(df))
# 
# input_dir <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535/pseudobulk/SC-3799/results/count_haps/sample_SA535X7XB03305/'
# df1 <- data.table::fread(paste0(input_dir,'allele_counts.tsv'))
# dim(df1)
# View(head(df1))