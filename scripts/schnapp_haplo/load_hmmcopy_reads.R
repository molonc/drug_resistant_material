
suppressPackageStartupMessages({
  require("optparse")
  require("data.table")
  require("feather")
})
option_list <- list(make_option(c("-l", "--library_ids"), type="character", 
                                default=NULL, help="library_ids", metavar="character"),
                    make_option(c("-s", "--samples"), type="character",
                                default=NULL, help="sample_id", metavar="character"),
                    make_option(c("-i", "--input_dir"), type="character", 
                                default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-o", "--output_file"), type="character", 
                                default=NULL, help="output_file", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
# print(opt$library_ids)
print(opt$input_dir)
print(opt$output_file)
print('Samples are: ')
# print(opt$samples)

library(parallel)
cores_use <- 8
nbcores <- detectCores()
if(cores_use > nbcores){
  cores_use <- nbcores
}

get_sample_id <- function(cell_ids, cores_use=8) {
  
  labels <- mclapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  }, mc.cores = cores_use)
  return(as.character(labels))
}

# initial.options <- commandArgs(trailingOnly = FALSE)
# file.arg.name <- "--file="
# script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
# script.basename <- dirname(script.name)
# print(script.basename)
# source(paste0(script.basename, "/utils.R"))

load_hmmcopy_reads <- function(samples, library_ids, input_dir, output_file){
  library_ids = strsplit(library_ids, ",")[[1]]
  samples = strsplit(samples, ",")[[1]]
  
  meta_df <- data.frame(library_id=library_ids, sample=samples,
                        row.names = library_ids,
                        stringsAsFactors=F, check.names=F)
  # rownames(meta_df) <- meta_df$library_id
  print(meta_df)
  save_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  # segment_file <- list.files(download_dir, pattern = "*_segments.csv", recursive = TRUE)
  # length(segment_file)  #should be 17 files
  # segment_files <- lapply(segment_file, function(f) {
  #   c(paste0(download_dir,f))
  # })
  # head(segment_file)
  
  for(f in library_ids){
    f1 = paste0(input_dir,'/',f,'/hmmcopy/',f,'_reads.csv')
    if(!file.exists(f1)){
      print(paste0('Double check library id: ',f1))
      stop("ERROR, check downloaded data in hmmcopy folder")
    }
  }
 
  
  
  raw_segments <- dplyr::bind_rows(lapply(library_ids, function(f) {
    f1 <- paste0(input_dir,'/',f,'/hmmcopy/',f,'_reads.csv')
    print(f1)
    seg_tmp <- fread(f1)
    print(paste0("Before filtering: ", dim(seg_tmp)[1], dim(seg_tmp)[2]))
    # print(colnames(seg_tmp))
    seg_tmp$sample <- get_sample_id(seg_tmp$cell_id, cores_use=6)
    print(seg_tmp$sample[1])
    print(seg_tmp$cell_id[1])
    seg_tmp <- seg_tmp[seg_tmp$sample==as.character(meta_df[f,'sample']),]
    seg_tmp <- seg_tmp[ ,c("chr", "start", "end", "cell_id", "state", "copy","reads")]
    print(paste0("After filtering: ",dim(seg_tmp)))
    print(paste0("f: ",f,' coln:  '))
    print(colnames(seg_tmp))
    saveRDS(seg_tmp,file = paste0(save_dir,f,'_cn.rds'))
    seg_tmp
  }))
  print(colnames(raw_segments))
  
  # input_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/SA535'
  # f <- 'A96201B'
  # save_dir <- '/home/htran/storage/datasets/metastasis_results/pseudobk_SA535/ascn_SA535/'
  # seg_tmp <- readRDS(paste0(save_dir,f,'_cn.rds'))
  # dim(seg_tmp)
  # colnames(raw_segments)[which(names(raw_segments) == "cell_id")] <- "cell_names"
  # colnames(raw_segments)[which(names(raw_segments) == "state")] <- "copy_number"
  # grouping_df <- read.csv(library_grouping, header=T, check.names=F, stringsAsFactors = FALSE)
  
  
  saveRDS(raw_segments,file = output_file)
  # split_raw_segments_by_sample(save_dir, raw_segments)
  # return(raw_segments)
}

load_hmmcopy_reads(opt$samples, opt$library_ids, opt$input_dir, opt$output_file)


# save_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_filtered/clonealign/'
