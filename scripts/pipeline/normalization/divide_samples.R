suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  require("argparse")
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  
})

option_list <- list(make_option(c("-s", "--mouse_id"), type="character", default=NULL, help="library_ids", metavar="character"),
                    make_option(c("-i", "--input_file"), type="character", default=NULL, help="input_dir", metavar="character"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="output_file", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
# print(opt$mouse_id)
print(opt$input_file)
print(opt$output_file)

divide_normalized_data_by_samples <- function(input_file, output_file){
  save_dir <- paste0(dirname(output_file),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  normalized_sce <- readRDS(input_file)
  samples <- unique(normalized_sce$mouse_id)
  for(s in samples){
    normalized_sce_tmp <- normalized_sce[,normalized_sce$mouse_id==s]
    print(dim(normalized_sce_tmp))
    saveRDS(normalized_sce_tmp,file = paste0(save_dir,s,'.rds'))
  }
  saveRDS(normalized_sce_tmp,file = output_file)
}

divide_normalized_data_by_samples(opt$input_file, opt$output_file)



