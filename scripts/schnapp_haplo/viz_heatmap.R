suppressPackageStartupMessages({
  require("optparse")
  require("data.table")
  require("feather")
  require("dplyr")
})
# option_list <- list(make_option(c("-l", "--library_ids"), type="character", 
#                                 default=NULL, help="library_ids", metavar="character"),
#                     make_option(c("-i", "--input_dir"), type="character", 
#                                 default=NULL, help="input_dir", metavar="character"),
#                     make_option(c("-o", "--output_file"), type="character", 
#                                 default=NULL, help="output_file", metavar="character"))

option_list <- list(make_option(c("-a", "--allele_counts_fn"), type="character", 
                                default=NULL, help="allele_counts_fn", metavar="character"),
                    make_option(c("-i", "--total_bincnv_fn"), type="character", 
                                default=NULL, help="total_bincnv_fn", metavar="character"),
                    make_option(c("-f", "--filtered_cells_fn"), type="character", 
                                default=NULL, help="filtered_cells_fn", metavar="character"),
                    make_option(c("-o", "--output_fn"), type="character", 
                                default=NULL, help="output_file", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
print(script.basename)
source(paste0(script.basename, "/schnapp_utils.R"))
