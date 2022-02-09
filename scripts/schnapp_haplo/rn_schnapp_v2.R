suppressPackageStartupMessages({
  require("optparse")
  require("data.table")
  # require("feather")
  require("dplyr")
  require("schnapps")
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
print('Hello World')
filter_data_v2(opt$allele_counts_fn, opt$total_bincnv_fn, opt$filtered_cells_fn, opt$output_fn,
                           return_data=F, tag='filtered')

# run_alleleSC_v2(opt$allele_counts_fn, opt$total_bincnv_fn, opt$output_fn, prefix='ascn')
# allele_counts_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_pseudobulk/ascn_SA1035/raw_allele_data.rds'
# total_bincnv_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_pseudobulk/ascn_SA1035/raw_cnbins.rds'
# filtered_cells_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/corrupt_grow/filtered_cells.txt'
# output_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_pseudobulk/ascn_SA1035/filtered_allele_data.rds'
# filter_data_v2(allele_counts_fn, total_bincnv_fn, 
#                filtered_cells_fn, output_fn,
#                return_data=F, tag='filtered')
