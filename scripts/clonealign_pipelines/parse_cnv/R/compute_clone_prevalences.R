#' Formats segment file

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)
library(feather)

parser <- ArgumentParser(description = "Select genes")
parser$add_argument('--segment_files', metavar='FILE', default=NULL, type='character', nargs ='+',
                    help="Segment files")
parser$add_argument('--feather_file', metavar='FILE', default=NULL, type='character',
                    help="Feather file")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Clone prevalences")
args <- parser$parse_args()

# MA: 20 June 2020: taking the info from the feather file

if(!is.null(args$segment_files)) {
  raw_segments <- dplyr::bind_rows(lapply(args$segment_files, function(f) {
    fread(f)
  }))
  
  segments <- raw_segments %>%
    dplyr::mutate(chr=paste0("chr", chr))
}

if (!is.null(args$feather_file)) {
  segments <- read_feather(args$feather_file)
}
  
prevalences <- segments %>%
  dplyr::select(cell_names, time, cluster) %>%
  unique %>%
  dplyr::group_by(cluster, time) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq=n/sum(n))

write.table(prevalences, file = args$outfname, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

cat("Completed.\n")

