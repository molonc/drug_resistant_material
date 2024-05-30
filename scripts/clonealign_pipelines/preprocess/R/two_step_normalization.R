library(SingleCellExperiment)
library(scater)
library(scran)

library(argparse)
library(yaml)

library(stringr)
library(tidyverse)


# this script from Nick for twice scran normalization of all samples together

parser <- ArgumentParser(description = "Two step scran normalization")
parser$add_argument('--input_csv', metavar='FILE', type='character',
                    help="A csv file with all the inputs")
parser$add_argument('--outdir', type = 'character', metavar = 'FILE',
                    help="Output path for normalized SCEs.")
parser$add_argument('--meta_file', type = 'character', metavar = 'FILE',
                    help="Meta file with series etc information")
args <- parser$parse_args()

meta <- read_yaml(args$meta_file)

paths <- read.csv(args$input_csv, sep=",", header=TRUE)
rdata <- as.character(paths$path)
samples <- as.character(paths$sample_id)
outdir <- args$outdir
print(paste0("Creating dir ", outdir))
dir.create(outdir, showWarnings = TRUE)

sces <- vector("list", length(rdata)-1)
saved_rowdata <- vector("list", length(rdata)-1)
i <- 1
for (i in 1:length(rdata)) {
  print(rdata[i])
  sample_level <- readRDS(rdata[i])
  sample_level <- normalize(sample_level)
  sample_level$sample <- samples[i]
  ## readjusting the series
  sample_level$series <- meta[[samples[i]]]$series
  sces[[i]] = sample_level
  print(sces[[i]])
}

new_sce <- sces[[1]]
saved_rowdata[[1]] <- rowData(new_sce)
rowData(new_sce) <- NULL
colData(new_sce)$tenx <- NULL
for (j in 2:length(rdata)) {
  print(j)
  saved_rowdata[[j]] <- rowData(sces[[j]])
  rowData(sces[[j]]) <- NULL
  # MA: 24 Nov 2020: if there were more than 1 tenx libraries, then there will be one extra colData called tenx
  #   Removing it here, so I can cbind
  colData(sces[[j]])$tenx <- NULL
  #colData(sces[[j]])$treat <- NULL
  colData(sces[[j]])$label <- NULL
  new_sce <- cbind(new_sce, sces[[j]])
}
#saveRDS(new_sce,file=paste0(outdir,"/bound.rdata"))

qclust <- quickCluster(new_sce, min.size = 100)
print("clustered")
new_sce <- computeSumFactors(new_sce, clusters = qclust)
new_sce$size_factor <- sizeFactors(new_sce)
new_sce <- normalize(new_sce)
#saveRDS(new_sce, file=paste0(outdir,"/renormed.rdata"))

j <- 1
for (sample in unique(new_sce$sample)) {
  print(sample)
  sample_level <- new_sce[,new_sce$sample==sample]
  rowData(sample_level) <- saved_rowdata[[j]]
  print(sample_level)
  saveRDS(sample_level, file=paste0(outdir,"/",sample,".rdata"))
  j <- j+1
}
