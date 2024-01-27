suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(parallel)
  library(dplyr)
  library(ggplot2)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
})
get_library_id <- function(cell_ids, cores_use=2) {
  labels <- mclapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  }, mc.cores = cores_use)
  return(as.character(labels))
}

get_reads_per_clone <- function(input_dir, cellclones_fn, output_file, 
                                obs_clones=NULL,tag=NULL, cores_use=NULL){
  if(is.null(cores_use)){
    cores_use <- 5
  }
  nbcores <- detectCores()
  if(cores_use > nbcores){
    cores_use <- nbcores
  }
  output_dir <- paste0(dirname(output_file),'/')
  if(!file.exists(output_dir)){
    dir.create(output_dir)
  }
  print("Read clone labels file")
  # input_dir <- paste0(input_dir,'/')
  print(input_dir)
  # input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA535/'
  cell_clones <- read.delim(cellclones_fn, stringsAsFactors = FALSE, sep = ",")
  print("DEBUG 1")
  print(dim(cell_clones))
  print("Remove unassigned cells")
  cell_clones <- cell_clones[!cell_clones$clone_id %in% c('Un','unassigned','Unassigned','None','un'),]
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  obs_libs <- get_library_id(cell_clones$cell_id)
  obs_libs <- unique(obs_libs)
  print(obs_libs)
  summary_wgd <- list()
  for(l in obs_libs){
    # lib_fn1 <- paste0(input_dir,l,'/',l,'/annotation/',l,'_metrics.csv.gz') # csv.gz or csv
    # lib_fn2 <- paste0(input_dir,l,'/',l,'/annotation/',l,'_metrics.csv') # csv.gz or csv
    # lib_fn3 <- paste0(input_dir,l,'/',l,'/annotation/','metrics.csv.gz') # csv.gz or csv
    # 
    lib_fn1 <- paste0(input_dir,l,'/annotation/',l,'_metrics.csv') # csv.gz or csv
    lib_fn2 <- paste0(input_dir,l,'/annotation/',l,'_metrics.csv.gz') # csv.gz or csv
    lib_fn3 <- paste0(input_dir,l,'/annotation/','metrics.csv.gz') # csv.gz or csv
    lib_fn <- ifelse(file.exists(lib_fn1),lib_fn1,
                     ifelse(file.exists(lib_fn2),lib_fn2,
                      ifelse(file.exists(lib_fn3),lib_fn3,'')))
    print(lib_fn)
    if(lib_fn!=''){
      tmp <- data.table::fread(lib_fn)
      tmp <- tmp %>%
        dplyr::filter(quality>=0.75)
      pct_wgd <- round(100 * sum(tmp$state_mode==4)/dim(tmp)[1],2)
      print(paste0(l, ' - ',pct_wgd, '%'))
      summary_wgd[[l]] <- pct_wgd
    }
    else{
      print('!!!!!!!!!Error, double check library: ')
      print(l)
    }
  } 
  
  stat_df <- tibble::tibble(library_id=names(summary_wgd), percentage_wgd=as.numeric(unlist(summary_wgd)))
  print(dim(stat_df))
  stat_df$tag <- tag
  print(round(mean(stat_df$percentage_wgd),2))
  print(round(sd(stat_df$percentage_wgd),2))
  data.table::fwrite(stat_df, paste0(output_dir,tag,'_percentage_WGD.csv.gz'))
  # stat_df <- data.table::fread(paste0(output_dir,tag,'_percentage_WGD.csv.gz'))
  # Number of whole genome duplication cells per library in each series: 
  ## avg=2.18%, sd=2.02% for SA609, avg 2.22%, sd=1.18% for SA535,
  ## avg=4.23%, sd=1.89% for SA1035
  
}  


cellclones_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cell_clones/SA609_cell_clones_metadata.csv'
# cellclones_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cell_clones/SA609_cell_clones.csv'
input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/'
output_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones_v3/'
output_file <- paste0(output_dir,'out.csv')
tag <- 'SA609'
cores_use=NULL
stat1 <- data.table::fread(paste0(output_dir,tag,'_percentage_WGD.csv.gz'))

cellclones_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cell_clones/SA535_cell_clones.csv'
input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA535_segment/'
output_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA535_segment/reads_clones_v2/'
output_file <- paste0(output_dir,'out.csv')
tag='SA535'
cores_use=NULL
stat2 <- data.table::fread(paste0(output_dir,tag,'_percentage_WGD.csv.gz'))


cellclones_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cell_clones/SA1035_cell_clones.csv'
input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA1035/'
output_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA1035/reads_clones_v2/'
output_file <- paste0(output_dir,'out.csv')
tag='SA1035'
cores_use=NULL
stat3 <- data.table::fread(paste0(output_dir,tag,'_percentage_WGD.csv.gz'))

dim(stat1)
dim(stat2)
dim(stat3)
stat <- dplyr::bind_rows(stat1, stat2)
stat <- dplyr::bind_rows(stat, stat3)
dim(stat)
print(round(mean(stat$percentage_wgd),2))
print(round(sd(stat$percentage_wgd),2))
