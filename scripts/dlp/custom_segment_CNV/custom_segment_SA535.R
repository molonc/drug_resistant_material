# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  # library(scater)
  library(data.table)
  library(methods)
  # library(scran)
  library(parallel)
  # library(feather)
  # library(annotables)
  library(dplyr)
  library(ggplot2)
  # library(snpEnrichment)
  # library(scMerge)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene) # old version
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(IRanges)
})
source('/home/htran/Projects/farhia_project/dlp/align_clones/custom_segment_CNV/utils_custom_mapping.R')
# initial.options <- commandArgs(trailingOnly = FALSE)
# file.arg.name <- "--file="
# script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
# script.basename <- dirname(script.name)
# # print(script.basename)
# source(paste0(script.basename, "/utils.R"))
# 
# 
# parser <- ArgumentParser(description = "Get reads data for each clone")
# parser$add_argument('--output_file', type = 'character', metavar = 'character',
#                     help="Output path for gene CN table (feather file).")
# parser$add_argument('--cellclones', type = 'character', default=NULL, metavar = 'FILE',
#                     help="cellclones file")
# parser$add_argument('--obs_clones', type = 'character', default=NULL, metavar = 'character',
#                     help="obs_clones")
# parser$add_argument('--input_dir', type = 'character', default=NULL, metavar = 'character',
#                     help="input_dir")
# parser$add_argument('--datatag', type = 'character', default=NULL, metavar = 'character',
#                     help="datatag")
# parser$add_argument('--library_grouping', type = 'character', default=NULL, default='library_groupings.csv', metavar = 'FILE',
#                     help="library_groupings.csv file")
# 
# 
# args <- parser$parse_args()
# # 
# output_file <- args$output_file
# print(output_file)
# 
# print(args$library_grouping)
# print(args$cellclones)
# # print(args$input_dir)
# if(!file.exists(args$segment_file) | !file.exists(args$library_grouping) 
#    | !file.exists(args$cellclones)){
#   stop('Input files do not exist, pls double check')
# }
# cellclones_fn <- args$cellclones
# input_dir <- args$input_dir
# obs_clones <- args$obs_clones

# ex: obs_clones <- c('I,'B-D') 
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
  
  cell_clones <- read.delim(cellclones_fn, stringsAsFactors = FALSE, sep = ",")
  print("DEBUG 1")
  print(dim(cell_clones))
  print("Remove unassigned cells")
  cell_clones <- cell_clones[!cell_clones$clone_id %in% c('Un','unassigned','Unassigned','None','un'),]
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  # unique(cell_clones$clone_id)
  # ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  # reference_genes_df <- data.table::fread(paste0(ref_dif,'Symbol_ensembl.csv')) %>% as.data.frame()
  # reference_genes_df <- reference_genes_df %>%
  #   dplyr::rename(ensembl_gene_id=Ensembl)
  # obs_clones <- 'T'
  chrs <- c(as.character(1:22), "X")  # Remove Y from analysis
  r38 <- annotables::grch38 %>% 
    dplyr::select(ensembl_gene_id = ensgene, gene_symbol = symbol, chr) %>%
    dplyr::filter(chr %in% chrs) %>% 
    dplyr::select(-chr)
  # dim(r38)
  
  cls <- unique(cell_clones$clone_id)
  if(is.null(obs_clones)){
    print("Get reads data for each clone")
    obs_clones <- cls
  }else{
    obs_clones <- union(obs_clones, cls) # take into account all clones in datasets. 
  }
  if(is.null(tag)){
    tag <- paste(obs_clones, sep='_')
  }  
  print(obs_clones)
  # cols_use <- c('chr', 'start', 'end', 'width', 'gc', 'map', 'reads','cell_id')
  total_reads <- tibble::tibble()
  mapped_clones <- tibble::tibble()
  for(c in obs_clones){
    cs <- as.character(unlist(strsplit(c, '-')))
    if(length(cs)!=sum(cs %in% cls)){
      print('Double check input data!!!')
    }
    cs <- cs[cs %in% cls]
    
    cells_id <- cell_clones %>%
      dplyr::filter(clone_id %in% cs) %>%
      dplyr::pull(cell_id)
    print('__________________________________')
    print(paste0('Observed clone: ',c,' and nb cells: ', length(cells_id)))
    obs_libs <- get_library_id(cells_id)
    obs_libs <- unique(obs_libs)
    print(paste0('Clone: ',c,' and list of observed libraries: '))
    print(paste(obs_libs, collapse = ' '))
    obs_reads <- list()
    for(l in obs_libs){
      lib_fn <- paste0(input_dir,l,'/hmmcopy/',l,'_reads.csv.gz') # csv.gz or csv
      if(file.exists(lib_fn)){
        tmp <- data.table::fread(lib_fn)
        # View(head(tmp))
        # dim(tmp)
        tmp <- tmp %>%
          dplyr::filter(cell_id %in% cells_id)
        tmp <- tmp %>%
          dplyr::select(chr, start, end, width, gc, map, reads, cell_id, state) # 
        print(l)
        print(dim(tmp))
        if(dim(tmp)[1]>0){
          obs_reads[[l]] <- tmp  
        }else{
          print(paste0('***Warning:  Do not have any output for clone: ',c))
        }
        
      }else{
        print(paste0('***ERROR:  Clone: ',c, ' and library: ', l))
        # print('Do not exist hmmcopy reads for library, double check input data')
        stop('Do not exist hmmcopy reads for library, double check input data')
      }
      
    }
    reads_df <- as.data.frame(dplyr::bind_rows(obs_reads))
    print(dim(reads_df))
    # data.table::fwrite(reads_df,file=paste0(output_dir,'clone_',c,'.csv'), row.names=F, quote=F, sep=',')
    ## Version 1
    # stat_reads_df <- reads_df %>% 
    #   dplyr::group_by(chr, start, end) %>% 
    #   dplyr::summarise(reads=sum(reads),
    #                    gc=mean(gc),
    #                    map=mean(map))%>% 
    #   dplyr::ungroup()%>%
    #   dplyr::select(chr, start, end, gc, map, reads)

    ## Version 2
    stat_reads_df <- reads_df %>% 
      dplyr::group_by(chr, start, end) %>% 
      dplyr::summarise(reads=sum(reads, na.rm=T),
                       gc=round(mean(gc, na.rm=T), 4),
                       avg_mappability=round(mean(map, na.rm=T), 4),
                       cn_mean=round(mean(state, na.rm=T), 4),
                       cn_min=min(state, na.rm=T),
                       cn_max=max(state, na.rm=T),
                       cn_median=median(state, na.rm=T),
                       cn_mode=calc_mode(state),
                       pct_pure=round(sum(cn_mode == state)/n(), 4),
                       n_cells=n()
                       )%>% 
      dplyr::ungroup()#%>%
      # dplyr::select(chr, start, end, cn_median, cn_mode, avg_mappability, cn_mean, cn_min, cn_max, reads, gc)
    
    print(dim(stat_reads_df))
    # print(summary(stat_reads_df$pct_pure))
    # print(sum(stat_reads_df$pct_pure>=0.6))
    # print(sum(stat_reads_df$avg_mappability<0.9))
    # colnames(stat_reads_df)
    # stat_reads_df <- stat_reads_df %>%
    #   dplyr::select(chr, start, end, width, gc, map, reads)
    
    # datatag <- args$datatag #'SA535'
    stat_reads_df$clone_id <- c
    data.table::fwrite(stat_reads_df,file=paste0(output_dir,'reads_clone_',c,'.csv'))
    total_reads <- dplyr::bind_rows(total_reads, stat_reads_df)
    gene_cn <- get_mapping_genes(stat_reads_df)
    print(dim(gene_cn))
    data.table::fwrite(gene_cn,file=paste0(output_dir,'mapped_gene_clone_',c,'.csv.gz'))
    gene_cn <- gene_cn %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::summarise(avg_mappability=mean(avg_mappability),
                       cn_median=mean(cn_median),
                       pct_pure=mean(pct_pure)) 
    gene_cn$clone_id <- c
    mapped_clones <- dplyr::bind_rows(mapped_clones, gene_cn)
    gene_cn2 <- gene_cn %>%
      dplyr::select(-clone_id) %>%
      dplyr::rename(!!sym(paste0('avg_mappability_',c)):=avg_mappability,
                    !!sym(paste0('cn_median_',c)):=cn_median,
                    !!sym(paste0('pct_pure_',c)):=pct_pure)
    # dim(gene_cn2)
    # data.table::fwrite(gene_cn2,file=paste0(output_dir,'filtered_mapped_gene_clone_',c,'.csv.gz'))
    # r38 <- merge(r38, gene_cn2, by='ensembl_gene_id', all=T)
    r38 <- r38 %>% left_join(gene_cn2, by='ensembl_gene_id')
    # dim(r38)
  }  
  
  # summary <- data.frame(obs_clone=obs_clones)
  # write.csv(summary, file=output_file) # run keep it here to run snakemake program
  # data.table::fwrite(total_reads, file=output_file)
  data.table::fwrite(total_reads,file=paste0(output_dir,'total_reads_',tag,'.csv.gz'))
  data.table::fwrite(mapped_clones,file=paste0(output_dir,'mapped_clones_',tag,'.csv.gz'))
  r38 <- r38 %>% tidyr::drop_na()
  data.table::fwrite(r38,file=paste0(output_dir,'mapped_wholedata_',tag,'_grch38.csv.gz'))
  
  # mapped_clones <- data.table::fread('/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones_v2/mapped_wholedata_SA609_v2.csv.gz') %>% as.data.frame()
  # dim(mapped_clones)
  # head(mapped_clones)
  mapped_clones <- mapped_clones %>%
    dplyr::select(-avg_mappability, -pct_pure) %>%
    tidyr::pivot_wider(names_from='clone_id', values_from='cn_median')
  
  data.table::fwrite(mapped_clones,file=paste0(output_dir,'mapped_wholedata_',tag,'_v2.csv.gz'))

}

# get_reads_per_clone(input_dir, cellclones_fn, output_file, obs_clones, cores_use=NULL)


# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total_v2/'
# base_name <- 'SA535_total'
# input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA535_segment/'
# output_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA535_segment/reads_clones/testing2/'
# output_file <- paste0(output_dir,'out.csv')
# cell_clones <- read.delim(paste0(results_dir,'cell_clones.csv'), stringsAsFactors = FALSE, sep = ",")
# obs_clones <- 'T'

cellclones_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cell_clones/SA535_cell_clones.csv'
input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA535_segment/'
output_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA535_segment/reads_clones_v2/'
output_file <- paste0(output_dir,'out.csv')
# obs_clones <- c('A','G','A-E-J','D')
obs_clones <- c("A","G","A-E-J","D", "B","E", "F", "H", "I", "J", "C")
obs_clones <- unique(obs_clones)
tag='SA535'
get_reads_per_clone(input_dir, cellclones_fn, output_file, obs_clones, tag, cores_use=NULL)


# CNbins <- read.csv(paste0(output_dir,'segment_cloneT/0/segs.csv'), stringsAsFactors = F, check.names = F)
# dim(CNbins)
# head(CNbins)
# /home/htran/anaconda3/envs/sisyphus/bin/python /home/htran/Projects/farhia_project/rscript/dlp/corrupt_tree/src/download_data-scgenome.py --library_id A96244B --download_dir /home/htran/storage/datasets/drug_resistance_DLP/testing/A96244B

# python /home/htran/Projects/farhia_project/rscript/dlp/custom_segment/hmmcopy/scripts/correct_read_count_v2.py 
# --reads /home/htran/storage/datasets/drug_resistance_DLP/SA535_segment/reads_clones/SA535_reads_clone_T.csv --output /home/htran/storage/datasets/drug_resistance_DLP/SA535_segment/reads_clones/corrected_SA535_reads_clone_T.csv


# output_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/SA535_custom_segment/'
# cellclones_fn <- '/home/htran/storage/datasets/metastasis_results/SA535_new_encoding/SA535_wholedata_v2/cell_clones.csv'
# input_dir <- '/home/htran/storage/datasets/hakwoo_metastasis/SA535/'
# datatag <- 'SA535'