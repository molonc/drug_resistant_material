
suppressPackageStartupMessages({
  library(argparse)
  # library(snpEnrichment)
  library(SingleCellExperiment)
  library(dplyr)
  # library(parallel)
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)
})

parser <- ArgumentParser(description = "Get reads data for each clone")
parser$add_argument('--results_dir', type = 'character', metavar = 'character',
                    help="results_dir")
parser$add_argument('--save_dir', type = 'character', metavar = 'character',
                    help="save_dir")
parser$add_argument('--output_file', type = 'character', default=NULL, metavar = 'FILE',
                    help="output_file file")
parser$add_argument('--input_dir', type = 'character', metavar = 'character',
                    help="input_dir")
parser$add_argument('--datatag', type = 'character', default=NULL, metavar = 'character',
                    help="datatag")
parser$add_argument('--pct_pure', type = 'character', default='0.6', metavar = 'character',
                    help="pct CN pure in each clone")
parser$add_argument('--cn_diff', type = 'integer', default=3,
                    help="CN different within a clone")
# 
# 
args <- parser$parse_args()
# 







# input_dir <- '/home/htran/storage/datasets/drug_resistance_DLP/SA535_segment/reads_clones'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total_v2'
# tenx_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/SA535_cisplatin'
# source("/home/htran/Projects/farhia_project/dlp/align_clones/utils.R")
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/clonealign/custom_segmentCNV/'
# 
# cores_use=3
# genes_use=NULL
# pct_pure=0.6
# cn_diff=3
# 
# input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA535/reads_clones/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/dlp_results/'
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/clonealign/'
# base_name='SA535_cisplatin'
# output_file=NULL
# get_segment_CNV_profiles(input_dir, results_dir, save_dir, output_file=NULL,
#                          base_name, cores_use=3, genes_use=NULL, pct_pure=0.6,  #pct_pure=0.6
#                          cn_diff=3)
# 
# input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA1035/reads_clones/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/dlp_results/'
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/clonealign/'
# base_name='SA1035_cisplatin'
# output_file=NULL
# 
# get_segment_CNV_profiles(input_dir, results_dir, save_dir, output_file=NULL,
#                          base_name, cores_use=3, genes_use=NULL, pct_pure=0.6,  #pct_pure=0.6
#                          cn_diff=3)
# 
# 
# input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/dlp_results/'
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/clonealign/'
# base_name='SA609_cisplatin'
# output_file=NULL
# 
# get_segment_CNV_profiles(input_dir, results_dir, save_dir, output_file=NULL,
#                          base_name, cores_use=3, genes_use=NULL, pct_pure=0.6,  #pct_pure=0.6
#                          cn_diff=3)

# input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA530/reads_clones/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA530_rna/dlp_results/'
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA530_rna/clonealign/'
# base_name='SA530'
# output_file=NULL
# 
# get_segment_CNV_profiles(input_dir, results_dir, save_dir, output_file=NULL,
#                          base_name, cores_use=3, genes_use=NULL, pct_pure=0.6,  #pct_pure=0.6
#                          cn_diff=3)
source("/home/htran/Projects/farhia_project/dlp/align_clones/utils.R")


get_segment_CNV_profiles <- function(input_dir, results_dir, save_dir, output_file=NULL,
                                     base_name='SA', 
                                     cores_use=3, genes_use=NULL, pct_pure=0.6,  #pct_pure=0.6
                                     cn_diff=3){
  if(!dir.exists(save_dir)){
    dir.create(save_dir)  
  }
  clone_fn <- paste0(results_dir, base_name,'_cell_clones.csv')
  if(!file.exists(clone_fn)){
    stop('cell_clones file is not exist, double check input data!!!')
  }
  cell_clones <- data.table::fread(clone_fn) %>% as.data.frame()
  print(dim(cell_clones))
  print(summary(as.factor(cell_clones$clone_id)))
  segments <- data.table::fread(paste0(input_dir,'/',base_name,'_segments_total.csv'))
  # View(head(segments))
  fn <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA609/reads_clones/reads_clone_A.csv'
  segments <- data.table::fread(fn) %>% as.data.frame()
  # dim(segments)  
  
  # summary(segments$pct_pure)
  # sum(segments$pct_pure>=0.6)
  # print(colnames(segments)) ## To Do: add median cn value for each clone, each region, check state is median value?
  # print(length(unique(segments$clone_id))) # cluster id
  # print(summary(as.factor(segments$clone_id)))
  # corrected_sce <- readRDS(file.path(tenx_dir,'SA535_norm_batch_corrected_sce.rds'))
  # genes_use <- rownames(corrected_sce)
  # genes_use[1]
  # dim(corrected_sce)
  
  
  ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
  reference_genes_df <- data.table::fread(paste0(ref_dif,'Symbol_ensembl.csv')) %>% as.data.frame()
  rownames(reference_genes_df) <- reference_genes_df$Ensembl
  if(is.null(genes_use)){
      genes_use <- reference_genes_df$Ensembl
      print(length(genes_use))
  }
  if(!grepl('chr',segments$chr[1])){
    segments <- segments %>%
      dplyr::mutate(chr=paste0("chr", chr))
  }
  unassigned_clone_lbs <- paste0('clone_',c('None','Unassigned', 'Un','un'))
  segments <- segments %>%
    # dplyr::filter(clone_id =='clone_D') %>%
    # dplyr::filter(!clone_id %in% unassigned_clone_lbs) %>%
    dplyr::rename(cluster=clone_id, copy_number=state)%>%
    dplyr::select(chr,start,end,cluster,copy_number,pct_pure)
  
  segments$pct_pure <- round(segments$pct_pure, 2)
  # head(segments)
  # test <- annotables::grch38 %>% 
  #   # dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
  #   dplyr::select(entrez)
  # length(unique(test$entrez))  
  # print(dim(segments))
  # print(unique(segments$cluster))
  dim(segments)
  gene_cn <- get_gene_coordinates(segments)
  # head(gene_cn)
  # gene_cn <- gene_cn %>% 
  #   dplyr::filter(ensembl_gene_id %in% genes_use)
  print(dim(gene_cn))
  # sum(is.na(gene_cn$percent_overlap))
  unique_genes <- unique(gene_cn$ensembl_gene_id)
  print(length(unique_genes))
  
  
  excluded_genes <- reference_genes_df %>%
    dplyr::filter(!Ensembl %in% unique_genes) %>%
    dplyr::pull()
  
  output_genecn <- paste0(save_dir,'gene_cn_summarized.feather')
  start_time <- Sys.time()
  gene_cn_summarized <- get_gene_cn_summary_custom(unique_genes, base_name, save_dir, 
                                                   gene_cn, output_genecn, cores_use)
  
  end_time <- Sys.time()
  
  print(end_time - start_time)
  
  print(dim(gene_cn_summarized))
  View(head(gene_cn_summarized))
  # gene_cn_fn <- output_file
  # # output_file <- paste0(save_dir,'cnv.csv')
  # gene_cn_summarized <- feather::read_feather(gene_cn_fn)
  # colnames(gene_cn_summarized)
  
  # gene_cn_summarized <- gene_cn_summarized %>% 
  #   dplyr::filter(!cluster %in% c('None','Unassigned'))
  
  print(dim(gene_cn_summarized))
  
  # View(gene_cn_summarized[gene_cn_summarized$ensembl_gene_id==bad_genes$ensembl_gene_id[1],])
  # bad_genes <- gene_cn_summarized %>%
  #   dplyr::filter((copy_number_max - copy_number_min) > cn_diff)
  
  
  # gene_cn_summarized <- gene_cn_summarized %>%
  #   dplyr::filter(!ensembl_gene_id %in% bad_genes$ensembl_gene_id)
  # print(dim(gene_cn_summarized))
  # print(paste0("Length of bad genes: ",length(bad_genes$ensembl_gene_id)))
  
  # bad_genes <- bad_genes %>% as.data.frame()
  # bad_genes$gene_symbol <- reference_genes_df[bad_genes$ensembl_gene_id,'Symbol']
  # # View(bad_genes)
  # data.table::fwrite(bad_genes, paste0(save_dir,'bad_genes.csv'))
  # reference_genes_df$bad_genes <- 'F'
  # reference_genes_df[bad_genes$ensembl_gene_id,'bad_genes'] <- 'T'

  gene_cn_summarized <- gene_cn_summarized %>%  # there is no time here
    dplyr::group_by(cluster, ensembl_gene_id) %>%
    dplyr::summarise(mean_cnmean=mean(copy_number_mean, na.rm = T),
                     mean_cnmode=mean(copy_number_mode, na.rm = T),
                     median_cnmode=median(copy_number_mode, na.rm = T),
                     mode_cnmode=calc_mode(copy_number_mode, na.rm = T),
                     pct_pure=sum(copy_number_mode == mode_cnmode & copy_number_min == copy_number_max)/n(),
                     n_cell_time=n())%>%
    ungroup()
  
  
  # cluster_variability <- gene_cn_summarized %>%
  #   dplyr::mutate(l_median_cnmode=log(median_cnmode+1)) %>%
  #   dplyr::group_by(ensembl_gene_id) %>%
  #   dplyr::summarise(mad=mad(l_median_cnmode),
  #                    sd=sd(l_median_cnmode),
  #                    max=max(median_cnmode),
  #                    min=min(median_cnmode),
  #                    range=max-min,
  #                    range_log=max(l_median_cnmode)-min(l_median_cnmode))
  # 
  # print(paste0("cluster_variability: ",dim(cluster_variability)))
  
  # impure_genes <- unique(cluster_time_cn_summarized$ensembl_gene_id[cluster_time_cn_summarized$pct_pure < args$pct_pure])
  print(summary(gene_cn_summarized$pct_pure))
  print(sum(gene_cn_summarized$pct_pure <= pct_pure))
  impure_genes <- unique(gene_cn_summarized$ensembl_gene_id[gene_cn_summarized$pct_pure < pct_pure])
  # sum(gene_cn_summarized$pct_pure >= pct_pure)
  # intertimepoint_variable_genes <- unique(timepoint_variability$ensembl_gene_id[timepoint_variability$diff > 1])
  # exclude_genes <- union(impure_genes, intertimepoint_variable_genes)
  
  exclude_genes <- impure_genes
  exclude_genes_df <- data.frame(ensembl_gene_id=exclude_genes, stringsAsFactors=F)
  exclude_genes_df$gene_symbol <- reference_genes_df[exclude_genes_df$ensembl_gene_id,'Symbol']
  # View(bad_genes)
  data.table::fwrite(exclude_genes_df, paste0(save_dir,'impure_genes.csv'))
  reference_genes_df$impure_genes <- F
  reference_genes_df[exclude_genes_df$ensembl_gene_id,'impure_genes'] <- T
  data.table::fwrite(reference_genes_df, file=paste0(save_dir, base_name, '_genes_infos.csv'))
  # View(exclude_genes_df)
  print(paste0("Nb exclude_genes: ",length(exclude_genes)))
  # selected_cn <- cluster_variability %>% 
  #   dplyr::filter(range > 0,
  #                 !ensembl_gene_id %in% exclude_genes) %>% 
  #   dplyr::arrange(-mad)
  
  
  
  # gene_cn_summarized <- gene_cn_summarized %>%
  #     dplyr::mutate(use_gene = !ensembl_gene_id %in% exclude_genes)
  # gene_cn_summarized <- gene_cn_summarized %>%
  #   dplyr::mutate(use_gene = ensembl_gene_id %in% selected_cn$ensembl_gene_id)
  # write.table(cluster_cn_summarized, args$outfname, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  # print(summary(gene_cn_summarized$use_gene))
  # write.table(cluster_cn_summarized, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  # saveRDS(gene_cn_summarized,paste0(save_dir, base_name,'.rds'))
  output_cnv <- paste0(save_dir, base_name, '_cnv.csv')
  data.table::fwrite(gene_cn_summarized, file=output_cnv)
  
  # paste0(save_dir,'cluster_cn_summarized.txt'),
  print(paste0("Save output: cluster_cn_summarized: ",
               dim(gene_cn_summarized)[1],' ',dim(gene_cn_summarized)[2]))
  
  print("Completed")
  
  # Get segment CNV profile
  # raw_cnvs <- gene_cn_summarized
  # cnv <- dplyr::filter(raw_cnvs, use_gene) %>%
  cnv <- gene_cn_summarized %>%
    dplyr::rename(clone = cluster,
                  copy_number=median_cnmode) %>% 
    dplyr::select(ensembl_gene_id, clone, copy_number) %>% 
    # dplyr::filter(clone %in% present_clones) %>%
    spread(clone, copy_number)
  
  cnv_mat <- cnv %>%
    as.data.frame %>%
    column_to_rownames("ensembl_gene_id") 
  # %>%as.matrix
  
  print("Debug: cnv_mat ")
  print(colnames(cnv_mat))
  if(grepl('clone_',colnames(cnv_mat)[1])){
    colnames(cnv_mat) <- gsub('clone_','',colnames(cnv_mat))
  }
  
  # View(head(cnv))
  # print(dim(cnv_mat))
  # t <- rowVars(as.matrix(cnv_mat))
  if(is.null(output_file)){
    output_file <- paste0(save_dir,base_name,'_cnv_mat.csv')
  }  
  
  data.table::fwrite(cnv, paste0(save_dir,base_name,'_cnv.csv'), quote=F)
  data.table::fwrite(cnv_mat, output_file, row.names = T, quote=F)
  # View(head(cnv))
  print("Completed.\n")
  # View(head(cnv_mat))
}





# de_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/SA535-v6/SA535_UUTTTT_T_UUUUUU_J/'
# de_genes <- data.table::fread(paste0(de_dir,'signif_genes.csv'))
# summary(as.factor(de_genes$classified_gene_dlp))
# # Check out
# cnv_mat <- cnv %>%
#   as.data.frame %>%
#   column_to_rownames("ensembl_gene_id")
# cnv_mat <- cnv_mat[,c('T','J')]
# head(cnv_mat)
# dim(cnv_mat)
# t <- rowVars(as.matrix(cnv_mat))
# cnv_mat <- cnv_mat[rowVars(as.matrix(cnv_mat))>0,]
# length(t)
# sum(de_genes$ensembl_gene_id %in% rownames(cnv_mat))

# input_dir <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA501/reads_clones/'
# results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA501_rna/dlp_results/'
# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA501_rna/clonealign/'
# base_name='SA501'
# output_file=NULL

get_segment_CNV_profiles(args$input_dir, args$results_dir, args$save_dir, args$output_file,
                         args$datatag, cores_use=3, genes_use=NULL, as.numeric(args$pct_pure),  #pct_pure=0.6
                         args$cn_diff)



