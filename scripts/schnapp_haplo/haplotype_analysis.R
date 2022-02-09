# https://www.biostars.org/p/254848/
# https://www.science.gov/topicpages/a/allele-specific+rna+transcription
# https://www.ncbi.nlm.nih.gov/pubmed/19768609
# In this context, the “B” allele is the non-reference allele observed 
# in a germline heterozygous SNP, i.e. in the normal/control sample. 
# Since the tumor cells’ DNA originally derived from normal cells’ DNA, 
# most of these SNPs will also be present in the tumor sample. 
# But due to allele-specific copy number alterations, loss of heterozygosity 
# or allelic imbalance, the allelic frequency of these SNPs may be different
# in the tumor, and that’s evidence that one (or both) of the germline copies 
# was gained or lost during tumor evolution. The shift in b-allele frequency
# is calculated relative to the expected heterozygous frequency 0.5, and 
# minor allele frequencies are “mirrored” above and below 0.5 so that it does 
# not matter which allele is considered the reference – the relative shift 
# from 0.5 will be the same either way. 
# (Multiple alternate alleles are not considered here.)

# The B-Allele Frequency is a normalized measure of the allelic intensity ratio
# of two alleles (A and B), such that a BAF of 1 or 0 indicates the complete 
# absence of one of the two alleles (e.g. AA or BB), 
# and a BAF of 0.5 indicates the equal presence of both alleles (e.g. AB).

# It is a usefull measure when you are studying CNVs (LOHs or also SVs):
# " detection of allelic imbalances such as those caused by duplications (e.g. AAB/BBA) 
# or mosaic deletions in the sample. Such imbalances can be identified on a BAF plot by 
# the presence of SNPs at frequencies between 0.5 and 0 or 1. For example, the theoretical
# BAF values of triploid regions (AAA, AAB, ABB or BBB) are 0, 0.33, 0.66 and 1 respectively."

# read libids, get total allele_counts


require(data.table)
source_dir <- ''
load(paste0(source_dir,'schnapp_utils.R'))

download_dir = '/home/htran/storage/datasets/hakwoo_metastasis/'
save_dir <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/schnapps_out/"

if (!file.exists(save_dir)){
  dir.create(save_dir)
}

# library_ids = c('A96204B','A96200A','A96229A','A96117B','A96117A',
#                 'A98300B','A98290A','A98299A',
#                'A98232B','A98306A','A98261B',
#                'A98217B','A98217A','A98293A')
meta_sample <- read.csv(paste0(download_dir,'SA919X7_metadata_Hoa.csv'),header=T,check.names = F)
# View(head(meta_sample))
dim(meta_sample)
meta_sample <- meta_sample[meta_sample$DLP_analysis=='TRUE',]
# sum(library_ids %in% meta_sample$`Library ID_DLP+`)
library_ids <- meta_sample$`Library ID_DLP+`
prefix <- 'total'



# Read filtered file 
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'


total_filtered_cn <- read.delim(paste0(results_dir,'total_merged_filtered_states_original.csv'),
                                check.names = FALSE, stringsAsFactors = FALSE, sep = ",")


dim(total_filtered_cn)
print(rownames(total_filtered_cn)[1:3])
cells_use <- colnames(total_filtered_cn)
length(cells_use)
genes_use <- rownames(total_filtered_cn)
length(genes_use)
# get_allele_counts(cells_use, library_ids, download_dir, save_dir, prefix)

get_hmmcopy_reads(cells_use, genes_use, library_ids, download_dir, save_dir, prefix)












visualize_tsne_clustering <- function(){
  clusters_dbscan <- readRDS(paste0(save_dir,prefix,'_clusters_hdbscan.rds'))
  View(head(clusters_dbscan))
  rownames(clusters_dbscan) <- clusters_dbscan$cell_id
  length(clones_labels_df$cell_id)
  cells_use <- unique(CNbins$cell_id)
  meta_data <- data.frame(ct_clone=clones_labels_df[cells_use,'clone_id'],
                          hdbscan_clone=clusters_dbscan[cells_use,'clone_id'],
                          row.names = cells_use)
  # pcs <- prcomp(t(total_filtered_cn[,cells_use]), scale. = TRUE)
  library(Seurat)
  srt_cnv <- CreateSeuratObject(counts = total_filtered_cn[,cells_use], meta.data = meta_data)
  srt_cnv <- FindVariableFeatures(object = srt_cnv, verbose = FALSE)
  srt_cnv <- ScaleData(object = srt_cnv, verbose = FALSE)
  srt_cnv <- RunPCA(object = srt_cnv,verbose = FALSE,npcs = 50)
  srt_cnv <- FindNeighbors(srt_cnv, dims = 1:20, verbose = FALSE)
  srt_cnv <- RunTSNE(object = srt_cnv, verbose = FALSE,reduction = "pca",dims = 1:50)
  p1 <- DimPlot(srt_cnv, reduction = "tsne", group.by = 'ct_clone')
  p2 <- DimPlot(srt_cnv, reduction = "tsne", group.by = 'hdbscan_clone')
  p <- CombinePlots(list(p1, p2),nrow=1)
  png(paste0(save_dir,prefix,"_tsne_clones_labels.png"), height = 2*400, width=2*900,res = 2*72)
  print(p)
  dev.off()
  
}




