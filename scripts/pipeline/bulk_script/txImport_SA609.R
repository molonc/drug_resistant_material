
suppressMessages({
  # library("sleuth")
  library("dplyr")
  library("tximport")
  # library("tximportData")
  # library("TxDb.Hsapiens.UCSC.hg19.knownGene")
  library("edgeR")
  library("csaw")
})
input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_resistance/SA609_kallisto/'


meta_df <- data.table::fread(paste0(input_dir,'bulk_info_SA609.csv'))
meta_df <- meta_df %>%
  dplyr::filter(bulk_alignment_kalliso=='DONE')
dim(meta_df)

txi <- load_kallisto_data(meta_df)


