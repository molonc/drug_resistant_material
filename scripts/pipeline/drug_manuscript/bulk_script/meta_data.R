input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_resistance/SA609_kallisto/'
input_dir <- '/home/htran/storage/rnaseq_datasets/bulk_resistance/SA535_kallisto/'
sids <- list.files(input_dir)
sids <- sids[grepl('SA609',sids)]
sids <- sids[grepl('SA535',sids)]


files <- list.files(input_dir)
files <- files[grepl('abundance.tsv',files)]
files <- paste0(input_dir,files)






sample_ids <- c('SA535X6XB03099','SA535X6XB03101','SA535X9XB03617',
                'SA535X9XB03616','SA535X8XB03664','SA535X8XB03548','SA535X5XB02891')



sample_ids <- c('SA604X7XB02089','SA604X9XB02425','SA609X4XB03080',
                'SA609X4XB03083','SA609X5XB03230','SA609X6XB03404','SA609X7XB03505',
                'SA609X7XB03510','SA609X7XB03554')

ids1 <- stringr::str_sub(sample_ids, nchar(sample_ids)-3, nchar(sample_ids))
ids1 <- paste0('SA609-',ids1)
ids1 <- paste0('SA535-',ids1)
ids1 <- ids1[!ids1 %in% sids]
meta_df <- data.frame(sample_id = sample_ids, bulk_output_id=ids1)
meta_df$bulk_alignment_kalliso <- ifelse(meta_df$bulk_output_id %in% sids,'DONE','TO DO')
meta_df$note <- ifelse(meta_df$bulk_alignment_kalliso =='DONE','','3617_R2-201620438 ending with _R2 instead of _R')
meta_df$sequencing_center <- ifelse(meta_df$bulk_alignment_kalliso =='DONE','BRC','GSC')
meta_df$sequencing_center <- 'BRC'
data.table::fwrite(meta_df, paste0(input_dir, 'bulk_info_SA609.csv'))
data.table::fwrite(meta_df, paste0(input_dir, 'bulk_info_SA535.csv'))
