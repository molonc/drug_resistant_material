
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/'
hvg_genes <- data.table::fread(paste0(input_dir,'SA609_3000_hvg_genes.csv')) %>% as.data.frame()
dim(hvg_genes)
hvg_genes$ens_gene_id <- hvg_genes$bc
meta_genes <- annotables::grch38 %>% 
  dplyr::select(ens_gene_id = ensgene, gene_symbol = symbol)

hvg_genes <- hvg_genes %>% left_join(meta_genes,by='ens_gene_id')
hvg_genes$bc <- NULL
data.table::fwrite(hvg_genes, paste0(input_dir,'SA609_3000_hvg_genes_v2.csv.gz'))
head(hvg_genes)

rm(hvg_genes)
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/'
hvg_genes <- data.table::fread(paste0(input_dir,'SA535_3000_hvg_genes.csv')) %>% as.data.frame()
dim(hvg_genes)
hvg_genes$ens_gene_id <- hvg_genes$bc
hvg_genes <- hvg_genes %>% left_join(meta_genes,by='ens_gene_id')
hvg_genes$bc <- NULL
data.table::fwrite(hvg_genes, paste0(input_dir,'SA535_3000_hvg_genes_v2.csv.gz'))
head(hvg_genes)


rm(hvg_genes)
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/'
hvg_genes <- data.table::fread(paste0(input_dir,'SA1035_3000_hvg_genes.csv')) %>% as.data.frame()
dim(hvg_genes)
hvg_genes$ens_gene_id <- hvg_genes$bc
hvg_genes <- hvg_genes %>% left_join(meta_genes,by='ens_gene_id')
hvg_genes$bc <- NULL
data.table::fwrite(hvg_genes, paste0(input_dir,'SA1035_3000_hvg_genes_v2.csv.gz'))
head(hvg_genes)

