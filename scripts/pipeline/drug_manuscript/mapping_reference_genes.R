# First get csv files for 4 data tables 

# Get avg, sd per series incis intrans 

results_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/SA535-v6/mapped_genes/'

pancancer_df <- read.csv(paste0(results_dir,'ref_fitness_genes_pct.csv'), check.names=F, stringsAsFactors=F)
head(pancancer_df)
dim(pancancer_df)
unique(pancancer_df$PDX)

depmap_df <- read.csv(paste0(results_dir,'broad_sanger_983_essential_gene_pct.csv'), check.names=F, stringsAsFactors=F)
dim(depmap_df)
depmap_df$desc

cosmic_df <- read.csv(paste0(results_dir,'reference_cosmic_gene_pct.csv'), check.names=F, stringsAsFactors=F)
dim(cosmic_df)


resistance_df <- read.csv(paste0(results_dir,'reference_cisplatin_resistance_gene_pct.csv'), check.names=F, stringsAsFactors=F)
dim(resistance_df)
