

# Generate evaluation plots


# load raw data - pre-filtered genes 
# Load HK genes, plot HK genes variance 
# Load scSEG genes, plot sc genes variance 
# plot PCA




# plttitle <- 'HK genes'
# p_hk <- plot_pca(sce_raw, res$ref_stable_genes$meta_genes_HK, meta_cells, save_dir, plttitle, 20)
# plttitle <- 'ref scMerge scSEG genes'
# p_sc <- plot_pca(sce_raw, res$ref_stable_genes$meta_genes_scSEG, meta_cells, save_dir, plttitle, npcs=20)
# p_total_pca <- cowplot::plot_grid(p_hk, p_sc, align = 'h')
# png(paste0(save_dir,"HK_genes_versus_scSEG_PC1_PC2_rawdata.png"), height = 2*430, width=2*1200,res = 2*72)
# print(p_total_pca)
# dev.off()



# Rscript compute_scSEG_stable_genes.R 
# -i /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA535-v6/SA535_CX5461_total_sce_v3.rds 
# -o /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/SA535_CX5461/segIndx_total.csv

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA535'
tag <- 'CX5461'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag, '_',tag,'/')
dir.create(save_dir)





input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
datatag <- 'SA609'
save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/',datatag,'/')
dir.create(save_dir)
# Rscript compute_scSEG_stable_genes.R 
# -i /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/total_sce_v3.rds 
# -o /home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/SA609/segIndx_total.csv
