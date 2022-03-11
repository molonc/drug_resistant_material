suppressPackageStartupMessages({
  require(DOSE)
  require(enrichplot)
  require(ggraph)
  require(ggplot2)
  require(igraph)
})


input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609_rna/deg_analysis/'

  
generate_genes_network(paste0(input_dir,'SA609_UTTTT_UUUUUUUU_R_E/'), 
                       datatag='SA609_R-UTTTT_vs_E-UUUUUUUU', 
                       de_desc='SA609_UTTTT_UUUUUUUU_R_E',
                       genes_set='hallmark',
                       pattern_use='^hallmark_'
)

generate_genes_network(paste0(input_dir,'SA609_UTT_UUUUUUUU_R_H/'), 
                       datatag='SA609_R-UTT_vs_H-UUUUUUUU', 
                       de_desc='SA609_UTT_UUUUUUUU_R_H',
                       genes_set='hallmark',
                       pattern_use='^hallmark_'
)

generate_genes_network(paste0(input_dir,'SA609_UTTTT_UUUUUUUU_R_H/'), 
                       datatag='SA609_R-UTTTT_vs_H-UUUUUUUU', 
                       de_desc='SA609_UTTTT_UUUUUUUU_R_H',
                       genes_set='hallmark',
                       pattern_use='^hallmark_'
)

generate_genes_network(paste0(input_dir,'SA609_R_UTTT_H_UUUUUUUU/'), 
                       datatag='SA609_R-UTTT_vs_H-UUUUUUUU', 
                       de_desc='SA609_R_UTTT_H_UUUUUUUU',
                       genes_set='hallmark',
                       pattern_use='^hallmark_'
)


generate_genes_network(input_dir, datatag='SA609_R-UTTTT_vs_E-UUUUUUUU', 
                       obs_pw='top_up_pathway_cls_',
                       de_desc='SA609_UTTTT_UUUUUUUU_R_E',
                       genes_set='hallmark',
                       pattern_use='^hallmark_'
)


input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/deg_analysis/'

generate_genes_network(paste0(input_dir,'SA535_A_B/'), 
                       datatag='SA535_A_vs_B', 
                       de_desc='SA535_A_B',
                       genes_set='hallmark',
                       pattern_use='^hallmark_'
)


generate_genes_network(paste0(input_dir,'SA535_A_G/'), 
                       datatag='SA535_A_vs_G', 
                       de_desc='SA535_A_G',
                       genes_set='hallmark',
                       pattern_use='^hallmark_'
)

generate_genes_network(paste0(input_dir,'SA535_DE_B/'), 
                       datatag='SA535_D-E_vs_B', 
                       de_desc='SA535_DE_B',
                       genes_set='hallmark',
                       pattern_use='^hallmark_'
)

generate_genes_network(paste0(input_dir,'SA535_DE_G/'), 
                       datatag='SA535_D-E_vs_G', 
                       de_desc='SA535_DE_G',
                       genes_set='hallmark',
                       pattern_use='^hallmark_'
)


# genes <- colnames(pathway_df)
# for(p in rownames(pathway_df)){
#   p_tmp <- c()
#   for(g in genes){
#     if(!is.na(pathway_df[p,g])){
#       p_tmp <- c(p_tmp, g)
#     }
#   }
#   p <- gsub(pattern_use,'',tolower(p))
#   geneSets[[p]] <- p_tmp
# }