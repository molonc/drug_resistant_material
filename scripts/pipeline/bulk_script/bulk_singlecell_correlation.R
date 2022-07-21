base_dir <- '/home/htran/Projects/farhia_project/drug_resistance/rnaseq/differential_expression/results/'
datatag <- 'SA609'
datatag <- 'SA535'
res_dir <- paste0(base_dir,datatag,'-v6/comps/')

# bulk_input_dir <- '/home/htran/storage/datasets/drug_resistance_RNAseq/SA609_bulk/'
# bulk_base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_bulk/'
bulk_base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_bulk/'


# de_genes_fn <- paste0(res_dir,'ressens11_SA535_UUTTTT_T_UUUUUU_J_logfc_results.csv')
de_genes_fn <- paste0(res_dir,'ressens18_SA535_UXXXX_U_UUUUU_J_logfc_results.csv')

# de_genes_fn <- paste0(res_dir,'pathway9_SA609_UTTTT_R_UUUUU_H_logfc_results.csv')
rm(de_genes)
de_genes <- data.table::fread(de_genes_fn) %>% as.data.frame()
print(dim(de_genes))
head(de_genes)
de_genes <- de_genes %>%
  dplyr::filter(abs(logFC)>0.25 & FDR<0.01 & PValue < 0.05)%>%
  dplyr::rename(singlecell_logFC=logFC)
print(dim(de_genes))
# subtag <- '3617_3616_3664'
subtag <- '3099_3664_3548'
# subtag <- '3505_3510_3554'
de_genes_fn_bulk <- paste0(bulk_base_dir, subtag,'/edgeR_significant_genes_2_Rx_1_UnRx.csv')
de_genes_fn_bulk <- paste0(bulk_base_dir, subtag,'/edgeR_significant_genes_1_UnRx_2_Rx.csv')

de_genes_bulk <- data.table::fread(de_genes_fn_bulk) %>% as.data.frame()
print(dim(de_genes_bulk))
head(de_genes_bulk)
de_genes_bulk <- de_genes_bulk %>%
  dplyr::filter(abs(logFC)>=1) %>%
  dplyr::select(gene_symbol,logFC)%>%
  dplyr::rename(bulk_logFC=logFC)

de_genes <- de_genes %>% inner_join(de_genes_bulk, by=c('gene_symbol'))
dim(de_genes)
save_dir <- paste0(bulk_base_dir, subtag, '/')

cr <- cor(de_genes$singlecell_logFC, de_genes$bulk_logFC, method = "spearman")
pltttitle <- paste0(datatag,': SingleCell-Bulk \n Correlation: ',round(cr,3))
de_genes$DE_gene <- ifelse(de_genes$singlecell_logFC>0,'Up-regulated','Down-regulated')
# stat <- data.frame(logFC_edgeR=markers_ls$logFC,SCTransform_norm_exp=mean_df$exp_mean, 
#                    gene_type=markers_ls$DE_gene,
#                    ensembl_gene_id = markers_ls$ensembl_gene_id,
#                    gene_symbol = markers_ls$gene_symbol,
#                    stringsAsFactors = F)

# de_desc_sc <- 'Res(X10:T) vs. Sen(X9:J)'
# de_desc_bulk <- 'Res(X9:Rx,RxH) - Sen(X8:UnRx)'

de_desc_sc <- 'Res(X8:U) - Sen(X9:J)'
de_desc_bulk <- 'Res(X8:Rx) - Sen(X8,X6:UnRx)'
dim(de_genes)
p <- ggplot(de_genes, aes(singlecell_logFC, bulk_logFC, color=DE_gene, shape=DE_gene)) + 
  geom_point(size=2.3, alpha=1) +
  scale_shape_manual(values=c('Up-regulated'=1, 'Down-regulated'=2)) +
  scale_color_manual(values=c('Up-regulated'='#E69F00', 'Down-regulated'='#56B4E9')) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.3) + 
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) 
p <- p + labs(x=paste0('SingleCell log2FC ',de_desc_sc),y=paste0('Bulk log2FC ',de_desc_bulk), title = pltttitle) + 
  theme(plot.title = element_text(color="black", size=13, hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size=12),
        axis.title = element_text(color="black", size=12),
        legend.title = element_text(color="black", size=12), 
        legend.text = element_text(color="black", size=12))
# p
# saveRDS(p, paste0(output_dir,"edgeR_corr_eval_",de_desc,".rds"))
png(paste0(save_dir,"singlecell_bulk_corr_eval.png"), height = 2*350, width=2*500,res = 2*72)
print(p)
dev.off()
