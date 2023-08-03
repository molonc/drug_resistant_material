


save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/'
genes_df <- data.table::fread(paste0(save_dir, 'pseudotime_genes_modules_Pt4_Pt5_Pt6.csv.gz')) %>% as.data.frame()
dim(genes_df)

View(head(genes_df))
unique(genes_df$gene_type)
genes_df$gene_type <- NULL
genes_df$description <- NULL
series=c("SA501","SA530","SA604","SA609","SA535","SA1035")
pts=paste0("Pt",rep(1:6,1))
names(pts) <- series
unique(genes_df$datatag)
genes_df$patient <- pts[genes_df$patient]
unique(genes_df$patient)

dim(genes_df)
data.table::fwrite(genes_df, paste0(save_dir, 'pseudotime_genes_modules_Pt4_Pt5_Pt6.csv.gz'))
treated_series <- c('SA609','SA535','SA1035')
total <- tibble::tibble()
for(datatag in treated_series){
  df <- data.table::fread(paste0(save_dir, datatag,'_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
}  
datatag <- 'SA609'
dim(df)
View(head(df))
