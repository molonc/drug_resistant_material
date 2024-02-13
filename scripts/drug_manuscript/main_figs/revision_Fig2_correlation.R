
library(dplyr)
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
save_dir <- base_dir

# Change clonealign summary file
summary_stat <- data.table::fread(paste0(save_dir, 'summary_stat_10x_dlp_subsetSA609.csv')) %>% as.data.frame()
'SA1035X4XB02879' %in% unique(summary_stat$sample_id)
unique(summary_stat$data_type1)
excluded <- c('SA609X5XB03231',
              'SA1035X4XB02879',
              'SA1035X5XB03015',
              'SA1035X5XB03021',
              'SA604X7XB02089')
summary_stat1 <- summary_stat %>%
  dplyr::filter(!sample_id %in% excluded & data_type1=='10x')
# unique(summary_stat$data_type1)

clonealign_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/clonealign_plot/clonealign/temp/'
fns <- list.files(clonealign_dir)    
fns <- fns[grepl('*csv*',fns)]
print(fns)
clonealign_df <- tibble::tibble()
tags <- c('SA1035','SA1035','SA1035','SA604','SA609')
for(i in seq(1:5)){
  f <- fns[i]
  if(!file.exists(paste0(clonealign_dir,f))){
    stop('File does not exist: ')
    print(f)
  }
  
  # First clonealign output
  tmp <- data.table::fread(paste0(clonealign_dir,f)) %>% as.data.frame()
  # tmp$datatag <- stringr::str_sub(f,1,5)
  tmp$datatag <- tags[i]
  if(!'library_id' %in% colnames(tmp)){
    lib <- gsub('.cache/','',tmp$Sample)
    lib <- gsub('/filtered_feature_bc_matrix','',lib)
    tmp$library_id <- lib
    tmp$Sample <- NULL  
  }
  clones_props <- tmp %>%
    dplyr::filter(!clone %in% c('unassigned','Unassigned','None')) %>%
    dplyr::mutate(clone=get_unique_clone_id(clone))%>%
    dplyr::group_by(clone) %>%
    dplyr::summarise(clone_prevalence=n(), freq=n()/dim(tmp)[1])
  
  
  clonealign_df <- clonealign_df %>%
    dplyr::rename(sample_id=id, PDX=datatag) %>%
    dplyr::mutate(data_type="scRNA-seq (10x)", data_type1="10x")
  
  clonealign_stat <- dplyr::bind_rows(clonealign_stat, clones_props)
  clonealign_df <- dplyr::bind_rows(clonealign_df, tmp)
}

clonealign_df <- clonealign_df %>%
  dplyr::filter(!unique_clone %in% c('unassigned','Unassigned','None')) %>%
  dplyr::select(-clone) %>%
  dplyr::rename(clone=unique_clone)
unique(clonealign_df$datatag)
dim(clonealign_df)
colnames(summary_stat1)
summary_stat1$sample_id[1]
colnames(clonealign_df)
clonealign_df$id[1]
unique(summary_stat1$data_type1)        
# [2] "clone_prevalence"
# [3] "freq"            



summary_stat <- data.table::fread(paste0(save_dir, 'summary_stat_10x_dlp_subsetSA609.csv')) %>% as.data.frame()
unique(summary_stat$data_type1)
unique(summary_stat$PDX)
colnames(summary_stat)
summary_stat_10x <- summary_stat %>%
  dplyr::filter(data_type1=='10x')
summary_stat_dlp <- summary_stat %>%
  dplyr::filter(data_type1=='DLP')
sum(unique(summary_stat_10x$sample_id) %in% unique(summary_stat_dlp$sample_id))
dim(summary_stat_10x)
unique(summary_stat_10x$clone)
unique(stat1$clone)
length(unique(summary_stat_10x$sample_id))
length(unique(summary_stat_dlp$sample_id))
length(unique(stat1$sample_id))
ls_obs_samples <- unique(summary_stat_10x$sample_id)

clonealign_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/clonealign_revision/total_alignment/'

fns <- list.files(clonealign_dir)    
fns <- fns[grepl('*csv*',fns)]
print(fns)

clonealign_df <- tibble::tibble()

for(i in seq(1:length(fns))){
  f <- fns[i]
  if(!file.exists(paste0(clonealign_dir,f))){
    stop('File does not exist: ')
    print(f)
  }
  
  # First clonealign output
  tmp <- data.table::fread(paste0(clonealign_dir,f)) 
  # print(dim(tmp))
  # print(colnames(tmp))
  tmp <- tmp %>%
    dplyr::select(-clone, -Barcode) %>%
    dplyr::rename(clone=unique_clone) %>%
    dplyr::rename(sample_id=id, PDX=datatag) %>%
    dplyr::mutate(data_type="scRNA-seq (10x)", data_type1="10x")
  unique(tmp$sample_id)
  if(f=="SA609_clonealign_labels_total.csv.gz"){
    obs_samples <- c('SA609X3XB01584','SA609X4XB03080','SA609X5XB03223',
                     'SA609X6XB03447','SA609X7XB03554','SA609X4XB003083',
                     'SA609X5XB03230','SA609X6XB03404','SA609X7XB03505',
                     'SA609X5XB03231','SA609X6XB03401','SA609X7XB03510')
    tmp <- tmp %>%
      dplyr::filter(sample_id %in% obs_samples)
  }
  clonealign_df <- dplyr::bind_rows(clonealign_df, tmp)
}  
dim(clonealign_df)
clonealign_df <- clonealign_df %>%
  dplyr::filter(!clone %in% c('unassigned','Unassigned','None'))
  
stat1 <- clonealign_df %>%
  dplyr::filter(sample_id %in% ls_obs_samples) %>%
  dplyr::group_by(clone, sample_id, data_type, data_type1, PDX) %>%
  dplyr::summarise(clone_prevalence=n())
dim(stat1)
stat2 <- clonealign_df %>%
  dplyr::filter(sample_id %in% ls_obs_samples) %>%
  dplyr::group_by(sample_id, PDX) %>%
  dplyr::summarise(nb_cells=n())
dim(stat2)
dim(summary_stat_10x)
length(unique(stat1$sample_id))

cols_use <- colnames(summary_stat_dlp)[colnames(summary_stat_dlp) !='library_id']
stat1 <- stat1 %>%
  dplyr::left_join(stat2, by=c('sample_id','PDX')) %>%
  dplyr::mutate(freq=clone_prevalence/nb_cells) %>%
  dplyr::select(all_of(cols_use))
dim(stat1)
colnames(stat1)
colnames(summary_stat_dlp)

summary_stat <- dplyr::bind_rows(summary_stat_dlp, stat1)
# t <- summary_stat_dlp %>%
#   dplyr::filter(PDX=='SA609') %>%
#   dplyr::group_by(clone) %>%
#   dplyr::summarise(nb_cells=sum(clone_prevalence))
# t1 <- summary_stat_10x %>%
#   dplyr::filter(PDX=='SA609')%>%
#   dplyr::group_by(clone) %>%
#   dplyr::summarise(nb_cells=sum(clone_prevalence))
# t2 <- stat1 %>%
#   dplyr::filter(PDX=='SA609')%>%
#   dplyr::group_by(clone) %>%
#   dplyr::summarise(nb_cells=sum(clone_prevalence))
# 
# unique(t$clone)
# unique(t1$clone)
# unique(t2$clone)
# t$nb_cells
# table(t2$clone, t2$nb_cells)
dim(summary_stat)
clonealign_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/clonealign_revision/total_alignment/'
data.table::fwrite(summary_stat, paste0(clonealign_dir, 'summary_clonal_prevalence_10x_dlp_all_series.csv.gz'))
summary_stat <- data.table::fread(paste0(clonealign_dir, 'summary_clonal_prevalence_10x_dlp_all_series.csv.gz'))
colnames(summary_stat)
dim(summary_stat)
