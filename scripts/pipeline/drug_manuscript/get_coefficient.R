
# install.packages("readxl")
library("readxl")
library(dplyr)
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/DLP_trees/'
df <- read_excel(paste0(base_dir, 'clone_fractions.xlsx'))
df <- df %>% 
  dplyr::filter(ncells>=50)%>% 
  dplyr::mutate(mean_1_plus_s=round(mean_1_plus_s,3),
                sd_1_plus_s=round(sd_1_plus_s,3),
                median_1_plus_s=round(median_1_plus_s, 3),
                fraction=round(fraction,2))
dim(df)
head(df)
colnames(df)
sa1035 <- df %>% 
  dplyr::filter(grepl('SA1035',datasetname) & clone_name %in% c('H','E'))
sa1035$patient_id <- 'Pt6'
dim(sa1035)

sa535 <- df %>% 
  dplyr::filter(grepl('SA535',datasetname) & clone_name %in% c('A','G'))
sa535$patient_id <- 'Pt5'
dim(sa535)
sa609 <- df %>% 
  dplyr::filter(grepl('SA609',datasetname) & clone_name %in% c('A','H'))
dim(sa609)
View(sa609)
sa609$patient_id <- 'Pt4'

stat <- dplyr::bind_rows(list(sa609, sa535, sa1035))
data.table::fwrite(stat, paste0(base_dir,'clone_coefficients.csv'))
