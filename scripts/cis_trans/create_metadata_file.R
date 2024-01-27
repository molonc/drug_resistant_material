input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/comparisons/'
df <- data.table::fread(paste0(input_dir,'input_files.csv'), header=T) %>% as.data.frame()
df2 <- data.table::fread(paste0(input_dir,'comparisons_drug_res.csv'), header=T) %>% as.data.frame()
View(df)
df$file_header
colnames(df2)
fns <- list.files(data_dir)
fns <- fns[grepl('.csv', fns)]
fns1 <- fns[!fns %in% df2$result_fn]
sum(fns %in% df2$result_fn)

df3 <- data.frame(result_fn=fns1)
pairs <- dplyr::bind_rows(df2, df3)
View(pairs)

pairs$file_header <- sapply(strsplit(pairs$result_fn, '_'), function(x){
  return(paste0(x[1],'_',x[2],'_',x[3]))
})

pairs$file_header <- ifelse(pairs$file_header=='scran-nocl-SA1035-1_SA1035_UT','scran-nocl-SA1035-1', pairs$file_header)
pairs <- pairs %>% left_join(df, by='file_header')
dim(pairs)
View(pairs)
pairs$V1 <- NULL
data.table::fwrite(pairs, paste0(input_dir,'comparisons_drug_res_v2.csv'))



p <- ggplot(df, aes(x=Gene_Type, y=plt_desc)) + 
  geom_point(aes(color = Gene_Type, size=sz), alpha=1) + #, size=2*log10(pct_genes), size=4.5
  scale_color_manual(values = keyvals_colour, guide = "none")

p <- p + facet_grid(series ~ ., scales="free_y", space='free')
p <- p + labs(title='Gene type', x=NULL, y=NULL) +
  # theme_bw(base_size = 10) + 
ct_theme <-   
  theme(#legend.position = "bottom",
    legend.position = "none",
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.text  = element_blank(),
    text = element_text(size = 11, hjust = 0.5, family=my_font),
    # axis.text.x = element_text(size=11, hjust = 0.5, family=my_font, angle = 90),  #
    axis.text.x = element_blank(),
    # axis.text.x = element_text(size=11, hjust = 0.5, family=my_font, angle = 90, color='black'), ## for testing
    # axis.text.y = element_text(size=11, hjust = 0.5, family=my_font), ## for testing
    axis.text.y = element_blank(),
    plot.title = element_text(size=11, hjust=0.5, family=my_font), #, face="bold"
    axis.title.y = element_blank(),
    # axis.title.x = element_text(size=11, hjust = 0.5, family=my_font)
    # axis.title.x = element_blank()
    # strip.text.x = element_text(color="black",size=11, family=my_font),
    # strip.text.y = element_text(color="black",size=11, family=my_font),
    strip.text = element_text(size=11, color="black", family=my_font),
    strip.background = element_blank(),
    # strip.placement = "outside",
    panel.background = element_rect(fill = "white", colour = "grey50")
  ) #+

data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v2.csv')
pairs <- data.table::fread(pair_groups_fn) %>% as.data.frame()
dim(pairs)
unique(pairs$datatag)

mapping <- data.frame(datatag1=c('SA609','SA535','SA1035'), series=c("Pt4","Pt5","Pt6"))

maps <- c('SA609','SA535','SA1035')
names(maps) <- c("Pt4","Pt5","Pt6")
# pairs <- pairs %>% left_join(mapping, by='series')
pairs$datatag <- ifelse(pairs$datatag=="",maps[pairs$series],pairs$datatag)
View(pairs)
unique(pairs$comp_type)
pairs$comp_type <- ifelse(pairs$comp_type=='','treated_vs_untreated',pairs$comp_type)

pairs <- pair_groups
pairs$passage <- ifelse(pairs$labels_detailed!='',substr(pairs$labels_detailed, 1, 2),'')
data.table::fwrite(pairs, paste0(data_dir,'comparisons_drug_res_v3.csv'))




library(dplyr)
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
df1 <- data.table::fread(paste0(input_dir, 'comparisons_drug_res_v4.csv')) %>% as.data.frame()
colnames(df1)
df2 <- data.table::fread(paste0(input_dir, 'comparisons_final.csv'), header=T) %>% as.data.frame()

dim(df2)
View(head(df2))
View(head(df1))
unique(df1$comp_type)
df1$series[1]
# [1] "result_fn"       "datatag"        
# [3] "clone1"          "clone2"         
# [5] "comp_type"       "file_header"    
# [7] "labels_detailed" "labels_short"   
# [9] "group"           "series"         
# [11] "passage"         "order" 
df2 <- df2 %>% 
  dplyr::rename(result_fn=file_header)

df2$file_header <- sapply(strsplit(df2$result_fn, '_'), function(x){
  return(paste0(x[1],'_',x[2],'_',x[3]))
})

df2$datatag <- sapply(strsplit(df2$result_fn, '_'), function(x){
  return(as.character(x[2]))
})

df2$clone1 <- sapply(strsplit(df2$result_fn, '_'), function(x){
  return(as.character(x[6]))
})
df2$clone2 <- sapply(strsplit(df2$result_fn, '_'), function(x){
  return(as.character(x[8]))
})
df2$comp_type <- "treated_vs_untreated"
df2$labels_detailed <- gsub('/','vs.',df2$labels_detailed)

df2$labels_short <- sapply(strsplit(df2$labels_detailed, ' '), function(x){
  t <- as.character(x[2])
  t2 <- strsplit(t, ':')
  return(as.character(t2[[1]][1]))
})
df2$V1 <- NULL
df2$order <- 0
View(head(df2))
dim(df2)
dim(df1)
df11 <- df1 %>%
  dplyr::filter(comp_type=='untreated') %>%
  dplyr::select(all_of(colnames(df2)))
df3 <- rbind(df11, df2)
dim(df3)
df3$order
data.table::fwrite(df3, paste0(input_dir, 'comparisons_final_Hoa.csv'))

df3 <- data.table::fread(paste0(input_dir, 'comparisons_final_Hoa.csv')) %>% as.data.frame()
colnames(df3)
dim(df3)
summary(as.factor(df3$comp_type))
df3 <- df3 %>%
  dplyr::mutate(
    comp_type = case_when(
      grepl('RxH',labels_detailed) ~ 'drugHoliday_vs_untreated',
     TRUE ~ comp_type          
    ))

df4 <- df3 %>%
  dplyr::filter(comp_type %in% c('drugHoliday_vs_untreated','untreated'))
    
df5 <- df3 %>%
  dplyr::filter(comp_type =='treated_vs_untreated')
dim(df5)
dim(df4)
dim(df3)
# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/comparisons/'
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
# input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
data.table::fwrite(df3, paste0(input_dir, 'comparisons_drug_res_v5.csv'))
data.table::fwrite(df4, paste0(input_dir, 'comparisons_drug_res_v5_untreated_DH.csv'))
data.table::fwrite(df5, paste0(input_dir, 'comparisons_drug_res_v5_treated.csv'))
df3 <- data.table::fread(paste0(input_dir, 'comparisons_drug_res_v5.csv')) %>% as.data.frame()
df3$file_header[1]
tmp
df5$order <- NULL
df5$order_desc <- paste0(df5$series,"_",df5$file_header)
order_df <- data.frame(order_desc=unique(df5$order_desc), order=seq(1:length(unique(df5$order_desc))))
df5 <- df5 %>% left_join(order_df,by='order_desc')%>% 
        dplyr::select(-order_desc)
df5 <- df5 %>%
  dplyr::filter(file_header!='scrande_SA1035_3')

df4$order

View(df5)



library(dplyr)
input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
df <- data.table::fread(paste0(input_dir, 'comparisons_across_treatments_revision35.csv')) %>% as.data.frame()
colnames(df)
df$result_fn[1]
df$file_header <- sapply(strsplit(df$result_fn, '_'), function(x){
  return(paste0(x[1],'_',x[2],'_',x[3]))
})

df$datatag <- sapply(strsplit(df$result_fn, '_'), function(x){
  return(as.character(x[2]))
})
# unique(df$datatag)
df$clone1 <- sapply(strsplit(df$result_fn, '_'), function(x){
  return(as.character(x[6]))
})
unique(df$clone1)
df$clone2 <- sapply(strsplit(df$result_fn, '_'), function(x){
  return(as.character(x[8]))
})
unique(df$clone2)
df$clone2 <- ifelse(df$clone2=='B-C','B',df$clone2)
df$labels_detailed <- gsub('/','vs.',df$labels_detailed)
# View(df)
data.table::fwrite(df, paste0(input_dir, 'comparisons_across_treatments_revision35.csv'))

