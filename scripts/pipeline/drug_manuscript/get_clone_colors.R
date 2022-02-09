# Get clone colors 

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/dlp_results/'

cells_clone <- data.table::fread(paste0(input_dir,'SA609_cisplatin_cell_clones.csv')) %>% as.data.frame()
dim(cells_clone)
levels <- unique(cells_clone$clone_id)
none_clones <- c('Un','Unassigned','unassigned','None')
levels <- levels[!levels %in% none_clones]


input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA501_v2/tree_viz_dream/'
p <- readRDS(paste0(input_dir,'summary_tree.rds'))
col <- p$data$color
col[col=='white'] <- 'black'

# read.newick('/home/htran/storage/datasets/drug_resistance/dlp_results/SA501_v2/tree.newick')
library(ggplot2)
library(dplyr)
p1 <- p + scale_color_manual(values=col)



make_clone_palette <- function(levels) {
  # install.packages("inlmisc", dependencies = TRUE)  # TO DO: check this package
  # clone_names <- sort(levels)
  clone_names <- levels
  pal <- as.character(inlmisc::GetColors(length(clone_names)))  
  # if (length(levels) <= 12 & length(levels)>8) {
  #   pal <- brewer.pal(max(length(levels), 3), "Set3")
  # } else if (length(levels) <= 20 & length(levels) > 12) {
  #   pal <- clone_palette_20
  # } else {
  #   pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels))
  #   print("WARNING: more clones than palette can accomodate!")
  # }
  pal[pal=='#E7EBFA']<- '#ABB9ED'
  clone_names <- c('D','G','H','E','NULL','C','B','A')
  names(pal) <- clone_names
  pal <- pal[levels]
  pal <- pal[!is.na(pal)]
  t <- '#191919'
  names(t) <- 'R'
  pal <- c(pal,t)
  return(pal)
}

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
color_df <- data.table::fread(paste0(base_dir,'figs/colorcode_total.csv'))
dim(color_df)
head(color_df)
color_df <- color_df %>%
  dplyr::filter(datatag!='SA609')
df <- data.frame(clone_id=names(pal), colour=pal)
df$datatag <- 'SA609'
df <- df %>%
  dplyr::select(datatag, clone_id, colour)
dim(df)
color_df <- dplyr::bind_rows(color_df, df)
summary(as.factor(color_df$datatag))

data.table::fwrite(color_df, paste0(base_dir,'colorcode_total_v2.csv'))
colnames(color_df)


predefined_colors <- data.table::fread(paste0(base_dir,'all_clones_new_colors.csv')) %>% as.data.frame()

unique(predefined_colors$datatag)
# SA1035
df <- predefined_colors %>%
  dplyr::filter(grepl('SA1035',datatag))
df <- df[!duplicated(df$clone_id),]

df <- df %>%
  dplyr::select(clone_id, colour,family)%>%
  dplyr::rename(datatag=family)%>%
  dplyr::select(datatag, clone_id, colour)
df$datatag
df$colour <- stringr::str_sub(df$colour,1,7)
dim(df)
color_df <- color_df %>%
  dplyr::filter(datatag!='SA1035')
color_df <- dplyr::bind_rows(color_df, df)


unique(predefined_colors$datatag)
tag <- 'SA535'
df <- predefined_colors %>%
  dplyr::filter(grepl(tag,datatag))
df <- df[!duplicated(df$clone_id),]

df <- df %>%
  dplyr::select(clone_id, colour,family)%>%
  dplyr::rename(datatag=family)%>%
  dplyr::select(datatag, clone_id, colour)
df$datatag <- tag
df$colour <- stringr::str_sub(df$colour,1,7)
dim(df)
df$clone_id
color_df <- color_df %>%
  dplyr::filter(datatag!=tag)
color_df <- dplyr::bind_rows(color_df, df)

data.table::fwrite(color_df, paste0(base_dir,'colorcode_total_v2.csv'))


base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
color_df <- data.table::fread(paste0(base_dir,'figs/colorcode_total_v2.csv'))
summary(as.factor(color_df$datatag))

color_df <- color_df %>%
  dplyr::filter(datatag!='SA535')
dim(color_df)
cols_use <- c('J'='#E67F33','I'='#ABB9ED','H'='#CE2220','G'='#57A2AC',
              'F'='#7EB875','E'='#D0B440','D'='#4E79C4','A'='#521913','B'='#B997C6',
              'C'='#824D99')

length(cols_use)
color_df1 <- data.frame(clone_id=names(cols_use),colour=cols_use)
color_df1$datatag <- 'SA535'
color_df <- dplyr::bind_rows(color_df, color_df1)
summary(as.factor(color_df$datatag))
data.table::fwrite(color_df, paste0(base_dir,'figs/colorcode_total_v2.csv'))
