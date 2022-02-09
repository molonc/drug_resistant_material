# Resistant clone: A, D-E
# Sensitive clone: G, B

library(stringr)
pt_use <- 'hallmark/'
pattern_use <- '^HALLMARK_'
keeps <- c('de','pathway')
serie <- 'SA535_'
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/deg_analysis/'
save_dir <- input_dir
obs_pw <- c('top_up_pathway_cls_','top_down_pathway_cls_')
desc <- c('Up','Down')
obs_pw_use <- obs_pw[2]
desc_use <- desc[2]
# tags <- c('SA535_A_B','SA535_A_G','SA535_DE_B','SA535_DE_G')
tags <- c('SA535_A_UTTTT_UT','SA535_A_UTTTT_UTTT',
          'SA535_A_UTTT_UT','SA535_DE_UTTTT_UT',
          'SA535_DE_UTTTT_UTTT','SA535_DE_UTTT_UT')
pathway_ls <- list()

for(i in rep(1:length(tags),1)){
  pw_fn <- paste0(input_dir,tags[i],'/',pt_use, obs_pw_use, tags[i],'.txt')
  if(file.exists(pw_fn)){
    top_pathway <- read.table(pw_fn, sep = '\t', header=T, check.names=F)
    top_pathway$de <- tags[i] 
    
    top_pathway <- top_pathway[,colnames(top_pathway) %in% keeps]
    print(dim(top_pathway))
    # View(head(top_pathway))
    top_pathway$pathway <- str_replace_all(top_pathway$pathway, pattern_use, "")
    pathway_ls[[tags[i]]] <- top_pathway
  }
  
}

# pw_df <- do.call(rbind, pathway_ls)
# pw_df$Count <- 1
# dim(pw_df)

pw_df_SA535 <- do.call(rbind, pathway_ls)
pw_df_SA535$serie <- 1
pw_df_SA535$label <- 'SA535'
dim(pw_df_SA535)


# pw_df$Count <- as.factor(pw_df$Count)
# # library(ggplot2)
# col <- c("#006400","#4B0082")
# 
# names(col) <- desc
# col_use <- col[desc_use]
# 
# pw_df$de <- str_replace_all(pw_df$de, serie, "")
# p1 <- ggplot(pw_df, aes(x = de, y = pathway)) 
# p1 <- p1 + geom_point(colour = col_use, size = 5) #aes(colour = Count) "#4B0082" dark purple
# p1 <- p1 + theme(axis.text.x = element_text(size=9, angle=90), 
#                  axis.text.y = element_text(size=13),
#                  legend.position = "top",
#                  panel.grid.major = element_blank(), 
#                  panel.grid.minor = element_blank())
# # p1 <- p1  + theme_classic()
# p1 <- p1 + labs(x="", title=paste0(serie,'DE analysis resistance late vs early'),
#                 subtitle = paste0(desc_use,'-regulated pathways'))
# p1

# 2*9.5*nrow(pw_df)
# png(paste0(save_dir,paste0("PW_",desc_use,"_resistant_late_early_",serie,".png")), 
#     height = , width=2*700,res = 2*72)
# print(hm)
# dev.off()






pt_use <- 'hallmark/'
pattern_use <- '^HALLMARK_'
keeps <- c('de','pathway')
serie <- 'SA609_'
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609_rna/deg_analysis/'
save_dir <- input_dir
obs_pw <- c('top_up_pathway_cls_','top_down_pathway_cls_')
desc <- c('Up','Down')
obs_pw_use <- obs_pw[2]
desc_use <- desc[2]

# tags <- c('SA609_UTT_UUUUUUUU_R_H','SA609_UTTTT_UUUUUUUU_R_E',
#           'SA609_UTTTT_UT','SA609_R_UTTTT_UTT')

# tags <- c('SA609_R_UTTT_H_UUUUUUUU','SA609_UTT_UUUUUUUU_R_H',
# 'SA609_UTTTT_UUUUUUUU_R_H','SA609_UTTTT_UUUUUUUU_R_E',
# 'SA609_UTTTT_UT','SA609_R_UTTTT_UTT')
tags <- c('SA609_R_UTTTT_UTT','SA609_R_UTTTT_UTTT',
          'SA609_R_UTTT_UTT')
pathway_ls <- list()

for(i in rep(1:length(tags),1)){
  pw_fn <- paste0(input_dir,tags[i],'/',pt_use, obs_pw_use, tags[i],'.txt')
  if(file.exists(pw_fn)){
    top_pathway <- read.table(pw_fn, sep = '\t', header=T, check.names=F)
    top_pathway$de <- tags[i] 
    
    top_pathway <- top_pathway[,colnames(top_pathway) %in% keeps]
    print(dim(top_pathway))
    # View(head(top_pathway))
    top_pathway$pathway <- str_replace_all(top_pathway$pathway, pattern_use, "")
    pathway_ls[[tags[i]]] <- top_pathway
  } else{
    print(paste0('Dont exist pathway: ',desc_use,' :', tags[i]))
    print(pw_fn)
  }
  
}
pw_df_SA609 <- do.call(rbind, pathway_ls)
pw_df_SA609$serie <- 2
pw_df_SA609$label <- 'SA609'
dim(pw_df_SA609)

# SA1035
serie <- 'SA1035'
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_rna/deg_analysis/'
save_dir <- input_dir
obs_pw <- c('top_up_pathway_cls_','top_down_pathway_cls_')
desc <- c('Up','Down')
obs_pw_use <- obs_pw[2]
desc_use <- desc[2]
tags <- c('SA1035_H_E')
pathway_ls <- list()
for(i in rep(1:length(tags),1)){
  pw_fn <- paste0(input_dir,tags[i],'/',pt_use, obs_pw_use, tags[i],'.txt')
  if(file.exists(pw_fn)){
    top_pathway <- read.table(pw_fn, sep = '\t', header=T, check.names=F)
    top_pathway$de <- tags[i] 
    
    top_pathway <- top_pathway[,colnames(top_pathway) %in% keeps]
    print(dim(top_pathway))
    # View(head(top_pathway))
    top_pathway$pathway <- str_replace_all(top_pathway$pathway, pattern_use, "")
    pathway_ls[[tags[i]]] <- top_pathway
  } else{
    print(paste0('Dont exist pathway: ',desc_use,' :', tags[i]))
    print(pw_fn)
  }
  
}


pw_df_SA1035 <- do.call(rbind, pathway_ls)
pw_df_SA1035$serie <- 3
pw_df_SA1035$label <- 'SA1035'
dim(pw_df_SA1035)

pw_df <- do.call(rbind, list(pw_df_SA535, pw_df_SA609, pw_df_SA1035))
pw_df <- do.call(rbind, list(pw_df_SA535, pw_df_SA609))

dim(pw_df)
library(reshape2)
get_label <- function(de) {
  labels <- sapply(strsplit(de, "_"), function(x) {
    return(x[1])
  })
  return(labels)
}
stat_pw <- acast(pw_df, pathway~de, value.var="serie")
# class(stat_pw)
# View(stat_pw)
# dim(stat_pw)
colnames(stat_pw)
clusters <- rowSums(!is.na(stat_pw), na.rm=T) >= 3
clusters
clusters_label <- ifelse(clusters==TRUE,'common_pw',
                         ifelse(clusters==FALSE,'specific_pw',NA))
# length(clusters)
# common_pw <- rownames(stat_pw)[clusters]
# not_common_pw <- rownames(stat_pw)[!clusters]

col = c("#006400","#8B0000","#00008B")
col = c("#006400","#8B0000")
# heatmap_legend_param = list(at = c(1, 2, 3),
#                             labels = c('SA535','SA609','SA1035'),
#                             title = "Serie")
# lgd = Legend(labels = c('SA535','SA609','SA1035'), title = "Serie")
hm <- ComplexHeatmap::Heatmap(stat_pw, na_col = "white",
                              show_column_names=T,
                              show_row_names = T,
                              col = col,
                              cluster_rows=F,cluster_columns=F,
                              column_split = get_label(colnames(stat_pw)),
                              row_split = clusters_label,
                              name = paste0(desc_use,'-regulated pathways in 3 series'),
                              column_names_gp = grid::gpar(fontsize = 7),
                              row_names_gp = grid::gpar(fontsize = 10),
                              show_heatmap_legend = F)



# tag <- "_resistant_vs_sensitive_"
tag <- "_late_vs_early_"
serie <- 'SA535_SA609_SA1035'
png(paste0(save_dir,paste0("PW_",desc_use,tag,serie,".png")), 
    height = 850, width=2*800,res = 2*72)
print(hm)
dev.off()

# ComplexHeatmap::Heatmap(stat_pw, col = col)





length(pathway_ls)
pw_df <- do.call(rbind, pathway_ls)
pw_df$Count <- 1
dim(pw_df)

pw_df$Count <- as.factor(pw_df$Count)
# library(ggplot2)
col <- c("#006400","#4B0082")

names(col) <- desc
col_use <- col[desc_use]
serie <- 'SA535_SA609_SA1035'
serie <- 'SA609'
# pw_df$de <- str_replace_all(pw_df$de, serie, "")
p1 <- ggplot(pw_df, aes(x = de, y = pathway)) 
p1 <- p1 + geom_point(colour = col_use, size = 5) #aes(colour = Count) "#4B0082" dark purple
p1 <- p1 + theme(axis.text.x = element_text(size=10, angle=90), 
                 axis.text.y = element_text(size=11),
                 legend.position = "top",
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank())
# p1 <- p1  + theme_classic()
p1 <- p1 + labs(x="", title=paste0(serie,' DE analysis resistance vs sensitive to drug'),
                subtitle = paste0(desc_use,'-regulated pathways'))
# p1
png(paste0(save_dir,paste0("PW_",desc_use,"_resistant_sensitive_",serie,".png")), 
    height = 2*17*nrow(pw_df)+130, width=2*800,res = 2*72)
print(p1)
dev.off()







tag <- 'SA535_A_B'
up_SA535_A_B <- read.csv(paste0(input_dir,tag,'/',pt_use,'SA535_up_pathway.csv'),
                          row.names=1, check.names=F, stringsAsFactors=F)
pt_ls <- rownames(up_SA535_A_B)
pt_ls <- str_replace_all(tolower(pt_ls), "^hallmark_", "")
rownames(up_SA535_A_B)  <- pt_ls
# rownames(up_SA535_GF_D)  <- paste0('SA535_',pt_ls)
dim(up_SA535_A_B)
max(up_SA535_A_B, na.rm = T)
View(up_SA535_A_B)
top_pathway_SA535_A_B <- read.table(paste0(input_dir,tag,'/',pt_use,'top_up_pathway_cls_',tag,'.txt'),
                              sep = '\t', header=T, check.names=F)
top_pathway_SA535_A_B$de <- tag 

top_pathway_SA535_A_B <- top_pathway_SA535_A_B[,colnames(top_pathway_SA535_A_B) %in% keeps]
View(head(top_pathway_SA535_A_B))
top_pathway_SA535_A_B$pathway <- str_replace_all(top_pathway_SA535_A_B$pathway, pattern_use, "")


d_pathway_SA535_A_B <- read.table(paste0(input_dir,tag,'/',pt_use,'top_down_pathway_cls_',tag,'.txt'),
                                    sep = '\t', header=T, check.names=F)
d_pathway_SA535_A_B$de <- tag 

d_pathway_SA535_A_B <- d_pathway_SA535_A_B[,colnames(d_pathway_SA535_A_B) %in% keeps]
dim(d_pathway_SA535_A_B)

View(top_pathway_SA535_A_B)
dim(top_pathway_SA535_A_B)
tag <- 'SA535_A_G'
top_pathway_SA535_A_G <- read.table(paste0(input_dir,tag,'/',pt_use,'top_up_pathway_cls_',tag,'.txt'),
                                    sep = '\t', header=T, check.names=F)
top_pathway_SA535_A_G$de <- tag 

top_pathway_SA535_A_G <- top_pathway_SA535_A_G[,colnames(top_pathway_SA535_A_G) %in% keeps]
dim(top_pathway_SA535_A_G)

d_pathway_SA535_A_G <- read.table(paste0(input_dir,tag,'/',pt_use,'top_down_pathway_cls_',tag,'.txt'),
                                  sep = '\t', header=T, check.names=F)
d_pathway_SA535_A_G$de <- tag 
d_pathway_SA535_A_G <- d_pathway_SA535_A_G[,colnames(d_pathway_SA535_A_G) %in% keeps]
dim(d_pathway_SA535_A_G)


tag <- 'SA535_DE_B'
top_pathway_SA535_DE_B <- read.table(paste0(input_dir,tag,'/',pt_use,'top_up_pathway_cls_',tag,'.txt'),
                                    sep = '\t', header=T, check.names=F)
top_pathway_SA535_DE_B$de <- tag 

top_pathway_SA535_DE_B <- top_pathway_SA535_DE_B[,colnames(top_pathway_SA535_DE_B) %in% keeps]
dim(top_pathway_SA535_DE_B)

d_pathway_SA535_DE_B <- read.table(paste0(input_dir,tag,'/',pt_use,'top_down_pathway_cls_',tag,'.txt'),
                                  sep = '\t', header=T, check.names=F)
d_pathway_SA535_DE_B$de <- tag 
d_pathway_SA535_DE_B <- d_pathway_SA535_DE_B[,colnames(d_pathway_SA535_DE_B) %in% keeps]
dim(d_pathway_SA535_DE_B)

tag <- 'SA535_DE_G'
top_pathway_SA535_DE_G <- read.table(paste0(input_dir,tag,'/',pt_use,'top_up_pathway_cls_',tag,'.txt'),
                                     sep = '\t', header=T, check.names=F)
top_pathway_SA535_DE_G$de <- tag 

top_pathway_SA535_DE_G <- top_pathway_SA535_DE_G[,colnames(top_pathway_SA535_DE_G) %in% keeps]
dim(top_pathway_SA535_DE_G)

d_pathway_SA535_DE_G <- read.table(paste0(input_dir,tag,'/',pt_use,'top_down_pathway_cls_',tag,'.txt'),
                                   sep = '\t', header=T, check.names=F)
d_pathway_SA535_DE_G$de <- tag 
d_pathway_SA535_DE_G <- d_pathway_SA535_DE_G[,colnames(d_pathway_SA535_DE_G) %in% keeps]
dim(d_pathway_SA535_DE_G)

df_ls <- list(top_pathway_SA535_A_B, top_pathway_SA535_A_G,
              top_pathway_SA535_DE_B, top_pathway_SA535_DE_G)

df_ls <- list(d_pathway_SA535_A_B, d_pathway_SA535_A_G,
              d_pathway_SA535_DE_B, d_pathway_SA535_DE_G)
pw_df <- do.call(rbind, df_ls)
dim(pw_df)
ps <- table(pw_df)
View()
ps <- as.data.frame(ps)
table(ps$pathway, ps$Freq)
colnames(ps)
p1 <- ggplot(ps, aes(x = de, y = pathway)) +
  theme_dose(12)
p1 <- p1 + geom_point(aes(colour = Freq), size = 2)
p1

png(paste0(save_dir,"DE_down_sensitive_SA535.png"), height = 2*700, width=2*1000,res = 2*72)
print(p1)
dev.off()