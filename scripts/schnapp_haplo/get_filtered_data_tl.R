
save_dir <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/"
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'

st <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/bin_cnvs_corrupt_double_padding.csv'
states <- read.csv(st, check.names = F, stringsAsFactors = FALSE)
dim(states)
View(head(states[1:5,1:5]))


cnstates <- read.csv(paste0(save_dir,'total_merged_filtered_states_chr.csv'), check.names = F, stringsAsFactors = FALSE)
View(head(cnstates[1:5,1:5]))
cnstates

cnstates$chr_lb <- paste0(cnstates$chr,'_',cnstates$start,'_',cnstates$end)
states$chr_lb <- paste0(states$chr,'_',states$start,'_',states$end)

padding_chr <- setdiff(states$chr_lb, cnstates$chr_lb)

padding_chr_label <- paste0('locus_',padding_chr)

sum(syn_features$loci %in% padding_chr_label)  # 0

sum(filtered$loci %in% padding_chr_label)  # 0


sum(features$loci %in% padding_chr)  
ff <- unique(features$loci)
length(ff)

ft <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/corrupt_tree_straightened_features.csv'
features <- read.csv(ft, check.names = F, stringsAsFactors = FALSE)
View(head(features[1:5,]))

features <- read.csv(paste0(save_dir,'corrupt_tree_features.csv'), check.names = F, stringsAsFactors = FALSE)
dim(features)
View(head(features))

syn_features <- read.csv(paste0(save_dir,'corrupt_tree_sync-cnv_features.csv'), check.names = F, stringsAsFactors = FALSE)
dim(syn_features)
View(syn_features[1:5,])
syn_features_chr7
chr7_bins <- grep('locus_7_', syn_features$loci, value=F)
length(chr7_bins)


filtered <- read.csv(paste0(save_dir,'filtered.csv'), check.names = F, stringsAsFactors = FALSE)
dim(filtered)
View(filtered[1:5,1:2])
chr7_bins_filtered <- grep('locus_7_', filtered$loci, value=F)
length(chr7_bins_filtered)


# saveRDS(total_bincnv, file = paste0(save_dir,prefix_filtered,'_CNbins_data.rds'))
prefix_filtered='total_filtered'
total_bincnv <- readRDS(paste0(save_dir,'schnapps_out/',prefix_filtered,'_CNbins_data.rds'))
total_filtered_cn <- read.csv(paste0(results_dir,'total_merged_filtered_states.csv'),
                              header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)

cn_tyler <- read.csv(paste0(save_dir,'schnapps_out/total_merged_filtered_states_tl1.csv'),
                              header=T, check.names = F,stringsAsFactors = FALSE)

cn_tyler$chr_lb <- paste0(cn_tyler$chr,'_',cn_tyler$start,'_',cn_tyler$end)
sum(rownames(total_filtered_cn) %in% cn_tyler$chr_lb)

rownames(cn_tyler) <- cn_tyler$chr_lb
drops <- c("chr","start","end","width","chr_lb")
cn_tyler <- cn_tyler[ , !(names(cn_tyler) %in% drops)]
out_fn <- paste0(save_dir,'total_merged_filtered_states.csv')
write.csv(cn_tyler, out_fn,  row.names = T, quote=F)

features <- read.csv(paste0(save_dir,'bin_cnvs_corrupt_double_padding.csv'),
                     check.names = F,stringsAsFactors = FALSE)

filtered_origin <- read.csv(paste0(results_dir,'filtered.csv'),
                     check.names = F,stringsAsFactors = FALSE)
dim(filtered_origin)
loci_origin <- unique(filtered_origin$loci)
length(loci_origin) #64

filtered_tyler <- read.csv(paste0(save_dir,'filtered.csv'),
                            check.names = F,stringsAsFactors = FALSE)
dim(filtered_tyler) 
loci_tl <- unique(filtered_tyler$loci)
length(loci_tl) #86

i <- intersect(loci_origin, loci_tl) #60
length(i)
View(head(filtered_origin))


filtered$width[1:3]

dim(total_bincnv)
colnames(total_bincnv)

total_bincnv1 <- total_bincnv[,c('chr_lb','width')]
total_bincnv1$lb <- paste0(total_bincnv1$chr_lb,"_",total_bincnv1$width)

t <- total_bincnv[1:7000,]
View(t1)
t1 <- t[,c('chr_lb','width')]
# library(reshape)
# t1 <- cast(t, chr_lb ~ width, value = width)
View(t1)
t1 <- t1 %>%
  pivot_wider(names_from = width, values_from = width)
t1

unique(t1$chr_lb,t1$width)
save_dir <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local_tyler/schnapps_out/"
prefix <- 'total'
prefix_filtered <- 'total_filtered'
fn <- paste0(save_dir,prefix,'_alleleCNbins_out.rds')
total_bincnv <- readRDS(fn)
colnames(total_bincnv)

# Observe chr 7
prefix_filtered='total_filtered'
total_bincnv <- readRDS(paste0(save_dir,'schnapps_out/',prefix_filtered,'_CNbins_data.rds'))
results_dir <- '/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919X7_whole_local/'
obs_chr <- 20
visualize_chr_allele_states_phase <- function(total_bincnv, obs_chr, save_dir,results_dir){
  select <- c("chr","start","cell_id","state_phase")
  total_bincnv <- total_bincnv[ , (names(total_bincnv) %in% select)]
  total_filtered_cn <- read.csv(paste0(results_dir,'total_merged_filtered_states.csv'),
                                header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  
  dim(total_bincnv)
  colnames(total_bincnv)
  # View(head(total_bincnv))
  
  reads <- total_bincnv[total_bincnv$chr==obs_chr,]
  reads <- reads[reads$cell_id %in% colnames(total_filtered_cn),]
  dim(reads)
  # reads <- fread("*_reads.csv.gz")
  # metrics <- read.csv("*_metrics.csv.gz")
  reads$chr <- factor(reads$chr, levels = c(1:22, "X", "Y"))
  slice <- reads[, c("chr", "start", "cell_id", "state_phase")]
  # slice$pos <- paste0(slice$chr, ":", slice$start)
  # assign_val <- data.frame(state_phase=as.character(unique(slice$state_phase)),
  #                          state_val=c(1,2,3,4,5))
  # 
  # slice <- slice %>% inner_join(assign_val, by = "state_phase")
  # slice1 <- slice[, c("chr", "start", "cell_id", "state_val","pos")]
  # wide <- spread(slice1[, -c(1:2)], pos, state_val)
  # # wide <- spread(slice[, c("cell_id", "pos", "state_val")], pos, state_val)
  # rownames(wide) <- wide$cell_id
  # wide$cell_id <- NULL
  # cluster <- hclust(dist(wide), method = "ward.D")
  # slice$cell_id <- factor(slice$cell_id, level = rownames(wide)[cluster$ord])
  # cols <- c(rev(brewer.pal(n = 3, "Blues"))[1:2], "#CCCCCC", tail(brewer.pal(n = 8, "OrRd"), 6))
  # cols <- c(cols, cols[c(9, 9, 9)])
  # names(cols) <- 1:12
  
  
  cn_colours_phase <- schnapps:::scCNphase_colors
  cols <- cn_colours_phase
  
  g <- ggplot(slice, aes(cell_id, start, fill = as.factor(state_phase))) + 
    geom_tile() + 
    scale_fill_manual(values = cols) + 
    # facet_grid(~chr, scales = "free") + 
    facet_grid(~chr, scales = "free", space = "free", switch = "x") +
    scale_y_continuous(expand = c(0, 0), breaks = NULL)
  
  g <- g + labs(x='Cells',y=paste0('Chromosome ',obs_chr),title=paste0('Chr',obs_chr,' ASCN State'))
  
  # panel.spacing = unit(0.1, "lines"),
  g <- g + theme(plot.title = element_text(color="black", size=15,hjust = 0.5),
                 axis.title = element_text(color="black", size=11,hjust = 0.5),
                 axis.text = element_blank())
  # g
  
  png(file = paste0(save_dir,"ascn_heatmap_chr",obs_chr,".png"), height = 2*300, width=2*900,res = 2*72)
  print(g)
  dev.off()
  
  saveRDS(g,file=paste0(save_dir,"ascn_chr",obs_chr,"_plots.rds"))
  
  
}  



library(ggplot2)
visualize_chr_cp_states <- function(total_bincnv, obs_chr, save_dir,results_dir){
  select <- c("chr","start","cell_id","state")
  total_bincnv <- total_bincnv[ , (names(total_bincnv) %in% select)]
  total_filtered_cn <- read.csv(paste0(results_dir,'total_merged_filtered_states.csv'),
                                header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  
  dim(total_bincnv)
  dim(total_bincnv_ch7)
  colnames(total_bincnv)
  # View(head(total_bincnv))
  
  reads <- total_bincnv[total_bincnv$chr==obs_chr,]
  reads <- reads[reads$cell_id %in% colnames(total_filtered_cn),]
  dim(reads)
  # reads <- fread("*_reads.csv.gz")
  # metrics <- read.csv("*_metrics.csv.gz")
  reads$chr <- factor(reads$chr, levels = c(1:22, "X", "Y"))
  slice <- reads[, c("chr", "start", "cell_id", "state")]
  slice$pos <- paste0(slice$chr, ":", slice$start)
  wide <- spread(slice[, -c(1:2)], pos, state)
  rownames(wide) <- wide$cell_id
  wide$cell_id <- NULL
  cluster <- hclust(dist(wide), method = "ward.D")
  slice$cell_id <- factor(slice$cell_id, level = rownames(wide)[cluster$ord])
  # cols <- c(rev(brewer.pal(n = 3, "Blues"))[1:2], "#CCCCCC", tail(brewer.pal(n = 8, "OrRd"), 6))
  # cols <- c(cols, cols[c(9, 9, 9)])
  # names(cols) <- 1:12
  
  cn_colours <- structure(
    c(
      "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
      "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
    ),
    names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
  )
  cols <- cn_colours
  
  g <- ggplot(slice, aes(cell_id, start, fill = as.factor(state))) + 
    geom_tile() + 
    scale_fill_manual(values = cols) + 
    # facet_grid(~chr, scales = "free") + 
    facet_grid(~chr, scales = "free", space = "free", switch = "x") +
    scale_y_continuous(expand = c(0, 0), breaks = NULL)
  
  g <- g + labs(x='Cells',y=paste0('Chromosome ',obs_chr),title=paste0('Chr',obs_chr,' Copy Number State'))
  
  # panel.spacing = unit(0.1, "lines"),
  g <- g + theme(plot.title = element_text(color="black", size=15,hjust = 0.5),
                 axis.title = element_text(color="black", size=11,hjust = 0.5),
                 axis.text = element_blank())
  # g
  
  png(file = paste0(save_dir,"schnapps_out/hmmcopy_heatmap_chr",obs_chr,".png"), height = 2*300, width=2*900,res = 2*72)
  print(g)
  dev.off()
  
  saveRDS(g,file=paste0(save_dir,"schnapps_out/chr",obs_chr,"_plots.rds"))
}  

# lim <- max(quantile(tmp$copy, na.rm = TRUE, 0.99), 4)
# ggplot(tmp, aes(start, copy * 2, col = as.factor(state))) + 
#   geom_point(size = 0.5) + 
#   facet_grid(cell_id ~ chr, sWarning message:, space = "free_x", switch = "x") + 
#   scale_x_continuous(expand = c(0, 0), breaks = NULL) + 
#   theme(paRemoved 2097 rows containing missing values (geom_point). aks = seq(0, 18, by = 2), limit = c(0, lim)) + scale_col> _manual(values = cols)
# 
# 



# Number of reads per cells 
total_filtered_cn
total_reads_percell <- colSums(total_filtered_cn)
total_reads_df <- data.frame(read_val=total_reads_percell,cell_id=colnames(total_filtered_cn))

#, color = as.factor(read_val)
# p <- ggplot(total_reads_df, aes(x=cell_id, y=read_val)) +
#   geom_point()
# p
total_reads_df$read_val <- log(total_reads_df$read_val)
p <- ggplot(total_reads_df, aes(cell_id, read_val)) +
     geom_point()

p <- p + labs(x='Cells',y='Total reads per cell (log)',title='')

# panel.spacing = unit(0.1, "lines"),
p <- p + theme(plot.title = element_text(color="black", size=15,hjust = 0.5),
               axis.title = element_text(color="black", size=13,hjust = 0.5),
               axis.text.x = element_blank(),
               axis.title.y = element_text(color="black",size=11,hjust = 0.5))

png(file = paste0(save_dir,"total_reads_dlp.png"), height = 2*500, width=2*1000,res = 2*72)
print(p)
dev.off()