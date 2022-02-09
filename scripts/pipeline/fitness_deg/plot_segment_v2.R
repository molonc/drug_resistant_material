suppressPackageStartupMessages({
  require("optparse")
  require("data.table")
  require("dplyr")
  require(RColorBrewer)
  require(glue)
  require(ggrepel)
  require(tidyr)
  require(scales)
  require(ComplexHeatmap)
})

# source_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_Tyler_v3/library_groupings.csv'
options(scipen = 999)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

get_chrs <- function(loci) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(labels)
}
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/'
results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total/'
save_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_total/deg_analysis/SA535_GF_treated_D/'
segment_fn <- paste0(input_dir,'clonealign/raw_segments.rds')
cellclones_fn <- paste0(results_dir, '/cell_clones.csv')
obs_clusters <- c('G', 'D')

filtered_genes <- read.table(paste0(input_dir,'clonealign/SA535X4XB02498/SA535X4XB02498.txt'), 
                             sep='\t',check.names=F, header=T)
View(head(filtered_genes))
summary(as.factor(filtered_genes$cluster))


total_cn_states <- read.csv(paste0(results_dir, '/total_merged_filtered_states.csv'), header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)

plotCNV <- function(segment_fn, obs_clusters, cellclones_fn, 
                    input_dir, output_file){
  cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD','10'='#C196C4','11'='#D0BAD8')
  
  chr_levels <- c(as.character(1:23), "X", "Y")
  
  save_dir <- paste0(dirname(output_file),'/')
  input_dir <- paste0(input_dir,'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  cell_clones <- read.delim(cellclones_fn, stringsAsFactors = FALSE, sep = ",")
  print(dim(cell_clones))
  colnames(cell_clones)
  # rownames(cell_clones) <- cell_clones$cell_id
  segments <- readRDS(segment_fn)
  dim(segments)
  segments$sample[1:3]
  obs_cells <- cell_clones[cell_clones$clone_id %in% obs_clusters,'cell_id']
  print(length(obs_cells))
  segments <- segments[segments$cell_names %in% obs_cells,]
  df_cnv <- segments
  colnames(df_cnv)
  segments <- segments %>% inner_join(cell_clones, by=c('cell_names'='cell_id'))
  segments_backup <- segments
  
  t <- '1_4000001_10000000'
  View(segments[(segments$chr=='1') & (segments$start=='4000001') & (segments$end=='10000000'),])
  df_cnv_ls <- list()
  
  for(c in obs_clusters){
    df_cnv <- segments[segments$clone_id==c,]
    print(dim(df_cnv))
    df_cnv <- df_cnv %>% 
      dplyr::select(chr, start, end, copy_number) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(chr, start, end) %>% 
      summarize(copy_number=median(copy_number), .groups = 'drop')
    
    
    
    df_cnv <- drop_na(df_cnv)
    # levels(df_cnv$copy_number) <- 0:(length(cnv_colors)-1)
    # df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
    print(summary(as.factor(df_cnv$copy_number)))
    df_cnv$copy_number <- as.factor(df_cnv$copy_number)
    df_cnv$cluster <- c
    df_cnv_ls[[c]] <- df_cnv
  }
  
  
  df_cnv_total <- do.call(rbind, df_cnv_ls)
  dim(df_cnv_total)
  View(head(df_cnv_total))
  df_cnv_total <- dplyr::select(df_cnv_total, chr, start) %>% 
    group_by(chr) %>% 
    dplyr::mutate(start_order = rank(start)) %>% 
    ungroup() %>% 
    inner_join(df_cnv_total)
  
  
  dim(df_cnv_total)
  View(head(df_cnv_total))
  
  stat_cnv <- df_cnv_total
  dim(stat_cnv)
  colnames(df_cnv_total)
  df_cnv_total$cluster[1:3]
  df_cnv_total$copy_number[1:3]
  sum(is.na(df_cnv_total$chr_desc))
  df_cnv_total <- df_cnv_total[,colnames(df_cnv_total) %in% c('chr_desc','cluster','copy_number')]
  library(reshape)
  df_cnv_total <- cast(df_cnv_total,cluster ~ chr_desc)
  # loci=rownames(total_cn_states)
  # df <- parse_bin_names(loci)
  # df$start[1:3]
  # dim(df)
  # df_cnv <- df_cnv %>% inner_join(df, by=c('chr','start','end'))
  
  summary(df_cnv_total$chr)
  df_cnv_total$chr_desc <- paste0(df_cnv_total$chr,'_',df_cnv_total$start,'_',df_cnv_total$end)
  length(unique(df_cnv_total$chr_desc))
  dim(df_cnv_total)
  filtered_loci = rownames(total_cn_states)
  all_loci <- unique(df_cnv_total$chr_desc)
  ext <- filtered_loci[!(filtered_loci %in% all_loci)]
  length(ext)
  ext_all <- all_loci[!(all_loci %in% filtered_loci)]
  length(ext_all)
  length(all_loci)
  df_exclude <- data.frame(excluded_loci=ext_all)
  df_ext <- data.frame(filtered_loci=ext)
  write.csv(df_ext,file=paste0(save_dir,'ext_filtered_loci.csv'), quote=F, row.names=F)
  write.csv(df_exclude,file=paste0(save_dir,'exclued_all.csv'), quote=F, row.names=F)
  df_ext <- data.frame(filtered_loci=ext)
  t <- '1_4000001_10000000'
  View(df_cnv_total[df_cnv_total$chr_desc==t,])
  df_cnv$chr_desc[1:3]
  df_cnv <- df_cnv[df_cnv$chr_desc %in% rownames(total_cn_states),]
  dim(df_cnv)
 
  cnv_plot <- dplyr::filter(df_cnv, !chr %in% c("X","Y")) %>% 
    ggplot(aes(x = start_order, y = cluster, fill = copy_number)) +
    geom_raster() +
    facet_wrap(~ chr, scales = "free_x", nrow = 1, switch = "x") +
    #theme(legend.position = "top", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    #scale_fill_manual(values=cnv_colors, name = "Copy number", guide = 'legend',labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
    #          theme(legend.position = "bottom") +   
    scale_fill_manual(values = cnv_cols, name = "Copy number") +  #, labels = 0:(length(cnv_cols)-1),drop=FALSE) +
    labs(x = "chromosome") +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=14), 
          strip.placement = "outside",
          legend.position = "none",
          panel.spacing = unit(c(0.1), 'cm')) +   # MA: was 0.2
    theme(text = element_text(size = 18)) +
    labs(x = "Chromosome", y = "Clone")+
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_continuous(expand = c(0,0))
  

  main_plot <- cowplot::plot_grid(
    cnv_plot,
    cnv_plot,
    ncol = 1,
    rel_heights = c(2.5, 2.5),
    align = 'v'
  )
  
  png(paste0(save_dir,"cnv_profile.png"), height = 2*100*length(obs_clusters)+20, width=2*900,res = 2*72)
  print(main_plot)
  dev.off()
  
}

cnv_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/clonealign/whole_data/SA535_whole_data.csv'
df_cnv <- read.csv(cnv_fn, check.names=F, stringsAsFactors=F)
View(head(df_cnv))
dim(df_cnv)
summary(as.factor(df_cnv$cluster))
obs_clones <- c('A','B')

library(dplyr)
cnv <- dplyr::filter(df_cnv, use_gene) %>%
  dplyr::rename(clone = cluster,
                copy_number=median_cnmode) %>% 
  dplyr::select(ensembl_gene_id, clone, copy_number) %>% 
  dplyr::filter(clone %in% obs_clones) %>%
  spread(clone, copy_number)

dim(cnv_mat)

cnv_mat <- cnv %>%
  as.data.frame %>%
  column_to_rownames("ensembl_gene_id") 


class(cnv_mat)

third_quantile <- as.numeric(quantile(t, probs = 0.75))
default_threshold <- 0.25
if(third_quantile < default_threshold){
  default_threshold = third_quantile
}

var_genes <- apply(cnv_mat, 1, var)
cnv_mat <- cnv_mat[var_genes >= default_threshold,]
colnames(cnv_mat)
df_cnv <- cnv_mat



dim(cnv_mat)
df_cnv$ensembl_gene_id <- rownames(df_cnv)
colnames(df_cnv)
df_cnv <- df_cnv %>%
  pivot_longer(!ensembl_gene_id, names_to = "clone", values_to = "cnv")
dim(df_cnv)
View(head(df_cnv))
plot_CNV <- function(df_cnv, obs_clones=c('A','B'), sample_name='SA535X9XB03616',
                     additional_genes = NULL){
  
  annots <- annotables::grch37
  
  annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  
  df_cnv <- inner_join(df_cnv, annots)
  
  df_cnv <- dplyr::select(df_cnv, ensembl_gene_id, chr, start) %>% 
    group_by(chr) %>% 
    dplyr::mutate(start_order = rank(start)) %>% 
    ungroup() %>% 
    inner_join(df_cnv)
  
  chr_levels <- c(as.character(1:23), "X", "Y")
  
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  
  df_cnv <- drop_na(df_cnv)
  df_cnv$cnv <- as.character(round(df_cnv$mean_cnmode))
  
  # summary(as.factor(df_cnv$cnv))
  
  # MA: 11 Apr 2020, setting the copy number colors that were used in the heatmap
  # TO COME BACK
  #print("data frame")
  #cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  #cnv_cols <- data.frame(cn=0:11,color<-cnv_cols)
  #colnames(cnv_cols) <- c("cn","color")
  #levels(cnv_cols$color) <- cnv_colors
  
  cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364',
                  '#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  
  cnv_cols <- structure(
    c(
      "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
      "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
    ),
    names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
  )
  # cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
  #               '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD','10'='#C196C4','11'='#D0BAD8')
  #levels(cnv_cols) <- 0:11
  #  cnv_cols <- c("0" = "#2166ac",
  #                "1" = "#92c5de", 
  #                "2" = "grey80", 
  #                "3" = "#f4a582", 
  #                "4" = "#d6604d",
  #                "5" = "#b2182b",
  #                "6+" = "#67001f")
  
  # df_cnv$cnv <- as.character(round(df_cnv$median_cnmode))
  # summary(as.factor(df_cnv$cnv))
  df_cnv$cnv <- as.factor(df_cnv$cnv)
  levels(df_cnv$cnv) <- 0:(length(cnv_colors)-1)
  #print(levels(df_cnv$cnv))
  
  # MA: removing the 6+ restriction (clonealign still has that restriction though)
  # df_cnv$cnv[round(df_cnv$median_cnmode) >= 6] <- "6+"
  #!chr %in% c("X","Y")
  # View(head(df_cnv))

  cnv_plot <- dplyr::filter(df_cnv, cluster %in% obs_clones, !chr %in% c("X","Y")) %>% 
    ggplot(aes(x = start_order, y = cluster, fill = cnv)) +
    geom_raster() +
    facet_wrap(~ chr, scales = "free_x", nrow = 1, switch = "x") +
    #theme(legend.position = "top", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values=cnv_cols, name = "Copy number", guide = 'legend',
                      labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
             theme(legend.position = "bottom") +
    # scale_fill_manual(values = cnv_cols, name = "Copy Number") +  #, labels = 0:(length(cnv_cols)-1),drop=FALSE) +
    labs(x = "chromosome") +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=14), 
          strip.placement = "outside",
          legend.position = "none",
          panel.spacing = unit(c(0.1), 'cm')) +   # MA: was 0.2
    theme(text = element_text(size = 18)) +
    labs(x = "Chromosome", y = "Clone")+
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_continuous(expand = c(0,0))
  
  # cnv_plot
  main_plot <- cowplot::plot_grid(
    track_plot2 + theme(axis.title.x = element_blank(),
                        strip.text.x = element_blank()),
    cnv_plot,
    ncol = 1,
    rel_heights = c(4.5,1),
    align = 'v'
  )
  
  # main_plot
  
  png(paste0(save_dir,"DE_cnv_expression_SA535_A_vs_B.png"), height = 2*500, width=2*1200,res = 2*72)
  print(main_plot)
  dev.off()
  
  
}



input_deg_fn <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535_rna_cys/deg_analysis/SA535_A_B/de_significant_genes.txt'
plot_DE_genes_chr <- function(cnv_mat, input_deg_fn, obs_clones=c('E','H'), sample_name='',
                              additional_genes = NULL, n_genes_to_annotate = 30){
  
  # 
  deg_df <- read.table(input_deg_fn, check.names=F, stringsAsFactors=F, header=T, row.names = 1)
  print(dim(deg_df))
  # View(head(deg_df))
  deg_df <- deg_df %>% rownames_to_column("ensembl_gene_id")
  colnames(deg_df)[which(names(deg_df) == "gene_symb")] <- "gene_symbol"
  deg_df$is_genome <- deg_df$ensembl_gene_id %in% rownames(cnv_mat)
  
  deg_df$variant_gene <- ifelse(deg_df$is_genome==TRUE,'var_dna_rna',
                           ifelse(deg_df$is_genome==FALSE,'var_rna','unknown'))
  
  print(summary(as.factor(deg_df$is_genome)))
  
  clone_str <- paste(obs_clones, collapse = " vs ")
  print(clone_str)
  save_dir <- paste0(dirname(input_deg_fn),'/')
  if (!file.exists(save_dir)){
    dir.create(save_dir)
  }
  chrs <- c(as.character(1:23), "X", "Y")
  df_track <- annotables::grch37 %>% 
    dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    inner_join(deg_df) %>% 
    dplyr::filter(chr %in% c(as.character(1:23), "X", "Y")) %>% 
    dplyr::mutate(chr = factor(chr, levels = chrs),
                  position = (start + end) / 2)
  
  # MA: sometimes a genes appears more than once, keep only the unique ones
  df_track <- distinct(df_track)
  
  cols <- rev(brewer.pal(n = 3, name = "RdBu") )
  cols <- cols[1:length(unique(deg_df$is_genome))]
  
  track_plot <- ggplot(df_track, aes(x = position, y = avg_logFC)) +  #-log10(p_val_adj)
    geom_point(aes(colour = variant_gene), size = 1.5) +   # MA: size was 1
    facet_wrap(~ chr, scales = "free_x", 
               strip.position = "bottom",
               # switch = "x", 
               nrow = 1) +
    # scale_colour_gradientn(name = glue("is variant genome"),
    #                        colours = cols, 
    #                        values = c('TRUE','FALSE')) +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=10),        
          strip.placement = "outside",
          legend.position = "top",
          panel.spacing = unit(0.1, 'cm')) +   ## MA: was 0.2
    theme(text = element_text(size = 18)) +
    labs(x = "Chromosome",
         subtitle = clone_str) +
    scale_x_continuous(expand = c(0,0))
  
  
  n_genes_to_annotate <- 50
  df_annot <- top_n(df_track, n_genes_to_annotate, abs(avg_logFC))
  
  
  print("The labeled genes")
  print(df_annot)  
  
  if(!is.null(additional_genes)) {
    df_annot <- df_annot  %>% 
      bind_rows(dplyr::filter(df_track, gene_symbol %in% additional_genes))
  }
  
  ## this is where the genes are added
  track_plot2 <- track_plot +
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 3.0)
  
  base_name <- paste(obs_clones, collapse = "_vs_")
  png(paste0(save_dir,"DE_chr_",base_name,".png"), height = 2*450, width=2*1200,res = 2*72)
  print(track_plot2)
  dev.off()
}