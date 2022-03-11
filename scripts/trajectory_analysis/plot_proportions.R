library(ggplot2)
# data <- out

# Donut plot
get_proportion_plt <- function(data, plot_label=F, cols=NULL){
  data <- data %>%
    dplyr::filter(count>0)
  if(is.null(cols)){
    # cyan, corn flower, cobalt, dark blue, midnight blue
    rx_cols <- c('#00FFFF','#6495ED','#0047AB', '#00008B','#191970')
    names(rx_cols) <- c('1Rx','2Rx','3Rx','4Rx','5Rx')
    #light grey, darkgray, dimgray, dark grey, dark dark grey
    # unrx_cols <- c('#D3D3D3','#A9A9A9','#696969', '#505050','#383838')
    unrx_cols <- c('#989898','#676767','#363636', '#191919','#0f0f0f')
    names(unrx_cols) <- c('1UnRx','2UnRx','3UnRx','4UnRx','5UnRx')
    
    rxh_cols <- c('#D5FF00','#FFFF00','#FFD500', '#FFBF00','#767600')
    names(rxh_cols) <- c('1RxH','2RxH','3RxH','4RxH','5RxH')
    
    cols <- c(rx_cols, unrx_cols, rxh_cols)
  }
  
  # Compute percentages
  # data$label <- factor(data$label, levels = names(cols))
  rownames(data) <- data$label
  data <- data[names(cols)[names(cols) %in% as.character(data$label)],]
  data$fraction <- data$count / sum(data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax <- cumsum(data$fraction)
  
  # Compute the bottom of each rectangle
  data$ymin <- c(0, head(data$ymax, n=-1))
  
  # Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2
  
  # Compute a good label
  # data$label <- data$category
  # data$label <- c('1Rx','2Rx','3Rx','1RxH')
  cols_use <- cols[as.character(unique(data$label))]
  if(plot_label==F){
    p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
      geom_rect() +
      # geom_label( x=3.5, aes(y=labelPosition, label=label), size=3.1, label.size=0) +
      scale_fill_manual(values = cols_use)   
  }else{
    p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
      geom_rect() +
      geom_label(x=3.5, aes(y=labelPosition, label=label), size=3.1, label.size=0) +
      scale_fill_manual(values = cols_use)
  }
  p <- p + coord_polar(theta="y") +
            xlim(c(2, 4)) +
            theme_void() +
            theme(legend.position = "none")
            # panel.background = element_rect(fill='transparent'),
            # plot.background = element_rect(fill='transparent', color=NA),
            # panel.grid.major = element_blank(),
            # panel.grid.minor = element_blank()
          # p
  return(p)
  
}

metacells <- colData(sce) %>% as.data.frame()
library(stringr)
metacells$label <- ifelse(grepl('T$',metacells$treatmentSt),'Rx','UnRx')
metacells$label <- ifelse(grepl('TU$',metacells$treatmentSt),'RxH',metacells$label)
if(datatag=='SA535'){
  metacells$label <- ifelse(metacells$label %in% c('Rx','RxH'),paste0(str_count(metacells$treatmentSt,'T'),metacells$label),
                            paste0(str_count(metacells$treatmentSt,'U')-1,metacells$label))
  
}else{
  metacells$label <- ifelse(metacells$label %in% c('Rx','RxH'),paste0(str_count(metacells$treatmentSt,'T'),metacells$label),
                            paste0(str_count(metacells$treatmentSt,'U'),metacells$label))
  
}

save_dir_fg <- paste0(save_dir,'proportion/')
dir.create(save_dir_fg)
for(cl in unique(metacells$cluster_label)){
  tmp <- metacells %>% 
    dplyr::filter(cluster_label==cl)%>% 
    dplyr::select(label, cluster_label)
  print(dim(tmp))
  out <- table(tmp$label, tmp$cluster_label)  %>% as.data.frame()
  colnames(out) <- c('label','cluster', 'count')
  p <- get_proportion_plt(out)
  png(paste0(save_dir_fg,"treatmentst_cls_",cl,".png"),height = 2*50*log10(dim(tmp)[1]), 
      width=2*50*log10(dim(tmp)[1]),res = 2*72)
  print(p)
  dev.off()
}


grouping_df <- data.table::fread("/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/encoder_trajectory/grouping_SA535_v2.csv.gz") %>% as.data.frame()
dim(grouping_df)
dim(sce)
grouping_df$cell_id[1]
metacells$cell_id <- rownames(metacells)
metacells_backup <- metacells

grouping_df <- grouping_df %>%
  dplyr::select(cell_id,clone)
head(grouping_df)
library(dplyr)
metacells <- metacells %>% inner_join(grouping_df, by=c("cell_id"))
dim(metacells)
metacells$clone[1]
save_dir_fg <- paste0(save_dir,'proportion/')
dir.create(save_dir_fg)
metacells$clone[1]
metacells$clone <- gsub('_','',metacells$clone)
metacells$clone <- NULL
cls_ls <- unique(metacells$clone)
# cols <- clone_palette_20[1:length(cls_ls)]
# names(cols) <- cls_ls

res1 <- set_color_clones(metacells$clone)
lg1 <- plot_colors(res1$clone_palette, output_dir, legend_label='clone')
cols <- res1$clone_palette
for(cl in unique(metacells$cluster_label)){
  tmp <- metacells %>% 
    dplyr::filter(cluster_label==cl)%>% 
    dplyr::select(clone, cluster_label)
  print(dim(tmp))
  if(dim(tmp)[1]>0){
    out <- table(tmp$clone, tmp$cluster_label)  %>% as.data.frame()
    colnames(out) <- c('label','cluster', 'count')
    p <- get_proportion_plt(out, F, cols)
    png(paste0(save_dir_fg,"clone_cls_",cl,".png"),height = 2*50*log10(dim(tmp)[1]), width=2*50*log10(dim(tmp)[1]),res = 2*72)
    print(p)
    dev.off()
    
  }
}
