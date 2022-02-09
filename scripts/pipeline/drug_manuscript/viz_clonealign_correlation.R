suppressPackageStartupMessages({
  library(gtools)
  library(dplyr)
  library(ggplot2)
})
library(extrafont)
font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
fonts()

get_unique_clone_id <- function(clone_labels){
  set.seed(42) 
  cls <- sapply(strsplit(clone_labels,'_'), function(x){
    if(length(x)==1){
      return(x[1])
    }else if(length(x)==2){
      idx <- sample(c(1,2),1)
      return(x[idx])
    }else{
      idx <- sample(c(1,2,3),1)
      return(x[idx])
    }
    
  })
  return(as.character(cls))
}
get_ids <- function(descs, idx) {
  labels <- sapply(strsplit(descs, "_"), function(x) {
    return(x[idx])
  })
  return(as.character(labels))
}
my_font <- "Helvetica"
thesis_theme <- ggplot2::theme(
  text = element_text(size = 8, hjust = 0.5, family=my_font),
  # axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
  # axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
  # axis.text.x = element_text(color="black",size=7, hjust = 0.5, family=my_font, angle = 90),
  # axis.text.x = element_blank(),
  # axis.ticks = element_blank(),
  strip.placement = "outside",
  # axis.line = element_line(colour = "black"),
  # axis.line = element_blank(),
  axis.text = element_text(size=9, face="bold",family=my_font, hjust = 0.5),#color="black",
  axis.title = element_text(size=9, face="bold",family=my_font, hjust = 0.5),
  plot.title = element_text(size=10, face="bold",family=my_font), #, hjust = 0.5
  legend.title=element_text(size=9, hjust = 0.5, family=my_font),
  # legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
  strip.text.x = element_text(size=9, family=my_font),
  strip.text.y = element_text(size=9, family=my_font),
  # legend.position = lg_pos,
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(-2,-2,-2,-2),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "grey50", fill=NA, size=0.5),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank()
)

viz_clonealign_correlation <- function(){
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
  clonealign_dir <- paste0(base_dir,'clonealign/')
  dlp_dir <- paste0(base_dir,'dlp/')
  fns <- list.files(clonealign_dir)    
  fns <- fns[grepl('*csv',fns)]
  
  clonealign_stat <- tibble::tibble()
  dlp_stat <- tibble::tibble()
  for(f in fns){
    if(!file.exists(paste0(dlp_dir,f))){
      stop('File does not exist: ')
      print(f)
    }
    
    # First clonealign output
    sdf <- data.table::fread(paste0(clonealign_dir,f)) %>% as.data.frame()
    # print(dim(sdf))
    lib <- gsub('.cache/','',sdf$Sample[1])
    lib <- gsub('/filtered_feature_bc_matrix','',lib)
    clones_props <- sdf %>%
      dplyr::filter(!clone %in% c('unassigned','Unassigned','None')) %>%
      dplyr::group_by(clone) %>%
      dplyr::summarise(clone_prevalence=n(), freq=n()/dim(sdf)[1]) %>%
      as.data.frame() %>%
      dplyr::mutate(sample_id=gsub('.csv','',f), library_id=lib)
    
    clonealign_stat <- dplyr::bind_rows(clonealign_stat, clones_props)
    
    # DLP prevalence
    dlp_df <- data.table::fread(paste0(dlp_dir,f)) %>% as.data.frame()
    # print(dim(dlp_df))
    dlp_df$sample_id <- f
    dlp_stat <- dplyr::bind_rows(dlp_stat, dlp_df)
    
  }
  
  colnames(dlp_stat)
  colnames(clonealign_stat)
  dim(clonealign_stat)
  clonealign_stat$clone_prevalence[3]
  clonealign_stat$clone
  dlp_stat <- dlp_stat %>%
    dplyr::rename(clone=cluster, clone_prevalence=n, timepoint=time)
  dlp_stat$freq[1]
  clonealign_stat1 <- clonealign_stat
  clonealign_stat1$clone <- get_unique_clone_id(clonealign_stat1$clone)
  
  unique(clonealign_stat1$clone)
  head(clonealign_stat1)
  clonealign_stat1$data_type <- 'scRNA-seq (10x)'  
  
  dlp_stat <- dlp_stat %>%
    dplyr::filter(clone %in% unique(clonealign_stat1$clone))%>%
    dplyr::select(-timepoint)
  dlp_stat$data_type <- 'scWGS (DLP)'
  View(head(summary_stat))
  summary_stat <- dplyr::bind_rows(dlp_stat, clonealign_stat1)
  dim(summary_stat)  
  summary_stat$PDX <- stringr::str_sub(summary_stat$sample_id,1,6)
  summary_stat$PDX <- gsub('X','',summary_stat$PDX)
  unique(summary_stat$sample_id)
  # table(summary_stat$sample_id, summary_stat$data_type)
  summary_stat$sample_id <- gsub('.csv','',summary_stat$sample_id)
  summary_stat$sample_id <- gsub('XB0','_',summary_stat$sample_id)
  # sids <- unique(summary_stat$sample_id)
  # 
  # summary_stat$sample_id <- stringr::str_sub(summary_stat$sample_id,1,8)
  # summary_stat$sample_id <- gsub('X$','',summary_stat$sample_id)
  summary_stat$sample_id <- ifelse(summary_stat$data_type=='scRNA-seq (10x)', 
                                   paste0(summary_stat$sample_id, ': 10x'), paste0(summary_stat$sample_id, ': DLP'))
  pdx1 <- c('SA501','SA530','SA604','SA609')
  pdx2 <- c('SA535','SA1035')
  summary_stat1 <- summary_stat %>%
    dplyr::filter(PDX %in% pdx1)
  summary_stat1$PDX <- factor(summary_stat1$PDX, levels = pdx1)
  summary_stat2 <- summary_stat %>%
    dplyr::filter(PDX %in% pdx2)
  summary_stat2$PDX <- factor(summary_stat2$PDX, levels = pdx2)
  p1 <- viz_correlation(summary_stat1)
  p2 <- viz_correlation(summary_stat2)
  plg <- viz_legend(summary_stat1) 
  pmain <- cowplot::plot_grid(p1, p2, ncol = 2, rel_heights = c(1,1))
  ptotal <- cowplot::plot_grid(pmain, plg, ncol = 1, rel_heights = c(16,1))
  save_dir <- base_dir
  png(paste0(save_dir,"Fig2_part2_clonealign_correlation.png"), height = 2*600, width=2*650, res = 2*72)
  print(ptotal)
  dev.off()
  

  summary_stat <- dplyr::bind_rows(dlp_stat, clonealign_stat1)
  data.table::fwrite(summary_stat, paste0(save_dir, 'summary_stat_10x_dlp.csv'))
  
  
}

viz_correlation_SA609 <- function(){
  
  meta_data <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/metadata_drug_resistance/SA609_10x.csv') %>% as.data.frame()
  dim(meta_data)  
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
  save_dir <- base_dir
  # summary_stat <- data.table::fread(paste0(save_dir, 'summary_stat_10x_dlp.csv'))
  summary_stat$sample_id[1]
  unique(summary_stat$PDX)
  summary_stat <- summary_stat %>%
    dplyr::filter(PDX=='SA609')
  dim(summary_stat)
  unique(summary_stat$PDX)
  summary_stat$data_type <- gsub(' ','_',summary_stat$data_type)
  summary_stat$sample_id <- gsub('.csv','',summary_stat$sample_id)
  summary_stat$PDX <- stringr::str_sub(summary_stat$sample_id,1,6)
  summary_stat$PDX <- gsub('X','',summary_stat$PDX)
  
  
  summary_stat$sample_id <- gsub('.csv','',summary_stat$sample_id)
  # summary_stat$sample_id <- gsub('XB0','_',summary_stat$sample_id)
  datatag=='SA609'
  if(datatag=='SA609'){
    obs_samples <- c('SA609X3XB01584','SA609X4XB03080','SA609X5XB03223',
                     'SA609X6XB03447','SA609X7XB03554','SA609X4XB003083',
                     'SA609X5XB03230','SA609X6XB03404','SA609X7XB03505',
                     'SA609X5XB03231','SA609X6XB03401','SA609X7XB03510')
    # sids <- unique(cell_clones$mouse_id)
    # t <- sids[!sids %in% obs_samples]
    # x <- grep('SA609X3*',t,value=T)

    summary_stat <- summary_stat %>%
      dplyr::filter(sample_id %in% obs_samples)
    print('Replicates excluded!!!')
    length(unique(summary_stat$sample_id))
    length(obs_samples)
  }
  dim(summary_stat)
  
  sum(unique(summary_stat$sample_id) %in% meta_data$mouse_id)
  colnames(meta_data)
  meta_data <- meta_data %>%
    dplyr::select(passage, mouse_id, treatmentSt)
  
  summary_stat <- summary_stat %>% left_join(meta_data, by=c('sample_id'='mouse_id'))
  ts <- unique(summary_stat$treatmentSt)
  df <- tibble::tibble()
  for(obs_treatment in c('UnRx','Rx','RxH')){
    if(obs_treatment=='Rx'){
      conds <- grep('*T$', ts, value=T)
    } else if(obs_treatment=='RxH'){
      conds <- grep('*TU$', ts, value=T)
    }else{
      conds <- c('U',grep('*UU$', ts, value=T))
    }
    tmp <- data.frame(treatmentSt=conds, treatment_desc=rep(obs_treatment, length(conds)))
    df <- dplyr::bind_rows(df, tmp)
  }
  
  summary_stat$data_type1 <- ifelse(summary_stat$data_type=='scRNA-seq (10x)', '10x', 'DLP')
  summary_stat$data_type <- ifelse(summary_stat$data_type=='10x','scRNA-seq (10x)', 'scWGS (DLP)')
  unique(summary_stat$data_type)
  dim(summary_stat)
  summary_stat <- summary_stat %>% left_join(df, by=c('treatmentSt'))
  summary_stat$desc <- paste0('Pt4:',summary_stat$passage,',',summary_stat$treatment_desc, ' (',summary_stat$data_type1,')')
  length(unique(summary_stat$desc))
  length(unique(summary_stat$sample_id))
  View(meta_data)
  summary_stat$p
  sids <- unique(summary_stat$sample_id)
  p1 <- viz_correlation(summary_stat, xstring = 'clone', ystring = 'desc')
  res <- viz_legend(summary_stat) 
  
  p1 <- p1 + labs(title='Pt4 - clonal proportions')
  pSA609_clonealign <- cowplot::plot_grid(p1, res$colplt,res$sizeplt, ncol = 1, rel_heights = c(10,0.9,0.9))
  save_dir <- base_dir
  png(paste0(save_dir,"Fig2_part22_clonealign_correlation.png"), height = 2*400, width=2*450, res = 2*72)
  print(ptotal)
  dev.off()
  saveRDS(p1, paste0(save_dir,"clonealign_correlation-SA609.rds"))
  p1 <- readRDS(paste0(save_dir,"clonealign_correlation-SA609.rds"))
  saveRDS(res, paste0(save_dir,"clonealign_correlation-SA609_legends.rds"))
  res <- readRDS(paste0(save_dir,"clonealign_correlation-SA609_legends.rds"))
}


viz_total_correlation <- function(){
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
  save_dir <- base_dir
  summary_stat <- data.table::fread(paste0(save_dir, 'summary_stat_10x_dlp.csv')) %>% as.data.frame()
  summary_stat$data_type1 <- ifelse(summary_stat$data_type=='scRNA-seq (10x)', '10x', 'DLP')
  summary_stat$sample_id <- gsub('.csv','',summary_stat$sample_id)
  summary_stat$PDX <- stringr::str_sub(summary_stat$sample_id,1,6)
  summary_stat$PDX <- gsub('X','',summary_stat$PDX)
  unique(summary_stat$PDX)
  
  #Replicates excluded!!!
  obs_samples <- c('SA609X3XB01584','SA609X4XB03080','SA609X5XB03223',
                   'SA609X6XB03447','SA609X7XB03554','SA609X4XB003083',
                   'SA609X5XB03230','SA609X6XB03404','SA609X7XB03505',
                   'SA609X5XB03231','SA609X6XB03401','SA609X7XB03510')  
  
  summary_stat_SA609 <- summary_stat %>%
    dplyr::filter(PDX=='SA609' & sample_id %in% obs_samples)
  dim(summary_stat_SA609)
  
  length(unique(summary_stat_SA609$sample_id))
  length(obs_samples)
  summary_stat_others <- summary_stat %>%
    dplyr::filter(PDX!='SA609')
  dim(summary_stat)
  summary_stat <- dplyr::bind_rows(summary_stat_SA609, summary_stat_others)
  dim(summary_stat)
  data.table::fwrite(summary_stat, paste0(save_dir, 'summary_stat_10x_dlp_subsetSA609.csv'))
  summary_stat <- data.table::fread(paste0(save_dir, 'summary_stat_10x_dlp_subsetSA609.csv')) %>% as.data.frame()
  
  # Remove X9 SA535
  summary_stat <- summary_stat %>%
    dplyr::filter(!sample_id %in% c('SA535X9XB03617', 'SA535X9XB03616'))
  # test_dir<- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/metadata_drug_resistance/'
  # metadata <- data.table::fread(paste0(test_dir, 'SA535_10x.csv')) %>% as.data.frame()
  # metadata <- metadata %>%
  #   filter(PDX=='SA535_CY' & passage=='X9')
  summary_stat$desc <- paste0(summary_stat$PDX,'_',summary_stat$sample_id,'_',summary_stat$clone)
  stat <- summary_stat %>% 
    dplyr::mutate(data_type1=replace(data_type1,data_type1=='10x','clonealign'))%>% 
    tidyr::pivot_wider(id_cols = 'desc',names_from = 'data_type1', values_from = 'freq')%>% 
    na.omit()
  # %>% tibble::column_to_rownames('desc')#
  dim(stat)
  stat$desc[1]
  # View(head(stat))
  stat$PDX <- get_ids(stat$desc, 1)
  unique(stat$PDX)
  pt_idx <- data.frame(PDX=c("SA501","SA530","SA604","SA609","SA535","SA1035"), 
                       pt=paste0("Pt",rep(1:6,1)))
  stat <- stat %>% inner_join(pt_idx, by='PDX')
  
  corr_stat <- tibble::tibble()
  for(pd in unique(summary_stat$PDX)){
    tmp <- stat %>%
      dplyr::filter(PDX==pd)
    cr <- round(cor(tmp$DLP,tmp$clonealign, method = "pearson"),2)
    corr_stat <- dplyr::bind_rows(corr_stat, tibble::tibble(PDX=pd, correlation=cr))
  }  
  corr_stat <- as.data.frame(corr_stat)
  dim(stat)
  stat <- stat %>% inner_join(corr_stat, by='PDX')
  # stat$PDX <- paste0(stat$PDX,' (r=',stat$correlation,')')
  stat$PDX <- paste0(stat$pt,' (r=',stat$correlation,')')
  unique(stat$PDX)  
  data.table::fwrite(stat, paste0(save_dir, '10x_dlp_stat_correlation.csv'))
  
  
  stat$Pt <- stat$PDX
  my_font <- "Helvetica"
  p <- ggplot(stat, aes(x = DLP, y = clonealign)) + 
    geom_point(aes(color = Pt, shape=Pt), alpha=0.9, size=3)  #, size=2*log10(pct_genes)size=4.5
  
  p <- p + thesis_theme
  lg_pos <- "bottom"
  p <- p + ggplot2::theme(legend.position = lg_pos) +
    labs(x='Clonal proportion DLP', y='Clonal proportion clonealign 10x', 
         title = '6 patients - Clonal proportions')
  
  p <- p + guides(color = guide_legend(nrow = 3, override.aes = list(size=3.5)))
  saveRDS(p, paste0(save_dir,'clonealign_correlation_all_patients.rds'))
  p_corr <- cowplot::plot_grid(p, pSA609_clonealign, ncol = 2, labels = c('b','c'))
  
  
  # p_corr
  
  png(paste0(save_dir,"Fig2_part2_clonealign_correlation_v2.png"), height = 2*420, width=2*650, res = 2*72)
  print(p_corr)
  dev.off()
  saveRDS(p_corr, paste0(save_dir,"clonealign_correlation-all_v2.rds"))
  
  data.table::fwrite(corr_stat, paste0(save_dir, 'correlation_10x_dlp.csv'))
}  


viz_correlation <- function(summary_stat, xstring='clone', ystring='sample_id'){
  colorcode <- c('darkgreen','#DA70D6')
  names(colorcode) <- c('scRNA-seq (10x)','scWGS (DLP)')
  my_font <- "Helvetica"
  
  p <- ggplot(summary_stat, aes_string(x = xstring, y = ystring)) + 
    geom_point(aes(color = data_type, size=freq), alpha=0.9) + #, size=2*log10(pct_genes)size=4.5
    scale_color_manual(values = colorcode) + 
    # facet_grid(PDX ~ ., scales="free_y", space='free',drop=T) + # PDX ~ ref_gene, . ~ PDX
    # geom_hline(yintercept = 1:6, col = "#e2e2e2") +
    # geom_text(aes(label=celltype_desc))+
    # annotate('text', x = df$success, y = df$index, label = df$success, size=3, colour = col[df$gender])+
    
    # scale_color_manual(values = col) +
    theme_bw() + 
    theme(strip.text = element_text(size=9, color="black", family=my_font),
          strip.background = element_blank(),
          legend.position = "none",
          # panel.grid.major.x = element_blank(),
          # panel.grid.minor.x = element_blank(),
          # axis.ticks.y = element_blank(),
          # axis.text  = element_blank(),
          axis.title.y = element_blank(),
          text = element_text(size = 7, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=9, hjust = 0.5, family=my_font),  #, angle = 90
          axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=10, face="bold", hjust=0, family=my_font)) #+
  
  
 
  return(p)
}

viz_legend <- function(summary_stat){
  colorcode <- c('darkgreen','#DA70D6')
  names(colorcode) <- c('scRNA-seq (10x)','scWGS (DLP)')
  my_font <- "Helvetica"
  summary_stat$Clone_Prevalence <- summary_stat$freq
  p1 <- ggplot(summary_stat, aes(x = clone, y = sample_id)) + 
    geom_point(aes(color = data_type), alpha=0.9) + #, size=2*log10(pct_genes)size=4.5
    scale_color_manual(values = colorcode, name=' ') +
    theme(legend.position ='bottom')
  p1 <- p1 + guides(color = guide_legend(override.aes = list(size=7, nrow = 1))) #shape = 0, 
  lg <- cowplot::get_legend(p1)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  
  p1 <- ggplot(summary_stat, aes(x = clone, y = sample_id)) + 
    geom_point(aes(size=Clone_Prevalence), alpha=0.9) + #, size=2*log10(pct_genes)size=4.5
    scale_color_manual(values = colorcode, name=' ') +
    theme(legend.position ='bottom')
  p1 <- p1 + guides(size = guide_legend(override.aes = list(nrow = 1))) #shape = 0, 
  lg <- cowplot::get_legend(p1)
  plg2 <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  return(list(colplt=plg,sizeplt=plg2))
}  
