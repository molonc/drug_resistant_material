# library(remotes)
# remotes::install_version("Rttf2pt1", version = "1.3.8")

suppressPackageStartupMessages({
  library(gtools)
  library(parallel)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggrepel)
})

library(extrafont)
font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
fonts()


plot_proportion_cistrans_reference_set <- function(pair_groups, 
                                                   reference_sets=c('CoreFitness','BroadSanger','cosmic','cisplatin_resistance'))
{
  minLogFC <- 0.5
  FDR_cutoff <- 0.01
  pValueThrs <- 0.05
  ref_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/cancer_reference_genes/'
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
  pair_groups <- data.table::fread(paste0(input_dir,'comparisons_drug_res.csv')) %>% as.data.frame()
  rownames(pair_groups) <- pair_groups$desc
  # View(pair_groups)
  dim(pair_groups)
  ref_stat <- tibble::tibble()
  for(s in reference_sets){
    print('Loading reference genes set: ')
    print(s)
    ref_df <- read.csv(paste0(ref_dir, s,'.csv'), stringsAsFactors=F, check.names = F)
    print(dim(ref_df))
    
    for(de in pair_groups$desc){
      signif_genes <- read.csv(paste0(base_dir,pair_groups[de,'datatag'],'/',pair_groups[de,'desc'],'/','signif_genes_FDR_0.01.csv'), check.names = F, stringsAsFactors=F)
      signif_genes <- signif_genes[abs(signif_genes$logFC)>minLogFC & signif_genes$FDR<FDR_cutoff  & signif_genes$PValue<pValueThrs,]
      signif_genes <- signif_genes %>%
        dplyr::filter(gene_symbol %in% ref_df$gene_symb)
      print(dim(signif_genes))
      # fitness_genes <- summary(as.factor(signif_genes$classified_gene_dlp))
      if(dim(signif_genes)[1]>0){
        nb_incis <- signif_genes %>%
          dplyr::filter(!is.na(classified_gene_dlp))%>%
          dplyr::summarise(n_cell=n())%>%
          dplyr::pull(n_cell)
        nb_intrans <- signif_genes %>%
          dplyr::filter(is.na(classified_gene_dlp))%>%
          dplyr::summarise(n_cell=n())%>%
          dplyr::pull(n_cell)
        # nb_incis <- sum(!is.na(signif_genes$classified_gene_dlp)==TRUE)
        # nb_intrans <- sum(is.na(signif_genes$classified_gene_dlp)==TRUE)
        # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
        tmp <- tibble(desc=de,ref_set=s,PDX=pair_groups[de,'datatag'],
                      incis=round(nb_incis*100/nrow(ref_df),2),
                      intrans=round(nb_intrans*100/nrow(ref_df),2))
        ref_stat <- dplyr::bind_rows(ref_stat, tmp)
      }
    }
    
    
    ref_stat <- as.data.frame(ref_stat)
    plt_desc_df <- data.table::fread(paste0(input_dir,'plt_desc_v2.csv')) %>% as.data.frame()
    dim(plt_desc_df)
    dim(pair_groups)
    plt_desc_df$ord <- rep(1:dim(plt_desc_df)[1],1)
    ref_stat <- ref_stat %>% left_join(plt_desc_df, by=c('desc'))
    dim(ref_stat)
    
    
    data.table::fwrite(ref_stat, paste0(base_dir,'summary_references_genes.csv'))
    
    
    
    save_dir <- paste0(base_dir,'reference_stat/')
    ref_stat <- data.table::fread(paste0(save_dir,'summary_references_genes.csv'))
    dim(ref_stat)
    output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/Fig4_prevalence_cistrans/'
    stat1 <- data.table::fread(paste0(output_dir,'summary_incis_intrans_genes_plots.csv')) %>% as.data.frame()
    ref_stat$desc[1]
    stat1$desc[1]
    pdxs = c('SA609','SA535','SA1035')
    prop_plt_ls <- list()
    for(pd in pdxs){
      df <- stat1 %>% 
        dplyr::filter(PDX==pd & comp_type=='treated_vs_untreated')
      dim(df)
      
      descs <- unique(df$desc)
      ref_stat_tmp <- ref_stat %>% 
        dplyr::filter(desc %in% descs)
      
      prop_plt <- plot_reference_set_proportion(ref_stat_tmp, save_dir)
      prop_plt_ls[[pd]] <- prop_plt$prevalence_plt
    }  
    prop_plt_ls$SA609
    res_prevalence_plt <- plot_reference_set_proportion(ref_stat, save_dir)
    res_prevalence_plt$prevalence_plt
    pathway_stat <- get_pathway_stat(pair_groups, input_dir, save_dir)
    
    pathway_stat <- data.table::fread(paste0(save_dir,'stat_fgsea_references_sets.csv')) %>% as.data.frame()
    pathway_stat <- pathway_stat %>% left_join(plt_desc_df, by=c('desc'))
    dim(pathway_stat)
    unique(pathway_stat$desc)
    res_pathway <- viz_reference_pathways(pathway_stat)
    
    pathway_plt_ls <- list()
    for(pd in pdxs){
      df <- stat1 %>% 
        dplyr::filter(PDX==pd & comp_type=='treated_vs_untreated')
      dim(df)
      
      descs <- unique(df$desc)
      pathway_stat_tmp <- pathway_stat %>% 
        dplyr::filter(desc %in% descs)
      # ref_stat_tmp$desc <- factor(ref_stat_tmp$desc, levels = descs)
      dim(pathway_stat_tmp)
      pathway_plt <- viz_reference_pathways(pathway_stat_tmp)
      pathway_plt_ls[[pd]] <- pathway_plt$pathway_sig_plt
    }  
    pathway_plt_ls$SA609
    pathway_plt$plg
    res_pathway$pathway_sig_plt
    res_pathway$pathway_sig_plt
    res_prevalence_plt <- readRDS(paste0(save_dir,'proportion_ref_genes.rds'))
    
    main_plt <- cowplot::plot_grid(res_prevalence_plt$prevalence_plt, res_pathway$pathway_sig_plt,
                                   rel_widths = c(6,1), ncol=2, align = 'h')
    lg_plt <- cowplot::plot_grid(res_prevalence_plt$plg, res_pathway$plg,
                                   rel_widths = c(6,1), ncol=2, align = 'h')
    
    ptotal <- cowplot::plot_grid(main_plt, lg_plt,
                                 ncol = 1, rel_heights = c(15,1), align = 'v')
    png(paste0(save_dir,"Fig2_part2.png"), height = 2*500, width=2*650, res = 2*72)
    print(ptotal)
    dev.off()
    
  }
}  

viz_reference_pathways <- function(pathway_stat){
  # sum(pathway_stat$padj<=0.05, na.rm = T)
  my_font <- "Helvetica"
  ref_genes=c('BS','CR','CF','CM')
  tmp <- data.frame(pathway=unique(pathway_stat$pathway), ref_gene=ref_genes)
  tmp <- tmp %>%
    dplyr::filter(ref_gene %in% c('CF','CR'))
  pathway_stat <- pathway_stat %>% inner_join(tmp, by=c('pathway'))
  
  # pathway_stat$ref_gene <- ifelse(pathway_stat$ref_gene=='Cisplt Resistance','Cisplt Res',
  #                                 pathway_stat$ref_gene)
  pathway_stat <- pathway_stat %>%
    dplyr::mutate(padj=replace(padj, padj>=0.05, 0))
  library(viridis)
  # pathway_stat$cell_size <- ifelse(pathway_stat$padj>0, log10(pathway_stat$size), NA)
  pathway_stat$cell_size <- ifelse(pathway_stat$padj>0, 5, NA)
  
  # pdx_level = c('SA501','SA530','SA604','SA609','SA535','SA1035')
  # pathway_stat$PDX <- factor(pathway_stat$PDX, levels = pdx_level)
  pathway_stat$ref_gene <- factor(pathway_stat$ref_gene, levels = c('CF','CR'))
  
  pathway_stat$ord <- factor(pathway_stat$ord)
  pathway_stat$type <- 'Stat Level'
  pathway_stat$plt_desc <- factor(pathway_stat$plt_desc, levels=unique(pathway_stat$plt_desc))
  p <- ggplot(pathway_stat, aes(x=ref_gene, y=plt_desc, color=padj)) +
    geom_point(aes(size=cell_size)) +
    viridis::scale_color_viridis(discrete=F, alpha=0.8)  + 
    # facet_grid(PDX ~ type, scales="free_y", space='free',drop=T) + # PDX ~ ref_gene, . ~ PDX
    theme_bw() + 
    theme(strip.text = element_text(size=9, color="black", family=my_font),
          strip.background = element_blank(),
          legend.position = "none",
          # panel.grid.major.x = element_blank(),
          # panel.grid.minor.x = element_blank(),
          # axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font, angle = 90),
          text = element_text(size = 7, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=9, hjust = 0.5, family=my_font, angle = 90),  #, angle = 90
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=10, hjust=0.5, family=my_font),
          axis.title = element_blank()) 
  p <- p + labs(title='Signif Level')
  # p <- p + guides(color=guide_legend(title="Sigf padj")) #size = FALSE,
  # p
  p1 <- ggplot(pathway_stat, aes(x=ref_gene, y=plt_desc, color=padj)) +
    geom_point(size=5) +
    viridis::scale_color_viridis(discrete=F, alpha=0.8, name = 'P-adjusted  ') +
    theme(legend.position ='bottom',
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.8, 'cm'), 
          legend.text = element_text(size = 7, hjust = 0.5, family=my_font))
  # p1 <- p1 + guides(fill = guide_legend(override.aes = list(shape = 0, size=8, nrow = 1)))
  
  lg <- cowplot::get_legend(p1)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  # plg
  res_pathway <- list(pathway_sig_plt=p, plg=plg)
  # saveRDS(res_pathway, paste0(save_dir,'pathway_ref_plt.rds'))
  # ptotal <- cowplot::plot_grid(p, plg, ncol = 1, rel_heights = c(15,1))
  # png(paste0(save_dir,"Fig2_part22.png"), height = 2*600, width=2*150, res = 2*72)
  # print(ptotal)
  # dev.off()
  
  return(res_pathway)
  
}
    
get_pathway_stat <- function(pair_groups, input_dir, save_dir){
  # source('/home/htran/Projects/farhia_project/rnaseq/pipeline/utils/pathway_utils.R')
  ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'  
  ref_genes_ls <- load_reference_genes_set(ref_dif)
  ref_genes_ls <- ref_genes_ls[!names(ref_genes_ls) %in% c('METASTASIS')]
  print(names(ref_genes_ls))
  minLogFC <- 0.25  # for pathway test, using this threshold
  FDR_cutoff <- 0.05
  pValueThrs <- 0.05
  
  pathway_stat <- tibble::tibble()
  for(de in pair_groups$desc){
    signif_genes <- read.csv(paste0(input_dir,pair_groups[de,'result_fn']), check.names = F, stringsAsFactors=F)
    signif_genes <- signif_genes[abs(signif_genes$logFC)>minLogFC & signif_genes$FDR<FDR_cutoff  & signif_genes$PValue<pValueThrs,]
    print(dim(signif_genes))
    if(dim(signif_genes)[1]>0){
      # Deal with overlapping and NA gene symbols
      signif_genes <- signif_genes %>%
        dplyr::select(gene_symbol, logFC) %>%
        na.omit() %>%
        distinct() %>%
        group_by(gene_symbol) %>%
        summarize(logFC=mean(logFC))
      
      deg_stat <- signif_genes$logFC
      names(deg_stat) <- signif_genes$gene_symbol
      gsea_out <- fgsea(pathways=ref_genes_ls, stats=deg_stat, nPermSimple = 10000) #, nperm=10000
      class(gsea_out)
      gsea_out$important_genes <- ''
      for(i in rep(1:length(gsea_out$leadingEdge),1)){
        gsea_out$important_genes[i] <- as.character(paste(unlist(gsea_out$leadingEdge[i]), collapse = ','))
      }
      colnames(gsea_out)
      gsea_out <- gsea_out %>%
        dplyr::select(-leadingEdge)%>%
        as.data.frame()
      
      gsea_out$desc=de
      gsea_out$PDX=pair_groups[de,'datatag']
      pathway_stat <- dplyr::bind_rows(pathway_stat, gsea_out)
      
    }
  }
  pathway_stat <- as.data.frame(pathway_stat)
  print(dim(pathway_stat))
  # View(head(pathway_stat[1:3,1:3]))
  # View(gsea_out)
  # sum(pathway_stat$padj<0.05, na.rm = T)
  data.table::fwrite(pathway_stat, paste0(save_dir,'stat_fgsea_references_sets.csv'))
  return(pathway_stat)
}




plot_reference_set_proportion <- function(stat, save_dir){
  reference_sets=c('CoreFitness','BroadSanger','cosmic','cisplatin_resistance')
  ref_genes=c('Core Fitness (CF)','BS Essential','COSMIC','Cispl Resistance (CR)')
  tmp <- data.frame(ref_set=reference_sets, ref_gene=ref_genes)
  stat <- stat %>% pivot_longer(cols = starts_with("in"),
                                names_to = 'Gene_Type', values_to = 'pct_gene')
  dim(stat)
  stat <- stat %>% left_join(tmp, by=c('ref_set'))
  stat <- stat %>% 
    dplyr::filter(ref_set %in% c('CoreFitness','cisplatin_resistance'))
  
  my_font <- "Helvetica"
  print(unique(stat$PDX))
  colorcode <- c('#FF3232','#40A0E0')
  names(colorcode) <- c('incis','intrans')
  
  stat$proportion_size <- ifelse(stat$pct_gene==0, 0, log10(stat$pct_gene))
  # pdx_level = c('SA501','SA530','SA604','SA609','SA535','SA1035')
  # stat$PDX <- factor(stat$PDX, levels = pdx_level)
  # lvs <- c('Core Fitness (CF)','Cispl Resistance (CR)')
  # stat$ref_gene <- factor(stat$ref_gene, levels = lvs)
  # stat$plt_desc <- paste0(stat$PDX,': ',stat$plt_desc)
  min(stat$pct_gene)
  stat$pct_gene <- ifelse(stat$pct_gene==0, 0.01, stat$pct_gene)
  stat$plt_desc <- factor(stat$plt_desc, levels = unique(stat$plt_desc))
  p <- ggplot(stat, aes(x = pct_gene, y = plt_desc)) + 
    geom_point(aes(color = Gene_Type, size=proportion_size), alpha=0.9) + #, size=2*log10(pct_genes)size=4.5
    scale_color_manual(values = colorcode) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    # scale_x_continuous(trans='log10') + 
    facet_grid(~ ref_gene, scales="free_y", space='free',drop=T) + # PDX ~ ref_gene, . ~ PDX
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
          axis.text.y  = element_blank(),
          text = element_text(size = 7, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=9, hjust = 0.5, family=my_font),  #, angle = 90
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=10, face="bold", hjust=0, family=my_font),
          axis.title = element_blank()) #+
    # geom_text_repel(data = stat, aes(label = pct_gene), nudge_x=0, nudge_y=0, max.overlaps = Inf,
    #                  min.segment.length = 0, segment.size=0.2,
    #                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
  # p
  
  p <- p + guides(color=guide_legend(title="Gene Type ", override.aes = list(shape = 0, size=6))) #size = FALSE,
  
  
  
  p1 <- ggplot(stat, aes(x = pct_gene, y = plt_desc)) + 
    geom_point(aes(color = Gene_Type), size=6.5, alpha=1) + #, size=2*log10(pct_genes)size=4.5
    scale_color_manual(name = 'Gene Type', values = colorcode, labels = c("In Cis", "In Trans")) + 
    theme(legend.position ='bottom')
  p1 <- p1 + guides(fill = guide_legend(override.aes = list(shape = 0, size=8, nrow = 1)))
  lg <- cowplot::get_legend(p1)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  # plg
  res <- list(prevalence_plt=p, plg=plg)
  # saveRDS(res, paste0(save_dir,'proportion_ref_genes.rds'))
  # ptotal <- cowplot::plot_grid(p, plg, ncol = 1, rel_heights = c(15,1.3))
  # png(paste0(save_dir,"Fig2_part21.png"), height = 2*600, width=2*520, res = 2*72)
  # print(ptotal)
  # dev.off()
  
  
  return(res)
}





# View(head(bs_df))
# stat1 <- dplyr::bind_rows(bs_df, fitness_df)
# dim(stat1)
# unique(stat1$PDX)
# stat1 <- stat1 %>%
#   filter(PDX!='SA535_CX5461')
# colnames(stat1)
# stat1$PDX <- gsub('_CISPLATIN','',stat1$PDX)
# head(stat1)
# stat1$Gene_Type <- ifelse(grepl('incis',stat1$desc),'in-cis','in-trans')
# stat1$desc <- gsub('_incis','',stat1$desc)
# stat1$desc <- gsub('_intrans','',stat1$desc)
# colorcode <- c('#FF3232','#40A0E0')
# names(colorcode) <- c('in-cis','in-trans')
# stat1$desc <- as.factor(stat1$desc)
# stat1$PDX <- as.factor(stat1$PDX)
# stat1$ref_gene <- as.factor(stat1$ref_gene)
# unique(stat1$desc)
# rv_desc <- c('SA609_UUUU_B_UUUU_C','SA609_UUUU_B_UUUU_H',
#              'SA1035_UTTTT_G_UUUUU_E','SA1035_UTTT_G_UUUU_E',
#              'SA535_UUTTT_R_UUUUU_J','SA535_UUTTT_T_UUUUU_Q','SA535_UUTTT_R_UUUUU_Q')
# stat1 <- stat1 %>%
#   filter(desc %in% rv_desc)
# pdxs = c('SA609','SA535','SA1035')
# stat1$PDX <- factor(stat1$PDX, levels = pdxs)
# stat1_backup <- stat1
# # stat1$desc <- paste0()
# unique(stat1$ref_gene)
# stat1$ref_gene <- ifelse(stat1$ref_gene=='DepMap BS','COSMIC','Cisplatin Resistance')
# stat1 <- stat1_backup
# p <- ggplot(stat1, aes(x = pct_gene, y = desc)) + 
#   geom_point(aes(color = Gene_Type, size=2*log10(pct_gene)), alpha=0.9) + #, size=2*log10(pct_genes)size=4.5
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
#   scale_color_manual(values = colorcode) + 
#   facet_grid(PDX ~ ref_gene, scales="free", space='free',drop=T) + # PDX ~ ref_gene, . ~ PDX
#   # geom_hline(yintercept = 1:6, col = "#e2e2e2") +
#   # geom_text(aes(label=celltype_desc))+
#   # annotate('text', x = df$success, y = df$index, label = df$success, size=3, colour = col[df$gender])+
#   
#   # scale_color_manual(values = col) +
#   theme_bw(base_size = 10) + 
#   theme(legend.position = "none",
#         # panel.grid.major.x = element_blank(),
#         # panel.grid.minor.x = element_blank(),
#         # axis.ticks.y = element_blank(),
#         # axis.text  = element_blank(),
#         strip.background =element_rect(fill="white", colour = NA), # for facet
#         strip.text = element_text(colour = 'black', size=11, face = "bold",family=my_font),  # for facet
#         text = element_text(size = 8, hjust = 0.5, family=my_font),
#         axis.text.x = element_text(size=10, hjust = 0.5, family=my_font),  #, angle = 90
#         axis.text.y = element_text(size=10, hjust = 0.5, family=my_font),
#         plot.title = element_text(size=10, face="bold", hjust=0, family=my_font),
#         axis.title = element_blank()) #+
#   # geom_text_repel( data = stat1, aes(label = pct_gene), nudge_x=0, nudge_y=0, max.overlaps = Inf,
#   #                  min.segment.length = 0, segment.size=0.2,
#   #                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
# p
# 
# p1 <- p1 + guides(color=guide_legend(title="Gene Type ", override.aes = list(size=6)))#size = FALSE,
# lg <- cowplot::get_legend(p1)
# plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
# 
# # plot_ls[['lg']] <- plg
# p_total <- cowplot::plot_grid(p, plg, nrow = 2, rel_heights = c(10,1))
# # p_total <- p_total + labs(x=NULL, y="(%) Gene Type ", title='Differentially Expressed Genes - Copy Number Driven Proportion')
# png(paste0(save_dir,"bs_fitness.png"), height = 2*300, width=2*650, res = 2*72)
# print(p_total)
# dev.off()
# 
# png(paste0(save_dir,"cosmic_cis.png"), height = 2*300, width=2*650, res = 2*72)
# print(p_total)
# dev.off()
# p_total_bs_fitness <- p_total
# p_total_cosmic_cis <- p_total
# save_dir1 <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/cis_trans/'
# p_cistran <- readRDS(paste0(save_dir1,"incis_intrans_genes_prevalence_dotplot.rds"))
# lbs <- c('Incis-intrans proportion','Broad-Sanger, Pancancer fitness reference','COSMIC, cisplatin reference genes')
# pfig2 <- cowplot::plot_grid(plotlist = list(p_cistran, p_total_bs_fitness, 
#                                             p_total_cosmic_cis), 
#                             ncol = 1, labels = lbs) # ,rel_heights = c(5,4,4)
# 
# pfig2
# pfig2 <- cowplot::plot_grid(p_cistran, p_total_bs_fitness, p_total_cosmic_cis)
# # p_total <- p_total + labs(x=NULL, y="(%) Gene Type ", title='Differentially Expressed Genes - Copy Number Driven Proportion')
# png(paste0(save_dir,"Fig2_mockup.png"), height = 2*13, width=2*650, res = 2*72)
# print(pfig2)
# dev.off()





# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/SA535-v6/mapped_genes/'
# 
# bs_df <- data.table::fread(paste0(save_dir,'broad_sanger_essential_gene_pct.csv')) %>% as.data.frame()
# fitness_df <- data.table::fread(paste0(save_dir,'ref_fitness_genes_pct.csv')) %>% as.data.frame()
# cosmic_df <- data.table::fread(paste0(save_dir,'reference_cosmic_gene_pct.csv')) %>% as.data.frame()
# cis_df <- data.table::fread(paste0(save_dir,'reference_cisplatin_resistance_gene_pct.csv')) %>% as.data.frame()
# 
# head(bs_df)
# bs_df <- bs_df %>%
#   rename(pct_gene=Broad_Sanger_gene)%>%
#   mutate(ref_gene='DepMap BS')
# # bs_df$
# head(fitness_df)
# fitness_df <- fitness_df %>%
#   rename(pct_gene=Fitness_genes_pct)%>%
#   mutate(ref_gene='Fitness PanCancer')
# 
# head(cosmic_df)
# cosmic_df <- cosmic_df %>%
#   rename(pct_gene=cosmic_genes_pct)%>%
#   mutate(ref_gene='COSMIC')
# 
# head(cis_df)
# cis_df <- cis_df %>%
#   rename(pct_gene=cisplatin_resistance_gene_pct)%>%
#   mutate(ref_gene='Cisplatin Resistance')
# 
# bs_df$desc <- paste0(bs_df$desc,'(*)')
# pbs <- plot_reference_set_proportion(bs_df)
# pfitness <- plot_reference_set_proportion(fitness_df)
# pcosmic <- plot_reference_set_proportion(cosmic_df)
# pcis <- plot_reference_set_proportion(cis_df)
# # sum(cosmic_df$pct_gene==0)
# pfig2 <- cowplot::plot_grid(pbs, pfitness, pcosmic, pcis)
# # p_total <- p_total + labs(x=NULL, y="(%) Gene Type ", title='Differentially Expressed Genes - Copy Number Driven Proportion')
# png(paste0(save_dir,"Fig2_mockup.png"), height = 2*500, width=2*650, res = 2*72)
# print(pfig2)
# dev.off()
